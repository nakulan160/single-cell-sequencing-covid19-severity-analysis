#!/usr/bin/env Rscript
# GSE196198 / Kim et al. (2023) — 01 metadata QC and input inventory
# Works with base R (Seurat is NOT needed for this stage).
# Run from the directory holding the 13 GSE196198 files:
#   Rscript 1_qc.R .
# Results: ./1_results/
# The GEO supplied metadata are already author-processed; this is an independent
# verification of reported QC criteria, NOT fresh QC on raw sequencing reads.

options(stringsAsFactors = FALSE, scipen = 99)
Sys.setenv(OMP_NUM_THREADS = "10", OPENBLAS_NUM_THREADS = "10", MKL_NUM_THREADS = "10")

run <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
  out_dir <- file.path(data_dir, "1_results")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  logfile <- file.path(out_dir, "1_qc_run_log.txt")
  notes <- character()
  say <- function(...) {
    msg <- paste0(..., collapse = "")
    cat(msg, "\n", sep = "")
    notes <<- c(notes, msg)
    writeLines(c(notes, "", capture.output(sessionInfo())), logfile)
  }

  say("GSE196198 paper reproduction — STEP 01 / metadata QC")
  say("Data directory: ", data_dir)
  say("R version: ", R.version.string)
  say("NOTE: counts and RNA values have NOT been renormalized or converted to counts.")

  collections <- c("pbmc", "T_cell", "B_cell", "NK_cell", "Monocyte")
  wanted <- c(
    paste0("GSE196198_", collections, "_meta.data.txt.gz"),
    paste0("GSE196198_", collections, "_normalized_count_integrated.txt.gz"),
    "GSE196198_TCR_seq_data.xlsx", "GSE196198_BCR_seq_data.xlsx",
    "GSE196198_Feature_barcode_CITE-seq.csv.gz"
  )
  locations <- file.path(data_dir, wanted)
  manifest <- data.frame(
    filename = wanted, exists = file.exists(locations),
    size_MB = round(as.numeric(file.info(locations)$size) / 1024^2, 2)
  )
  write.csv(manifest, file.path(out_dir, "file_inventory.csv"), row.names = FALSE)
  if (!all(manifest$exists)) {
    say("ERROR: missing input files: ", paste(manifest$filename[!manifest$exists], collapse = ", "))
    stop("Cannot continue: 13-file inventory incomplete.")
  }
  say("Found all ", length(wanted), " expected GEO supplementary files.")

  needed_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA",
                   "percent.mt", "seurat_clusters", "severity", "id.match", "severity2")

  read_metadata <- function(collection) {
    path <- file.path(data_dir, paste0("GSE196198_", collection, "_meta.data.txt.gz"))
    say("Reading metadata: ", basename(path))
    cn <- gzfile(path, open = "rt")
    on.exit(close(cn), add = TRUE)
    # GEO Seurat write.table() exports quote the column headings and IDs.
    # R must recognize those quotes rather than keeping literal " characters.
    dat <- read.table(cn, header = TRUE, row.names = 1, sep = "",
                      quote = "\"", comment.char = "", check.names = FALSE,
                      stringsAsFactors = FALSE)
    names(dat) <- trimws(sub('"$', '', sub('^"', '', names(dat))))
    missing <- setdiff(needed_cols, names(dat))
    if (length(missing)) {
      say("Parsed metadata columns: ", paste(names(dat), collapse = " | "))
      stop(collection, " missing metadata columns: ", paste(missing, collapse = ", "),
           ". Send the parsed-column line for review.")
    }
    if (anyDuplicated(rownames(dat))) stop("Duplicated cell row names: ", collection)
    if (anyNA(dat$id.match) || anyDuplicated(dat$id.match)) {
      stop("Missing/duplicated id.match within ", collection)
    }
    for (column in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
      dat[[column]] <- suppressWarnings(as.numeric(dat[[column]]))
      if (any(!is.finite(dat[[column]]))) stop("Invalid QC values in ", collection, ": ", column)
    }
    for (column in c("orig.ident", "severity", "severity2", "id.match")) {
      dat[[column]] <- as.character(dat[[column]])
      if (anyNA(dat[[column]]) || any(!nzchar(dat[[column]]))) {
        stop("Missing ", column, " in ", collection)
      }
    }
    # Exact inequalities reported in the publication supplement.
    dat$qc_pass_paper <- (
      dat$nFeature_RNA > 200 & dat$nFeature_RNA < 3500 &
      dat$nCount_RNA < 15000 & dat$percent.mt < 15
    )
    dat
  }

  lineage_rows <- list()
  sample_rows <- list()
  pbmc <- NULL
  for (collection in collections) {
    dat <- read_metadata(collection)
    by_sample <- split(seq_len(nrow(dat)), dat$orig.ident)
    sample_summary <- do.call(rbind, lapply(names(by_sample), function(sample) {
      idx <- by_sample[[sample]]
      vals <- dat[idx, , drop = FALSE]
      severity <- unique(vals$severity)
      severity2 <- unique(vals$severity2)
      if (length(severity) != 1 || length(severity2) != 1) {
        stop("Conflicting severity labels within sample ", sample)
      }
      data.frame(
        lineage = collection, sample = sample,
        severity = severity, severity2 = severity2,
        n_cells = length(idx), n_pass_paper_qc = sum(vals$qc_pass_paper),
        n_fail_paper_qc = sum(!vals$qc_pass_paper),
        pass_percent = round(100 * mean(vals$qc_pass_paper), 3),
        median_nCount_RNA = median(vals$nCount_RNA),
        median_nFeature_RNA = median(vals$nFeature_RNA),
        median_percent_mt = median(vals$percent.mt),
        stringsAsFactors = FALSE
      )
    }))
    rownames(sample_summary) <- NULL
    sample_rows[[collection]] <- sample_summary
    lineage_rows[[collection]] <- data.frame(
      lineage = collection, n_cells = nrow(dat), n_samples = length(by_sample),
      n_pass_paper_qc = sum(dat$qc_pass_paper),
      n_fail_paper_qc = sum(!dat$qc_pass_paper),
      min_nCount = min(dat$nCount_RNA), max_nCount = max(dat$nCount_RNA),
      min_nFeature = min(dat$nFeature_RNA), max_nFeature = max(dat$nFeature_RNA),
      min_percent_mt = min(dat$percent.mt), max_percent_mt = max(dat$percent.mt)
    )
    say(sprintf("%-9s cells=%d samples=%d QC-pass=%d QC-fail=%d",
                collection, nrow(dat), length(by_sample),
                sum(dat$qc_pass_paper), sum(!dat$qc_pass_paper)))
    if (collection == "pbmc") pbmc <- dat
    rm(dat)
    invisible(gc(verbose = FALSE))
  }
  lineage_summary <- do.call(rbind, lineage_rows)
  sample_summary <- do.call(rbind, sample_rows)
  write.csv(lineage_summary, file.path(out_dir, "qc_lineage_summary.csv"), row.names = FALSE)
  write.csv(sample_summary, file.path(out_dir, "qc_sample_summary.csv"), row.names = FALSE)

  # Keep the PBMC identifiers intact for later gene-expression matching.
  pbmc_out <- data.frame(cell_id = rownames(pbmc), pbmc,
                         check.names = FALSE, row.names = NULL)
  cn <- gzfile(file.path(out_dir, "pbmc_metadata_with_qc.csv.gz"), "wt")
  write.csv(pbmc_out, cn, row.names = FALSE)
  close(cn)
  rm(pbmc_out)
  invisible(gc(verbose = FALSE))

  # For longitudinal review, derive a DONOR CODE from the sample label only;
  # do NOT call the optional suffix actual weeks or clinical days.
  pbmc_samples <- subset(sample_summary, lineage == "pbmc")
  stripped <- sub("^Exp_", "", pbmc_samples$sample)
  donor <- sub("_.*$", "", stripped)
  pbmc_samples$donor_code_inferred <- donor
  pbmc_samples$visit_code_inferred <- ifelse(grepl("_", stripped), sub("^[^_]+_", "", stripped), "")
  write.csv(pbmc_samples, file.path(out_dir, "pbmc_sample_donor_review.csv"), row.names = FALSE)
  if (!all(grepl("^[HCMS][0-9]+$", donor))) {
    say("WARNING: some inferred donor labels are nonstandard. Review sample_donor_review.csv manually.")
  }

  # Three paper-QC distributions across PBMC sample labels.
  plot_box <- function(column, title, limit) {
    boxplot(split(pbmc[[column]], pbmc$orig.ident), outline = FALSE,
            las = 2, cex.axis = 0.66, col = "grey85", border = "black",
            ylab = column, main = title)
    abline(h = limit, col = "red3", lty = 2, lwd = 2)
  }
  png(file.path(out_dir, "qc_pbmc_by_sample.png"), width = 1900, height = 2100, res = 160)
  old <- par(mfrow = c(3, 1), mar = c(7, 5, 3, 1))
  plot_box("nFeature_RNA", "PBMC genes detected per cell (thresholds >200 and <3500)", 3500)
  abline(h = 200, col = "red3", lty = 2, lwd = 2)
  plot_box("nCount_RNA", "PBMC UMIs per cell (threshold <15000)", 15000)
  plot_box("percent.mt", "PBMC mitochondrial percentage (threshold <15%)", 15)
  par(old)
  dev.off()

  # A cell count is NOT an independent-patient count. This is descriptive only.
  size_table <- table(pbmc$severity)
  png(file.path(out_dir, "pbmc_cells_by_severity.png"), width = 1400, height = 950, res = 160)
  barplot(size_table, col = "grey75", border = "black", las = 2,
          ylab = "Number of PBMC cells (not patient replicates)",
          main = "Captured PBMC cells by supplied severity label")
  dev.off()

  # Make an independent donor-level review table (not a statistical test).
  donor_table <- unique(pbmc_samples[c("donor_code_inferred", "severity", "severity2")])
  write.csv(donor_table, file.path(out_dir, "donor_severity_review.csv"), row.names = FALSE)
  say("PBMC supplied severity labels: ", paste(names(size_table), size_table, sep = "=", collapse = "; "))
  say("Inferred unique donor labels: ", length(unique(pbmc_samples$donor_code_inferred)))
  say("WARNING: These files already passed much of the authors' filtering;")
  say("         the QC check may retain every exported cell. That is expected.")
  say("IMPORTANT: Metadata cannot reconstruct raw per-gene UMI counts.")
  say("IMPORTANT: No Seurat normalization, RPCA, clustering, UMAP, or DE was run.")
  say("SUCCESS. Check 1_results/qc_lineage_summary.csv, qc_sample_summary.csv,")
  say("         pbmc_sample_donor_review.csv and qc_pbmc_by_sample.png")
}

tryCatch(run(), error = function(e) {
  cat("STEP 01 FAILED: ", conditionMessage(e), "\n", sep = "", file = stderr())
  quit(save = "no", status = 1)
})
