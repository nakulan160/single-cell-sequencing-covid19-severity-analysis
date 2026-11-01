#!/usr/bin/env Rscript
# GSE196198 | Kim et al. (2023) | STEP 02: Seurat PBMC PCA / clusters / UMAP
#
# Run from the directory with the original GSE196198 files:
#   Rscript 2_umap.R .
#
# DEFAULT: reproducible 40,000-cell proportional stratified sample, rather than
# forcing the 399,044-cell dense processed matrix into RAM.
# Change only if your machine has enough memory:
#   GSE_MAX_CELLS=20000 Rscript 2_umap.R .
#   GSE_MAX_CELLS=80000 Rscript 2_umap.R .
#
# NOTE: the GEO integrated expression table contains signed, noninteger values.
# We import it as a Seurat Assay5 DATA layer, NOT raw counts. No NormalizeData,
# SCTransform, FindVariableFeatures(vst), RPCA or raw-count differential testing
# is performed here. This is an exploratory atlas reconstruction, NOT an exact
# reproduction of Supplementary Figure S2A.
#
# To regenerate figures from an already-saved Step 2 checkpoint (no reload):
#   GSE_REUSE_CHECKPOINT=1 Rscript 2_umap.R .
# Outputs: 2_results/ and a shareable 2_results.zip. The full Seurat RDS stays
# in 2_results/ on the user's machine and is deliberately excluded from the zip.

options(stringsAsFactors = FALSE, scipen = 99)
Sys.setenv(OMP_NUM_THREADS = "10", OPENBLAS_NUM_THREADS = "10",
           MKL_NUM_THREADS = "10", VECLIB_MAXIMUM_THREADS = "10")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  datadir <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
  outdir <- file.path(datadir, "2_results")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  logfile <- file.path(outdir, "2_umap_run_log.txt")
  notes <- character()
  logmsg <- function(...) {
    m <- paste0(..., collapse = "")
    cat(m, "\n", sep = "")
    notes <<- c(notes, paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), m))
    writeLines(notes, logfile)
  }
  logmsg("STEP 02 | GSE196198 | Seurat PBMC UMAP")
  logmsg("Data directory: ", datadir)
  logmsg("R version: ", R.version.string)

  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Seurat is not installed in THIS R environment. Fix Seurat installation first; then rerun. See 2_README.txt.")
  }
  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    stop("This processed-data-only import requires SeuratObject >= 5.0.0 (Assay5 data-only layer).")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  suppressPackageStartupMessages(library(Seurat))
  logmsg("Seurat version: ", as.character(utils::packageVersion("Seurat")))
  logmsg("SeuratObject version: ", as.character(utils::packageVersion("SeuratObject")))
  logmsg("Using prenormalized, signed integrated values as DATA, never as counts.")

  maxcells <- suppressWarnings(as.integer(Sys.getenv("GSE_MAX_CELLS", "40000")))
  if (is.na(maxcells) || maxcells < 1000L) stop("GSE_MAX_CELLS must be an integer >= 1000")
  expected_genes <- 1800L # Previously audited PBMC GEO export: 1,800 genes.
  seed <- 196198L
  set.seed(seed)
  setwd(datadir)

  # Optional recovery after a previous run has saved the checkpoint.
  # In a normal first run, GSE_REUSE_CHECKPOINT is not set.
  checkpoint <- file.path(outdir, "2_pbmc_seurat_sampled.rds")
  reuse_checkpoint <- identical(Sys.getenv("GSE_REUSE_CHECKPOINT", "0"), "1")
  if (reuse_checkpoint) {
    if (!file.exists(checkpoint)) {
      stop("Cannot resume: no Seurat checkpoint found at ", checkpoint,
           ". Run without GSE_REUSE_CHECKPOINT first.")
    }
    logmsg("RESUME MODE: reading previously saved Seurat object; skipping matrix/PCA/UMAP.")
    obj <- readRDS(checkpoint)
    if (!inherits(obj, "Seurat") || !"umap" %in% names(obj@reductions) ||
        !"analysis_clusters" %in% colnames(obj@meta.data)) {
      stop("The checkpoint is not a completed Step 2 Seurat UMAP object")
    }
    wanted <- ncol(obj)
    n_total <- obj@misc$gse_total_cells
    if (is.null(n_total) || length(n_total) != 1L ||
        !identical(as.integer(obj@misc$gse_seed), seed) ||
        !identical(as.integer(obj@misc$gse_max_cells), maxcells)) {
      stop("Checkpoint settings differ; run a fresh analysis without resume mode")
    }
    logmsg("Restored ", wanted, " cells and ", nrow(obj), " genes; ",
           length(unique(obj$analysis_clusters)), " clusters.")
  } else {

  metadata_file <- file.path(datadir, "GSE196198_pbmc_meta.data.txt.gz")
  expression_file <- file.path(datadir, "GSE196198_pbmc_normalized_count_integrated.txt.gz")
  if (!file.exists(metadata_file) || !file.exists(expression_file)) {
    stop("Missing PBMC metadata and/or expression file from GSE196198 in: ", datadir)
  }

  # The first column contains cell IDs, with GEO's quoted text headers.
  logmsg("Reading PBMC metadata...")
  meta_con <- gzfile(metadata_file, "rt")
  meta <- tryCatch(
    read.table(meta_con, header = TRUE, row.names = 1, sep = "",
               quote = "\"", comment.char = "", check.names = FALSE,
               stringsAsFactors = FALSE),
    finally = close(meta_con)
  )
  required <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                "seurat_clusters", "severity", "id.match", "severity2")
  missing <- setdiff(required, colnames(meta))
  if (length(missing)) stop("PBMC metadata missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(rownames(meta))) stop("Duplicated metadata row IDs")
  n_total <- nrow(meta)
  logmsg("PBMC metadata cells: ", n_total)

  # Check Seurat paper QC, already expected to be satisfied in the GEO export.
  numeric_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  for (nm in numeric_cols) meta[[nm]] <- suppressWarnings(as.numeric(meta[[nm]]))
  if (any(!is.finite(as.matrix(meta[, numeric_cols, drop = FALSE])))) {
    stop("Non-finite QC value in metadata")
  }
  qc_pass <- (meta$nFeature_RNA > 200 & meta$nFeature_RNA < 3500 &
              meta$nCount_RNA < 15000 & meta$percent.mt < 15)
  logmsg("QC eligible cells: ", sum(qc_pass), "/", n_total)

  # Read ONE header line, verify 1:1 cell match despite '-' vs '.' conventions.
  logmsg("Reading the huge expression-file header (no full matrix read)...")
  expr_con <- gzfile(expression_file, "rt")
  on.exit(try(close(expr_con), silent = TRUE), add = TRUE)
  header_line <- readLines(expr_con, n = 1L, warn = FALSE)
  if (!length(header_line)) stop("Expression matrix header is empty")
  # GEO write.table exports may be TAB, COMMA or whitespace-delimited.
  # This particular matrix was previously audited as having whitespace tokens.
  # Do not require a tab/comma; do not split single spaces if repeated spaces occur.
  sep_type <- if (grepl("\t", header_line, fixed = TRUE)) "tab" else
              if (grepl(",", header_line, fixed = TRUE)) "comma" else
              if (grepl("[[:space:]]", header_line)) "whitespace" else "unknown"
  if (sep_type == "unknown") {
    stop("No recognizable delimiters in PBMC expression header. First 100 characters: ",
         substr(header_line, 1L, 100L))
  }
  logmsg("Detected expression separator: ", sep_type)
  parse_fields <- function(line) {
    line <- sub("\r$", "", line)
    if (sep_type == "tab") return(strsplit(line, "\t", fixed = TRUE)[[1L]])
    if (sep_type == "comma") return(strsplit(line, ",", fixed = TRUE)[[1L]])
    strsplit(trimws(line), "[[:space:]]+")[[1L]]
  }
  strip_quotes <- function(x) sub('"$', "", sub('^"', "", x))
  col_headers <- strip_quotes(parse_fields(header_line))
  col_headers[[1L]] <- sub("^\ufeff", "", col_headers[[1L]])
  logmsg("Expression header parsed fields: ", length(col_headers),
         " (expected ", n_total, " or ", n_total + 1L, ")")
  if (length(col_headers) == n_total + 1L) col_headers <- col_headers[-1L]
  if (length(col_headers) != n_total) {
    stop("RNA header length ", length(col_headers), " != metadata cells ", n_total)
  }
  normalize_ids <- function(x) gsub("[^A-Za-z0-9_.]", ".", x)
  norm_meta <- normalize_ids(rownames(meta))
  norm_hdr <- normalize_ids(col_headers)
  if (anyDuplicated(norm_meta) || anyDuplicated(norm_hdr)) {
    stop("Collision after punctuation normalization of cell names")
  }
  corresponding <- match(norm_hdr, norm_meta)
  if (anyNA(corresponding) || length(corresponding) != nrow(meta) ||
      anyDuplicated(corresponding)) {
    stop("RNA cell IDs do not match metadata one-to-one. Matrix examples: ",
         paste(head(col_headers, 3), collapse = ", "), "; metadata examples: ",
         paste(head(rownames(meta), 3), collapse = ", "))
  }
  # Align cells to RNA-column order; never assume metadata order is the same.
  meta <- meta[corresponding, , drop = FALSE]
  eligible <- which(qc_pass[corresponding])
  if (!length(eligible)) stop("No QC-eligible cells")
  logmsg("RNA cell IDs matched to metadata: ", length(corresponding), " exactly")

  # Check gene-row formatting before allocating a large expression matrix.
  # This uses a separate connection so the main stream remains at the first gene.
  check_con <- gzfile(expression_file, "rt")
  preview <- tryCatch({
    readLines(check_con, n = 1L, warn = FALSE)  # header
    readLines(check_con, n = 1L, warn = FALSE)  # first gene row
  }, finally = close(check_con))
  if (!length(preview)) stop("Expression file has a header but no gene rows")
  first_fields <- parse_fields(preview)
  if (length(first_fields) != n_total + 1L) {
    stop("First gene row has ", length(first_fields), " fields; expected ",
         n_total + 1L, ". File format requires inspection before analysis.")
  }
  logmsg("First gene row: ", strip_quotes(first_fields[[1L]]),
         " | parsed values: ", length(first_fields) - 1L)
  rm(preview, first_fields)
  invisible(gc(verbose = FALSE))

  # Select cells proportionally from each of the 25 collections.
  # For exact whole-cohort UMAP set GSE_MAX_CELLS >= 399044 ONLY if RAM allows.
  wanted <- min(length(eligible), maxcells)
  groups <- split(eligible, meta$orig.ident[eligible])
  raw_quota <- lengths(groups) * wanted / length(eligible)
  quota <- floor(raw_quota)
  remaining <- wanted - sum(quota)
  if (remaining > 0L) {
    which_extra <- order(raw_quota - quota, decreasing = TRUE)[seq_len(remaining)]
    quota[which_extra] <- quota[which_extra] + 1L
  }
  sampled <- sort(unlist(Map(function(idx, n) {
    if (n == 0L) integer(0) else sample(idx, size = n, replace = FALSE)
  }, groups, quota), use.names = FALSE))
  if (length(sampled) != wanted || anyDuplicated(sampled)) stop("Internal sampling error")
  selected_metadata <- meta[sampled, , drop = FALSE]
  rownames(selected_metadata) <- col_headers[sampled]
  selected_metadata$author_clusters <- as.character(selected_metadata$seurat_clusters)
  selected_metadata$sample_id <- as.character(selected_metadata$orig.ident)
  selected_metadata$donor_code_inferred <- sub("_.*$", "", sub("^Exp_", "", selected_metadata$sample_id))
  selected_metadata$visit_code_inferred <- ifelse(
    grepl("_", sub("^Exp_", "", selected_metadata$sample_id)),
    sub("^[^_]+_", "", sub("^Exp_", "", selected_metadata$sample_id)), "")
  selected_metadata$severity_broad <- sub("_.*$", "", as.character(selected_metadata$severity))
  summary_sampling <- data.frame(
    sample = names(groups), total_qc_cells = as.integer(lengths(groups)),
    selected_cells = as.integer(quota),
    sampling_fraction = as.numeric(quota) / as.numeric(lengths(groups)),
    row.names = NULL
  )
  write.csv(summary_sampling, file.path(outdir, "2_sampling_by_sample.csv"), row.names = FALSE)
  logmsg("Proportionally stratified sampled cells: ", wanted, " / ", length(eligible),
         " (", sprintf("%.1f", 100 * wanted / length(eligible)), "%)")
  logmsg("WARNING: This sampled UMAP cannot be claimed as an exact paper Figure S2A reproduction.")

  # Stream all gene rows; store ONLY selected cells to control memory.
  # CSV/TSV entries are simple quoted/unquoted fields; no embedded tabs/commas
  # are expected in cell IDs or numeric expression values.
  memory_GiB <- expected_genes * wanted * 8 / 1024^3
  logmsg("Estimated ONE dense double matrix: ", sprintf("%.2f GiB", memory_GiB),
         "; PCA and Seurat need more RAM for working copies.")
  expr <- matrix(NA_real_, nrow = expected_genes, ncol = wanted)
  selected_gene_names <- character(expected_genes)
  step <- 0L
  logmsg("Streaming expression table (this can take several minutes)...")
  repeat {
    next_line <- readLines(expr_con, n = 1L, warn = FALSE)
    if (!length(next_line)) break
    step <- step + 1L
    if (step > expected_genes) stop("More than ", expected_genes, " genes found; check format")
    fields <- parse_fields(next_line)
    if (length(fields) != n_total + 1L) {
      stop("Wrong number of values on expression row ", step, ": ",
           length(fields), " vs expected ", n_total + 1L)
    }
    selected_gene_names[step] <- strip_quotes(fields[[1L]])
    values <- suppressWarnings(as.numeric(strip_quotes(fields[sampled + 1L])))
    if (length(values) != wanted || any(!is.finite(values))) {
      stop("Non-numeric / missing data in gene row ", step,
           " (", selected_gene_names[step], ")")
    }
    expr[step, ] <- values
    if (step %% 100L == 0L) logmsg("Read gene rows: ", step, "/", expected_genes)
  }
  close(expr_con)
  if (step != expected_genes) stop("Expected 1800 genes; read ", step)
  if (anyDuplicated(selected_gene_names)) stop("Duplicated gene identifiers in PBMC export")
  if (any(!nzchar(selected_gene_names))) stop("Empty gene identifiers")
  rownames(expr) <- gsub("_", "-", selected_gene_names, fixed = TRUE)
  if (anyDuplicated(rownames(expr))) stop("Gene identifiers collided after Seurat naming convention")
  colnames(expr) <- col_headers[sampled]
  if (!identical(colnames(expr), rownames(selected_metadata))) stop("Expression metadata columns misaligned")
  if (!any(expr < 0)) logmsg("WARNING: Sample had no negative values; check transform separately")
  logmsg("Loaded expression rows/columns: ", nrow(expr), " x ", ncol(expr))
  logmsg("Sampled expression range: ", sprintf("%.4f", min(expr)), " to ",
         sprintf("%.4f", max(expr)))

  # Seurat v5 permits a data-only assay, with no invented UMI count layer.
  # This is descriptive reclustering of the PROCESSED authors' expression values.
  logmsg("Creating Seurat data-only assay and object...")
  assay <- SeuratObject::CreateAssay5Object(data = expr)
  rm(expr)
  invisible(gc(verbose = FALSE))
  obj <- SeuratObject::CreateSeuratObject(
    counts = assay, assay = "RNA", meta.data = selected_metadata,
    project = "GSE196198_ProcessedPBMC"
  )
  rm(assay, meta, selected_metadata, groups, corresponding)
  invisible(gc(verbose = FALSE))
  DefaultAssay(obj) <- "RNA"
  # Preserve published QC fields; a data-only assay cannot recompute them.
  # Do NOT call NormalizeData, SCTransform, or FindVariableFeatures(vst).
  features <- rownames(obj)
  logmsg("Scaling ", length(features), " supplied genes for exploratory PCA...")
  obj <- Seurat::ScaleData(obj, features = features, do.center = TRUE,
                           do.scale = TRUE, scale.max = 10, verbose = FALSE)
  max_pcs <- min(30L, nrow(obj) - 1L, ncol(obj) - 1L)
  obj <- Seurat::RunPCA(obj, features = features, npcs = max_pcs,
                        seed.use = seed, verbose = FALSE)
  # The paper's exact subset clustering parameters are not fully verified.
  # These are explicitly ANALYST settings, not quoted as original settings.
  use_dims <- seq_len(min(30L, ncol(Seurat::Embeddings(obj, "pca"))))
  logmsg("PCA ready. Building neighbor graph using ", length(use_dims), " PCs...")
  obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = use_dims,
                               verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 0.5, algorithm = 1,
                              random.seed = seed, verbose = FALSE)
  obj$analysis_clusters <- as.character(SeuratObject::Idents(obj))
  logmsg("Running UMAP...")
  obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = use_dims,
                         n.neighbors = 30L, min.dist = 0.3,
                         seed.use = seed, verbose = FALSE)
  logmsg("Clusters found: ", length(unique(obj$analysis_clusters)))

  # Preserve finished computation BEFORE any plotting. A failed plot can be
  # retried without rereading 399,044 RNA column entries for all 1,800 genes.
  obj@misc$gse_total_cells <- as.integer(n_total)
  obj@misc$gse_max_cells <- as.integer(maxcells)
  obj@misc$gse_seed <- as.integer(seed)
  logmsg("Saving completed UMAP checkpoint BEFORE figures (may take time)...")
  saveRDS(obj, file = checkpoint, compress = "gzip")
  logmsg("Saved local Seurat checkpoint: ", checkpoint)
  } # End of fresh analysis; either branch now provides obj, wanted, n_total.

  save_png_plot <- function(plot, filename, width = 12, height = 8) {
    ggplot2::ggsave(filename = file.path(outdir, filename), plot = plot,
                    width = width, height = height, units = "in", dpi = 175,
                    limitsize = FALSE)
    logmsg("Saved figure: ", filename)
  }
  plot_umap <- function(field, title, filename, label = FALSE) {
    p <- Seurat::DimPlot(obj, reduction = "umap", group.by = field,
                         raster = requireNamespace("scattermore", quietly = TRUE),
                         label = label, repel = label,
                         pt.size = 0.02) + ggplot2::ggtitle(title)
    save_png_plot(p, filename)
  }
  plot_umap("analysis_clusters", "New Seurat clusters (sampled processed PBMC)",
            "2_umap_new_clusters.png", label = TRUE)
  plot_umap("severity_broad", "COVID-19 broad severity labels (sampled cells)",
            "2_umap_severity.png")
  plot_umap("donor_code_inferred", "Sampled PBMC cells by inferred donor",
            "2_umap_donors.png")
  plot_umap("author_clusters", "Author-supplied cluster IDs (not new clustering)",
            "2_umap_author_clusters.png")

  sdev <- SeuratObject::Stdev(obj, reduction = "pca")
  variance_df <- data.frame(PC = seq_along(sdev),
                            explained_pct = 100 * sdev^2 / sum(sdev^2))
  write.csv(variance_df, file.path(outdir, "2_pca_variance.csv"), row.names = FALSE)
  pvar <- ggplot2::ggplot(variance_df, ggplot2::aes(x = PC, y = explained_pct)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "PCA variance (30 computed components)",
                  y = "% of variance within computed PCs")
  save_png_plot(pvar, "2_pca_variance.png", 9, 5)

  requested_markers <- c("CD3D", "CD3E", "TRAC", "MS4A1", "CD79A",
                         "NKG7", "GNLY", "LYZ", "S100A8", "LILRA4")
  found <- intersect(requested_markers, rownames(obj))
  write.csv(data.frame(marker = requested_markers,
                       present = requested_markers %in% rownames(obj)),
            file.path(outdir, "2_requested_marker_availability.csv"), row.names = FALSE)
  if (length(found) >= 2L) {
    # Seurat's order=TRUE + raster=TRUE combination requires ggrastr.
    # Use non-raster ordered points to avoid that optional dependency.
    logmsg("Rendering gene markers without ordered rasterization (no ggrastr required)...")
    tryCatch({
      marker_plot <- Seurat::FeaturePlot(
        obj, features = head(found, 8L), reduction = "umap", ncol = 4L,
        pt.size = 0.02, raster = FALSE, order = TRUE
      )
      save_png_plot(marker_plot, "2_umap_marker_preview.png", 16, 8)
    }, error = function(e) {
      # The global atlas and checkpoint remain usable even if plotting fails.
      logmsg("WARNING: Optional marker plot could not be saved: ", conditionMessage(e))
      writeLines(conditionMessage(e), file.path(outdir, "2_marker_plot_warning.txt"))
    })
  } else {
    logmsg("Marker preview skipped: fewer than 2 requested markers in 1800-gene export")
  }

  emb <- as.data.frame(Seurat::Embeddings(obj, reduction = "umap"))
  emb$cell_id <- rownames(emb)
  emb$sample_id <- obj$sample_id
  emb$donor_code <- obj$donor_code_inferred
  emb$severity <- obj$severity
  emb$severity2 <- obj$severity2
  emb$severity_broad <- obj$severity_broad
  emb$analysis_clusters <- obj$analysis_clusters
  emb$author_clusters <- obj$author_clusters
  embed_con <- gzfile(file.path(outdir, "2_umap_coordinates_sampled.csv.gz"), "wt")
  utils::write.csv(emb, embed_con, row.names = FALSE)
  close(embed_con)

  counts <- as.data.frame(table(analysis_clusters = obj$analysis_clusters,
                                 severity_broad = obj$severity_broad))
  write.csv(counts, file.path(outdir, "2_cluster_by_severity_descriptive.csv"),
            row.names = FALSE)
  cell_counts <- as.data.frame(table(analysis_clusters = obj$analysis_clusters))
  names(cell_counts) <- c("analysis_cluster", "n_sampled_cells")
  write.csv(cell_counts, file.path(outdir, "2_cluster_sizes.csv"), row.names = FALSE)

  # Checkpoint was saved immediately after UMAP (fresh mode) or restored
  # from that file (resume mode). It is excluded from the upload ZIP.
  writeLines(capture.output(sessionInfo()), file.path(outdir, "2_session_info.txt"))
  writeLines(c(
    "METHOD: exploratory de novo Seurat reclustering and UMAP of GEO processed expression.",
    "INPUT: GSE196198_pbmc_normalized_count_integrated.txt.gz and matching metadata.",
    "NO RAW COUNT LAYER. No NormalizeData, SCTransform, FindVariableFeatures, RPCA, or DE.",
    paste("Cells selected:", wanted, "of", n_total, "| seed:", seed),
    paste("Genes read:", expected_genes),
    "Sampling: proportionate stratified random sample from each GEO collection.",
    "PCA: scaled all 1800 supplied genes; PCA 30 PCs; neighbor graph using 30 PCs.",
    "Clustering: Louvain algorithm 1, resolution 0.5, seed 196198.",
    "UMAP: 30 PCs, 30 neighbors, min.dist 0.3, seed 196198.",
    "These analysis parameters are NOT verified as an exact copy of Kim et al. figure-generation settings.",
    "New clustering is not an official or exact reproduction of Supplementary Figure S2.",
    "Group and cluster cell counts are DESCRIPTIVE and NOT donor-level significance tests.",
    "Author cluster IDs were preserved in author_clusters for comparison only."
  ), file.path(outdir, "2_methods_and_cautions.txt"))

  # Produce portable ZIP without repeating the large Seurat object.
  zip_exe <- Sys.which("zip")
  archive_path <- file.path(datadir, "2_results.zip")
  if (nzchar(zip_exe)) {
    shareable <- list.files("2_results", full.names = TRUE, recursive = FALSE)
    shareable <- shareable[!dir.exists(shareable)]
    shareable <- shareable[basename(shareable) != "2_pbmc_seurat_sampled.rds"]
    if (file.exists(archive_path)) unlink(archive_path)
    status <- suppressWarnings(system2(zip_exe,
         args = c("-q", "-9", shQuote(archive_path), shQuote(shareable))))
    if (!identical(as.integer(status), 0L) || !file.exists(archive_path)) {
      logmsg("WARNING: ZIP creation failed; manually zip 2_results except the RDS checkpoint")
    } else {
      logmsg("Created upload archive: ", archive_path)
    }
  } else {
    logmsg("WARNING: 'zip' executable not found. Install it or create 2_results.zip manually.")
  }
  logmsg("SUCCESS: next send 2_results.zip; keep 2_pbmc_seurat_sampled.rds locally.")
}

tryCatch(main(), error = function(e) {
  msg <- paste0("STEP 02 FAILED: ", conditionMessage(e))
  cat(msg, "\n", file = stderr())
  try({
    dd <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[[1]] else "."
    od <- file.path(dd, "2_results")
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    writeLines(c(msg, capture.output(traceback())), file.path(od, "2_error.txt"))
  }, silent = TRUE)
  quit(save = "no", status = 1)
})
