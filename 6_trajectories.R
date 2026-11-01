#!/usr/bin/env Rscript
# GSE196198 | STEP 06 | Donor-aware exploratory T/B cluster trajectories
# Input: the completed Step 4 and Step 5 per-sample CSV tables.
# Reuses saved analysis tables; does NOT reload large matrices or compute DE.
options(stringsAsFactors = FALSE, scipen = 99)
Sys.setenv(OMP_NUM_THREADS="10", OPENBLAS_NUM_THREADS="10", MKL_NUM_THREADS="10")

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
  setwd(root)
  out <- file.path(root, "6_results")
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  logpath <- file.path(out, "6_run_log.txt")
  if (file.exists(logpath)) unlink(logpath)
  msg <- function(...) {
    txt <- paste0(..., collapse="")
    cat(txt, "\n")
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt, "\n", file=logpath, append=TRUE)
  }
  msg("STEP 06 | GSE196198 | Patient-aware B/T-cell cluster trajectories")
  if (!requireNamespace("ggplot2", quietly=TRUE)) stop("ggplot2 missing")

  read_step <- function(n) {
    path <- file.path(root, sprintf("%d_results/%d_clusters_by_sample_descriptive.csv",n,n))
    if (!file.exists(path)) stop("Missing file: ",path)
    d <- read.csv(path, check.names=FALSE)
    wanted <- c("sample_id", "cluster", "n_cells", "sample_total", "fraction",
                "donor_inferred", "severity2")
    if (!all(wanted %in% names(d))) stop("Unexpected Step ",n," table columns")
    if (anyNA(d[,wanted]) || anyDuplicated(d[c("sample_id","cluster")]))
      stop("Missing or duplicate sample/cluster records for Step ", n)
    d$cluster <- as.character(d$cluster)
    if (any(d$n_cells < 0) || any(d$sample_total < 1) ||
        any(d$n_cells != round(d$n_cells)) ||
        any(d$sample_total != round(d$sample_total)) ||
        any(abs(d$fraction - d$n_cells/d$sample_total) > 1e-7))
      stop("Bad sample counts/fractions in Step ",n)
    sm <- unique(d[c("sample_id", "sample_total", "donor_inferred", "severity2")])
    if (anyDuplicated(sm$sample_id)) stop("Conflicting sample metadata in Step ",n)
    totals <- aggregate(n_cells ~ sample_id, data=d, sum)
    check <- match(sm$sample_id, totals$sample_id)
    if (anyNA(check) || any(totals$n_cells[check] != sm$sample_total))
      stop("Cluster counts do not sum to sample totals in Step ", n)
    msg("Step ",n,": ", sum(sm$sample_total)," cells in ",nrow(sm)," collections and ",
        length(unique(sm$donor_inferred))," donors")
    list(d=d, sm=sm)
  }
  T <- read_step(4)
  B <- read_step(5)
  if (!setequal(T$sm$sample_id, B$sm$sample_id))
    stop("B and T samples do not match")
  matched <- match(T$sm$sample_id, B$sm$sample_id)
  if (!identical(as.character(T$sm$severity2), as.character(B$sm$severity2[matched])) ||
      !identical(as.character(T$sm$donor_inferred), as.character(B$sm$donor_inferred[matched])))
    stop("B and T clinical/donor records disagree")

  group_info <- function(x) {
    x <- as.character(x)
    out <- ifelse(x=="Healthy", "Healthy", ifelse(x=="Asymptomatic", "Asymptomatic",
           ifelse(grepl("^Mild_[124]$", x), "Mild",
           ifelse(grepl("^Severe\\.m_[124]$", x), "Severe.m",
           ifelse(grepl("^Severe\\.s_[124]$", x), "Severe.s", NA_character_)))))
    if (anyNA(out)) stop("Unrecognized severity2 labels: ", paste(unique(x[is.na(out)]), collapse=", "))
    out
  }
  metadata <- T$sm[c("sample_id", "donor_inferred", "severity2")]
  metadata$disease_group <- group_info(metadata$severity2)
  metadata$visit_code <- ifelse(metadata$disease_group %in% c("Healthy", "Asymptomatic"),
                                0L, suppressWarnings(as.integer(sub(".*_([124])$", "\\1", metadata$severity2))))
  if (anyNA(metadata$visit_code)) stop("Cannot decode visit suffixes")
  if (anyDuplicated(metadata[c("donor_inferred","visit_code")]))
    stop("Duplicate donor and visit records")
  donor_groups <- aggregate(disease_group ~ donor_inferred, metadata, function(x) length(unique(x)))
  if (any(donor_groups$disease_group != 1)) stop("Inconsistent disease labels within donor")
  write.csv(metadata[order(metadata$donor_inferred,metadata$visit_code),],
            file.path(out,"6_sample_donor_clinical_index.csv"), row.names=FALSE)

  # These clusters were identified AFTER inspecting the Step 4/5 marker heatmaps.
  # They are exploratory candidate programs, not validated immunophenotypes.
  programs <- list(
    T_Cytotoxic_2_10 = list(source="T", clusters=c("2","10"), evidence="GNLY;GZMB;PRF1"),
    T_Cycling_13 = list(source="T", clusters="13", evidence="MKI67;CD38"),
    T_Interferon_16 = list(source="T", clusters="16", evidence="ISG15;IFIT1;MX1"),
    B_Secretory_13 = list(source="B", clusters="13", evidence="JCHAIN;MZB1;XBP1;PRDM1"),
    B_Atypical_6 = list(source="B", clusters="6", evidence="FCRL5;ITGAX"),
    B_Interferon_14 = list(source="B", clusters="14", evidence="ISG15;IFIT3;MX1")
  )
  definitions <- do.call(rbind,lapply(names(programs),function(p) {
    a <- programs[[p]]
    data.frame(program=p, source=a$source, cluster_ids=paste(a$clusters, collapse="+"),
               supporting_markers=a$evidence)
  }))
  write.csv(definitions,file.path(out,"6_exploratory_program_definitions.csv"),row.names=FALSE)

  make_program <- function(name) {
    spec <- programs[[name]]
    input <- if (spec$source=="T") T else B
    sm <- input$sm
    selected <- input$d[input$d$cluster %in% spec$clusters, c("sample_id","n_cells")]
    sums <- if (nrow(selected)) aggregate(n_cells ~ sample_id, selected, sum) else
      data.frame(sample_id=character(),n_cells=numeric())
    lookup <- match(sm$sample_id, sums$sample_id)
    n <- sums$n_cells[lookup]
    n[is.na(n)] <- 0
    d <- data.frame(sample_id=sm$sample_id, program=name, lineage=spec$source,
                    n_program=as.integer(n), n_lineage_sampled=as.integer(sm$sample_total))
    d$percentage <- 100*d$n_program/d$n_lineage_sampled
    merge(d,metadata,by="sample_id",sort=FALSE)
  }
  result <- do.call(rbind,lapply(names(programs),make_program))
  result <- result[order(result$program,result$donor_inferred,result$visit_code),]
  rownames(result) <- NULL
  if (nrow(result) != 25L*length(programs) || anyNA(result)) stop("Incomplete sample-program table")
  write.csv(result,file.path(out,"6_sample_level_program_fractions.csv"),row.names=FALSE)

  # Summarize DONOR/sample proportions, NEVER treat cells as independent donors.
  grp <- split(result,interaction(result$program,result$disease_group,result$visit_code,drop=TRUE))
  summaries <- do.call(rbind, lapply(grp,function(d) {
    v <- d$percentage
    data.frame(program=d$program[1],lineage=d$lineage[1],
      disease_group=d$disease_group[1],visit_code=d$visit_code[1],
      n_donors=length(unique(d$donor_inferred)),
      mean_percent=mean(v),median_percent=median(v),min_percent=min(v),
      q1_percent=as.numeric(quantile(v,.25)),q3_percent=as.numeric(quantile(v,.75)),
      max_percent=max(v))
  }))
  rownames(summaries) <- NULL
  write.csv(summaries,file.path(out,"6_by_group_visit_descriptive.csv"),row.names=FALSE)

  longitudinal <- result[result$visit_code>0,]
  paired <- split(longitudinal, interaction(longitudinal$donor_inferred,longitudinal$program,drop=TRUE))
  paired_rows <- do.call(rbind,lapply(paired,function(d) {
    d <- d[order(d$visit_code),]
    if (nrow(d)<2L) return(NULL)
    data.frame(donor_inferred=d$donor_inferred[1],disease_group=d$disease_group[1],
               program=d$program[1],lineage=d$lineage[1],
               first_visit=d$visit_code[1],last_visit=d$visit_code[nrow(d)],
               first_percent=d$percentage[1],last_percent=d$percentage[nrow(d)],
               delta_percentage_points=d$percentage[nrow(d)]-d$percentage[1],
               n_observed_visits=nrow(d))
  }))
  if (is.null(paired_rows)) stop("No paired longitudinal data")
  rownames(paired_rows)<-NULL
  write.csv(paired_rows,file.path(out,"6_within_donor_first_to_last_change.csv"),row.names=FALSE)

  save_plot <- function(p,n,width=12,height=6) {
    ggplot2::ggsave(file.path(out,n),p,width=width,height=height,dpi=160,units="in",limitsize=FALSE)
    msg("Figure: ",n)
  }
  programs_to_plot <- c("B_Secretory_13", "B_Atypical_6", "B_Interferon_14",
                        "T_Cytotoxic_2_10", "T_Cycling_13", "T_Interferon_16")
  for (name in programs_to_plot) {
    d <- longitudinal[longitudinal$program==name,]
    d$disease_group <- factor(d$disease_group,levels=c("Mild","Severe.m","Severe.s"))
    p <- ggplot2::ggplot(d,ggplot2::aes(x=visit_code,y=percentage,
                       group=donor_inferred,color=donor_inferred)) +
         ggplot2::geom_line(linewidth=.9,alpha=.75) +
         ggplot2::geom_point(size=2.5) +
         ggplot2::facet_grid(~disease_group) +
         ggplot2::scale_x_continuous(breaks=c(1,2,4),labels=c("1","2","4")) +
         ggplot2::scale_y_continuous(limits=c(0,NA)) +
         ggplot2::labs(title=gsub("_"," ",name),
                       subtitle="Each line represents one donor (exploratory clustering)",
                       x="Supplied collection code (NOT confirmed days)",
                       y="% of sampled cells within this lineage", color="Donor") +
         ggplot2::theme_classic(base_size=12)
    save_plot(p,paste0("6_trajectory_",name,".png"),width=13,height=5)
  }
  # The first-collection plot includes all donors with an available baseline.
  first <- result[result$visit_code %in% c(0,1) &
                    result$program %in% c("B_Secretory_13","T_Cytotoxic_2_10"),]
  first$disease_group <- factor(first$disease_group,
    levels=c("Healthy","Asymptomatic","Mild","Severe.m","Severe.s"))
  p_first <- ggplot2::ggplot(first,ggplot2::aes(x=disease_group,y=percentage)) +
    ggplot2::geom_point(ggplot2::aes(color=donor_inferred),
                        size=3,position=ggplot2::position_jitter(width=.09,height=0),alpha=.9) +
    ggplot2::facet_wrap(~program,scales="free_y",ncol=1) +
    ggplot2::labs(title="Initial captured cell fractions: donor-level observations",
                  subtitle="Not a significance test; cross-lineage percentages have different denominators",
                  x="Original clinical group labels",y="% of sampled cells within lineage",color="Donor") +
    ggplot2::theme_classic(base_size=12) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=22,hjust=1))
  save_plot(p_first,"6_first_collection_B_secretory_T_cytotoxic.png",12,9)

  writeLines(c(
    "GSE196198 STEP 06: PRELIMINARY DONOR-AWARE DESCRIPTIVE TRAJECTORIES",
    "Data sources: completed 4_results/ and 5_results/ per-sample cluster tables.",
    paste("B and T samples:",nrow(metadata),"| inferred donors:",length(unique(metadata$donor_inferred))),
    "All proportions have lineage-specific denominators (sampled T or sampled B cells).",
    "Within each donor, codes 1, 2, 4 are supplied collection identifiers, not confirmed days.",
    "Do not interpret Severe.m and Severe.s as specified clinical trajectories without published metadata confirmation.",
    "Healthy and asymptomatic donors each contribute one collection. Mild/severe donors contribute repeated collections.",
    "Candidate clusters were selected POST HOC from Step 4 and Step 5 RNA marker heatmaps; labels require validation.",
    "Group summaries use individual donor/sample fractions; NOT independent single-cell observations.",
    "No hypothesis tests, p-values, causal claims, normalized UMI counts, or differential-expression analysis.",
    "Processed GEO integrated expression values are NOT a suitable replacement for raw count-based DE.",
    "Differences may reflect donor heterogeneity, variable cell capture, sampling, and treatment/timing confounders.",
    "Plots and cell frequencies support exploratory replication, not an exact reproduction of the paper's main figures."
  ),file.path(out,"6_methods_and_cautions.txt"))
  writeLines(capture.output(sessionInfo()),file.path(out,"6_session_info.txt"))
  failure <- file.path(out,"6_error.txt")
  if (file.exists(failure)) unlink(failure)
  msg("All figure/table generation complete; packaging ZIP...")
  zip <- Sys.which("zip")
  if (nzchar(zip)) {
    dst <- file.path(root,"6_results.zip")
    if (file.exists(dst)) unlink(dst)
    files <- list.files("6_results",full.names=TRUE)
    status <- suppressWarnings(system2(zip,c("-q","-9",shQuote(dst),shQuote(files))))
    if (identical(as.integer(status),0L) && file.exists(dst)) msg("Created archive: ",dst)
    else msg("WARNING: ZIP creation failed; outputs remain in 6_results/")
  } else msg("WARNING: system 'zip' missing; use zip -r 6_results.zip 6_results/")
  msg("SUCCESS: Step 06 completed. Send 6_results.zip for review.")
}

tryCatch(main(),error=function(e){
  txt<-paste0("STEP 06 FAILED: ",conditionMessage(e))
  cat(txt,"\n",file=stderr())
  args<-commandArgs(trailingOnly=TRUE)
  root<-if(length(args)) args[[1]] else "."
  od<-file.path(root,"6_results")
  dir.create(od,recursive=TRUE,showWarnings=FALSE)
  writeLines(txt,file.path(od,"6_error.txt"))
  quit(save="no",status=1)
})
