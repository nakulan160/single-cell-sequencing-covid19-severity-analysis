#!/usr/bin/env Rscript
# Kim et al. (2023), GSE196198 | STEP 03: PBMC lineage annotation
#
# PURPOSE
#   Connect the Step 2 Seurat UMAP to the ORIGINAL AUTHOR-PROVIDED PBMC
#   subpopulation metadata (T/B/NK/Monocyte). Assess cluster purity and
#   inspect marker signal on the existing signed processed RNA data layer.
#   Also calculate full-cohort (399,044 PBMC) lineage counts per collection.
#
# THIS IS NOT independent discovery of cell types, a raw count analysis, or
# an exact reproduction of the paper's Supplementary Figure S2.
#
# RUN (from directory holding 13 GSE196198 files):
#   env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
#     conda run --no-capture-output -n ngs_seurat \
#     Rscript --vanilla 3_annotation.R .
#
# INPUT: 2_results/2_pbmc_seurat_sampled.rds plus five GEO metadata tables.
# OUTPUT: 3_results/ and 3_results.zip; large annotated RDS stays local.

options(stringsAsFactors = FALSE, scipen = 99)
Sys.setenv(OMP_NUM_THREADS = '10', OPENBLAS_NUM_THREADS = '10', MKL_NUM_THREADS = '10')

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  datadir <- normalizePath(if (length(args)) args[1] else '.', mustWork = TRUE)
  outdir <- file.path(datadir, '3_results')
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  olderr <- file.path(outdir, '3_error.txt')
  if (file.exists(olderr)) unlink(olderr)
  log_path <- file.path(outdir, '3_annotation_run_log.txt')
  writeLines(character(), log_path)
  logmsg <- function(...) {
    message_text <- paste0(..., collapse = '')
    cat(message_text, '\n', sep = '')
    cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ', message_text,
        '\n', file = log_path, append = TRUE, sep = '')
  }

  logmsg('STEP 03 | GSE196198 | PBMC source-lineage annotation + verification')
  logmsg('Data directory: ', datadir)
  for (p in c('Seurat', 'SeuratObject', 'ggplot2', 'Matrix')) {
    if (!requireNamespace(p, quietly = TRUE)) stop('Missing R package: ', p)
  }
  suppressPackageStartupMessages(library(Seurat))
  logmsg('Seurat version: ', as.character(packageVersion('Seurat')))

  checkpoint <- file.path(datadir, '2_results', '2_pbmc_seurat_sampled.rds')
  if (!file.exists(checkpoint)) stop('Missing Step 2 local checkpoint: ', checkpoint,
                                    '\nThis file is NOT in 2_results.zip. Keep it on your computer.')
  logmsg('Loading completed Step 2 Seurat object...')
  obj <- readRDS(checkpoint)
  if (!inherits(obj, 'Seurat')) stop('Step 2 checkpoint is not a Seurat object')
  if (!'umap' %in% names(obj@reductions)) stop('Step 2 checkpoint has no UMAP')
  required_obj_cols <- c('analysis_clusters','sample_id','severity_broad','author_clusters')
  if (!all(required_obj_cols %in% colnames(obj@meta.data)))
    stop('Step 2 object lacks: ', paste(setdiff(required_obj_cols, colnames(obj@meta.data)), collapse=', '))
  if (length(unique(as.character(obj$analysis_clusters))) < 2L) stop('Unexpectedly few clusters')
  logmsg('Loaded ', ncol(obj), ' sampled cells; ', nrow(obj), ' supplied genes; ',
         length(unique(obj$analysis_clusters)), ' analysis clusters.')

  # The author-derived subset metadata may use different row-name conventions
  # from the PBMC object. Use the provided canonical sample-aware id.match
  # key for cross-assay matching, never normalized metadata row names.
  checked_id_match <- function(dat, context) {
    if (!('id.match' %in% colnames(dat)))
      stop(context, ' missing canonical id.match column')
    key <- as.character(dat[['id.match']])
    if (anyNA(key) || any(!nzchar(key)))
      stop(context, ' has missing or blank id.match values')
    if (anyDuplicated(key))
      stop(context, ' has duplicate id.match values (first duplicate: ',
           key[anyDuplicated(key)], ')')
    key
  }
  read_gz_meta <- function(filename) {
    path <- file.path(datadir, filename)
    if (!file.exists(path)) stop('Missing input metadata: ', path)
    con <- gzfile(path, 'rt')
    dat <- tryCatch(
      utils::read.table(con, header=TRUE, row.names=1, sep='', quote='"',
                        comment.char='', check.names=FALSE, stringsAsFactors=FALSE),
      finally=close(con))
    if (!('orig.ident' %in% names(dat))) stop(filename, ' missing orig.ident')
    if (anyDuplicated(rownames(dat))) stop(filename, ' has duplicate cell names')
    dat
  }

  logmsg('Reading complete PBMC metadata (full cohort)...')
  all_pbmc <- read_gz_meta('GSE196198_pbmc_meta.data.txt.gz')
  pbmc_keys <- checked_id_match(all_pbmc, 'PBMC metadata')
  logmsg('PBMC canonical id.match identifiers: ', length(pbmc_keys),
         ' unique; row names intentionally not used for subset matching.')
  sample_ids <- as.character(all_pbmc$orig.ident)
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) stop('Missing PBMC sample IDs')
  ref_names <- c('T cells','B cells','NK cells','Monocytes')
  ref_prefix <- c('T_cell','B_cell','NK_cell','Monocyte')
  membership <- matrix(FALSE, nrow=nrow(all_pbmc), ncol=4L,
                       dimnames=list(NULL, ref_names))
  source_summary <- data.frame(lineage=character(), source_cells=integer(),
                               unique_matched_cells=integer(), missing_ids=integer())
  for (j in seq_along(ref_prefix)) {
    file <- paste0('GSE196198_', ref_prefix[j], '_meta.data.txt.gz')
    logmsg('Loading author-provided subset: ', file)
    submeta <- read_gz_meta(file)
    keys <- checked_id_match(submeta, file)
    ix <- match(keys, pbmc_keys)
    if (anyNA(ix)) {
      unmatched <- keys[is.na(ix)]
      utils::write.csv(data.frame(id.match=head(unmatched, 100L)),
        file.path(outdir, paste0('3_unmatched_', ref_prefix[j], '_first100.csv')),
        row.names=FALSE)
      stop(length(unmatched), ' canonical id.match values from ', file,
           ' absent in PBMC metadata; first example: ', unmatched[1],
           '. Do not force matching; inspect unmatched report.')
    }
    same_sample <- as.character(submeta$orig.ident) == sample_ids[ix]
    if (anyNA(same_sample) || !all(same_sample))
      stop('Sample labels disagree for ', sum(!same_sample, na.rm=TRUE),
           ' cells of ', file)
    for (field in c('severity', 'severity2')) {
      if (field %in% colnames(submeta) && field %in% colnames(all_pbmc)) {
        same_field <- as.character(submeta[[field]]) ==
                      as.character(all_pbmc[[field]][ix])
        if (anyNA(same_field) || !all(same_field))
          stop('Clinical label ', field, ' does not match PBMC for ', file)
      }
    }
    membership[ix, j] <- TRUE
    source_summary <- rbind(source_summary,
      data.frame(lineage=ref_names[j], source_cells=nrow(submeta),
                 unique_matched_cells=length(unique(ix)), missing_ids=0L))
    logmsg('  Matched ', nrow(submeta), ' of ', nrow(submeta), ' cells in correct sample')
    rm(submeta, keys, ix)
  }
  utils::write.csv(source_summary, file.path(outdir,'3_subset_match_audit.csv'), row.names=FALSE)

  source_count <- rowSums(membership)
  assigned <- rep('Unassigned PBMC', nrow(all_pbmc))
  exactly_one <- which(source_count == 1L)
  assigned[exactly_one] <- ref_names[max.col(membership[exactly_one, , drop=FALSE],
                                             ties.method='first')]
  assigned[source_count > 1L] <- 'Conflicting membership'
  all_categories <- c(ref_names, 'Unassigned PBMC','Conflicting membership')
  logmsg('Full PBMC source lineage counts: ',
         paste(names(table(assigned)), as.integer(table(assigned)), collapse='; '))
  logmsg('Conflicting author-subset memberships: ', sum(source_count > 1L))

  # All-cell composition uses COMPLETE metadata, unlike Step 2's 40k UMAP.
  full_sample_tab <- as.data.frame(table(
    sample_id=sample_ids, source_lineage=factor(assigned, levels=all_categories)))
  names(full_sample_tab)[3] <- 'n_cells'
  full_sample_tab$sample_total <- as.integer(table(sample_ids)[full_sample_tab$sample_id])
  full_sample_tab$proportion <- full_sample_tab$n_cells / full_sample_tab$sample_total
  full_sample_tab$donor_inferred <- sub('_.*$', '', sub('^Exp_', '', full_sample_tab$sample_id))
  severity_by_sample <- tapply(as.character(all_pbmc$severity), sample_ids,
                               function(z) paste(unique(z), collapse='|'))
  severity2_by_sample <- tapply(as.character(all_pbmc$severity2), sample_ids,
                                function(z) paste(unique(z), collapse='|'))
  full_sample_tab$severity <- unname(severity_by_sample[full_sample_tab$sample_id])
  full_sample_tab$severity2 <- unname(severity2_by_sample[full_sample_tab$sample_id])
  full_sample_tab$severity_broad <- sub('_.*$', '', full_sample_tab$severity)
  if (any(grepl('|', full_sample_tab$severity, fixed=TRUE)) ||
      any(grepl('|', full_sample_tab$severity2, fixed=TRUE)))
    stop('Inconsistent clinical group labels within at least one sample')
  utils::write.csv(full_sample_tab, file.path(outdir,'3_full_cohort_lineage_by_sample.csv'),
                   row.names=FALSE)
  logmsg('Full cohort composition summarized across ', length(unique(sample_ids)), ' collections.')

  # Use source metadata as author-assigned reference membership, not predictions.
  # Step 2 kept the original PBMC id.match in Seurat cell metadata.
  # Seurat column names may use a different punctuation/suffix convention.
  sampled_meta <- obj@meta.data[colnames(obj), , drop=FALSE]
  sampled_keys <- checked_id_match(sampled_meta, 'Step 2 Seurat metadata')
  ii <- match(sampled_keys, pbmc_keys)
  if (anyNA(ii) || anyDuplicated(ii))
    stop('Sampled checkpoint id.match values do not match full PBMC metadata; ',
         sum(is.na(ii)), ' unmatched')
  if (!identical(as.character(obj$sample_id), sample_ids[ii]))
    stop('A Step 2 sampled cell has a conflicting sample label')
  obj$source_lineage <- assigned[ii]
  obj$source_membership_count <- as.integer(source_count[ii])
  logmsg('Mapped author-source identities to all ', ncol(obj), ' sampled UMAP cells.')

  # Quantify how heterogeneous each new cluster is in the author source view.
  clusters <- as.character(obj$analysis_clusters)
  ctab <- as.data.frame(table(cluster=clusters,
                              source_lineage=factor(obj$source_lineage,
                                                    levels=all_categories)))
  names(ctab)[3] <- 'n_cells'
  ctab$cluster_n <- as.integer(table(clusters)[ctab$cluster])
  ctab$fraction_within_cluster <- ctab$n_cells / ctab$cluster_n
  utils::write.csv(ctab, file.path(outdir,'3_cluster_lineage_composition.csv'), row.names=FALSE)
  by_cluster <- split(ctab, ctab$cluster)
  cluster_calls <- do.call(rbind, lapply(by_cluster, function(t) {
    t <- t[order(-t$fraction_within_cluster, t$source_lineage), ]
    best <- t$source_lineage[1]
    best_frac <- t$fraction_within_cluster[1]
    second_frac <- if (nrow(t) >= 2) t$fraction_within_cluster[2] else 0
    informative <- !(best %in% c('Unassigned PBMC','Conflicting membership'))
    confident <- informative && best_frac >= 0.70 && (best_frac-second_frac) >= 0.20
    data.frame(cluster=t$cluster[1], n_cells=t$cluster_n[1],
               top_source_label=as.character(best), top_fraction=best_frac,
               second_fraction=second_frac,
               cluster_report_label=if (confident) as.character(best) else 'Mixed/unresolved',
               passes_descriptive_rule=confident, row.names=NULL)
  }))
  cluster_calls <- cluster_calls[order(suppressWarnings(as.numeric(cluster_calls$cluster))), ]
  utils::write.csv(cluster_calls, file.path(outdir,'3_cluster_annotation_review.csv'),row.names=FALSE)
  cm <- setNames(cluster_calls$cluster_report_label,cluster_calls$cluster)
  obj$cluster_review_label <- unname(cm[clusters])
  if (anyNA(obj$cluster_review_label)) stop('Missing cluster labels from review mapping')
  logmsg('Clusters with source-dominant lineage: ',
         sum(cluster_calls$passes_descriptive_rule), '/', nrow(cluster_calls),
         '; others stay Mixed/unresolved.')

  # Check independent *marker evidence* using selected genes present in this
  # limited 1,800-gene processed panel. Negative values are allowed.
  requested <- c('CD3D','CD3E','TRAC','CD4','CD8A','CD8B','IL7R',
                 'MS4A1','CD79A','CD79B','CD74','JCHAIN','MZB1',
                 'NKG7','GNLY','KLRD1','PRF1','GZMB',
                 'LYZ','S100A8','S100A9','FCGR3A',
                 'LILRA4','FCER1A','PPBP','GP1BB','HBB')
  found <- requested[requested %in% rownames(obj)]
  utils::write.csv(data.frame(marker=requested, present=requested %in% found),
                   file.path(outdir,'3_marker_availability.csv'), row.names=FALSE)
  if (length(found) < 3L) stop('Too few expected markers in Seurat data layer')
  expression_data <- SeuratObject::LayerData(obj[['RNA']], layer='data')
  marker_expr <- expression_data[found, , drop=FALSE]
  rm(expression_data)
  gene_means <- lapply(split(seq_len(ncol(obj)), clusters), function(cols) {
    as.numeric(Matrix::rowMeans(marker_expr[, cols, drop=FALSE]))
  })
  means_df <- do.call(rbind, lapply(names(gene_means), function(cl) {
    data.frame(cluster=cl, marker=found,
               mean_processed_value=gene_means[[cl]])
  }))
  means_df$gene_relative_cluster_z <- NA_real_
  for (gene in found) {
    idx <- which(means_df$marker == gene)
    zz <- means_df$mean_processed_value[idx]
    sdv <- stats::sd(zz)
    means_df$gene_relative_cluster_z[idx] <- if (is.finite(sdv) && sdv > 0)
      (zz-mean(zz))/sdv else 0
  }
  utils::write.csv(means_df, file.path(outdir,'3_marker_mean_by_cluster.csv'),
                   row.names=FALSE)
  rm(marker_expr, gene_means)
  invisible(gc(verbose=FALSE))

  # Save enriched Step 3 object BEFORE figures to protect against plot failures.
  obj@misc$gse_annotation_method <- paste(
    'GSE196198 author-provided T/B/NK/Monocyte subset membership,',
    'cluster majority reference rule (>=70% and gap>=20%),',
    'manual verification required; Step 2 sample.')
  saved_rds <- file.path(outdir, '3_pbmc_seurat_annotated.rds')
  logmsg('Saving local annotated Seurat object...')
  saveRDS(obj, saved_rds, compress='gzip')
  logmsg('Checkpoint saved: ', saved_rds)

  # Transparent, readable display of *actual sampled UMAP* points. Avoid the
  # excessively faint DimPlot rasterization seen in Step 2 output.
  coords <- as.data.frame(Seurat::Embeddings(obj, reduction='umap'))
  if (ncol(coords) < 2) stop('UMAP coordinates unavailable')
  colnames(coords)[1:2] <- c('UMAP1','UMAP2')
  coords$source_lineage <- factor(obj$source_lineage, levels=all_categories)
  coords$cluster_review_label <- obj$cluster_review_label
  coords$analysis_clusters <- obj$analysis_clusters
  coords$severity_broad <- obj$severity_broad
  coords$sample_id <- obj$sample_id
  coords$donor_inferred <- sub('_.*$', '', sub('^Exp_', '', coords$sample_id))
  utils::write.csv(head(coords, 100L), file.path(outdir, '3_umap_coordinate_preview_100.csv'),
                   row.names=FALSE)

  library(ggplot2)
  save_plot <- function(plot, filename, width=12, height=8) {
    ggplot2::ggsave(file.path(outdir,filename), plot, width=width,height=height,
                    units='in',dpi=170,limitsize=FALSE,bg='white')
    logmsg('Figure: ', filename)
  }
  plot_umap <- function(group, title, filename) {
    df <- coords
    df$group <- factor(as.character(df[[group]]))
    p <- ggplot(df, aes(x=UMAP1, y=UMAP2, color=group)) +
      geom_point(size=0.18, alpha=0.75, shape=16, stroke=0) +
      scale_color_discrete(drop=FALSE, guide=guide_legend(override.aes=list(size=3,alpha=1))) +
      coord_equal() + theme_classic(base_size=12) +
      labs(title=title, x='UMAP 1',y='UMAP 2',color=group) +
      theme(plot.title=element_text(face='bold'),
            legend.position='right')
    save_plot(p, filename, 12,8)
  }
  plot_umap('source_lineage', 'PBMC UMAP | Author-subset lineage membership',
            '3_umap_author_lineages.png')
  plot_umap('cluster_review_label', 'PBMC UMAP | Cluster-level tentative labels',
            '3_umap_reviewed_cluster_labels.png')
  plot_umap('analysis_clusters', 'PBMC UMAP | 27 de novo Seurat clusters (readable)',
            '3_umap_clusters_readable.png')
  plot_umap('severity_broad', 'PBMC UMAP | Clinical severity (readable)',
            '3_umap_severity_readable.png')

  p_ct <- ggplot(ctab, aes(x=factor(cluster, levels=sort(unique(clusters),
                    decreasing=FALSE)), y=source_lineage,
                    fill=100*fraction_within_cluster)) +
          geom_tile(color='white', linewidth=0.25) +
          scale_fill_gradient(low='white', high='#134a8e',limits=c(0,100)) +
          labs(title='Author-subset membership within each new cluster',
               subtitle='Percentage of cells per cluster; based on 40,000-cell UMAP sample',
               x='New Seurat cluster',y='Author-source lineage',fill='% of cluster') +
          theme_minimal(base_size=11) +
          theme(axis.text.x=element_text(angle=70,hjust=1),panel.grid=element_blank())
  save_plot(p_ct,'3_cluster_source_concordance_heatmap.png', 13,5)

  if (length(found)) {
    # Gene-wise standardized MEANS across clusters (not UMI counts or logFC).
    means_df$cluster_f <- factor(means_df$cluster,
      levels=as.character(sort(unique(suppressWarnings(as.integer(means_df$cluster))))))
    means_df$marker_f <- factor(means_df$marker, levels=rev(found))
    p_marker <- ggplot(means_df, aes(x=cluster_f,y=marker_f,
                                     fill=gene_relative_cluster_z)) +
      geom_tile(color='white', linewidth=0.1) +
      scale_fill_gradient2(low='#2364aa',mid='white',high='#a91e28',
                           midpoint=0,limits=c(-3,3),oob=scales::squish) +
      labs(title='Selected immune markers across Seurat clusters',
           subtitle='Per-gene Z-score of CLUSTER MEAN signed processed values; descriptive only',
           x='New Seurat cluster', y='Gene', fill='Relative\nmean Z') +
      theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=70,hjust=1),panel.grid=element_blank())
    save_plot(p_marker,'3_cluster_marker_mean_heatmap.png', 14,8)
  }

  # Full-cohort bar chart: all PBMC cells, but all source labels come from the
  # original source subsets; this is not independent CITE/TCR label validation.
  p_comp <- ggplot(full_sample_tab,
      aes(x=sample_id,y=proportion,fill=source_lineage)) +
      geom_col(width=0.83) +
      facet_grid(~severity_broad, scales='free_x', space='free_x') +
      scale_y_continuous(labels=function(z) paste0(round(100*z),'%')) +
      labs(title='PBMC lineage composition of all 25 original collections',
           subtitle='399,044 original PBMC cells; author-subset membership (not 40k UMAP sample)',
           x='Collection',y='Fraction of PBMC cells',fill='Author-source lineage') +
      theme_minimal(base_size=11) +
      theme(axis.text.x=element_text(angle=75,hjust=1),
            panel.grid.major.x=element_blank())
  save_plot(p_comp,'3_full_cohort_sample_composition.png',17,7)

  # Descriptive per-patient proportions, not an inferential COVID-severity test.
  donor_tab <- aggregate(n_cells ~ donor_inferred + severity_broad + source_lineage,
                         data=full_sample_tab,FUN=sum)
  donor_den <- aggregate(n_cells ~ donor_inferred, data=donor_tab, FUN=sum)
  names(donor_den)[2] <- 'donor_total_cells'
  donor_tab <- merge(donor_tab,donor_den,by='donor_inferred',sort=FALSE)
  donor_tab$fraction_of_all_visits <- donor_tab$n_cells/donor_tab$donor_total_cells
  utils::write.csv(donor_tab,file.path(outdir,'3_donor_descriptive_all_visits.csv'),row.names=FALSE)

  writeLines(c(
    'GSE196198 Kim et al. STEP 03: annotation method and cautions',
    paste('Input PBMC cells (full metadata):', nrow(all_pbmc)),
    paste('Existing Step 2 sampled UMAP cells:', ncol(obj)),
    paste('New analysis clusters:', length(unique(clusters))),
    'Cell-level source_lineage is transferred by EXACT canonical id.match keys from',
    'the original T_cell, B_cell, NK_cell and Monocyte metadata, which are',
    'author-processed subsets of the SAME PBMC cohort.',
    'This is a reference/source-membership label, NOT an independent unsupervised',
    'or SingleR-based annotation of every cell.',
    'Cells not in any source subset remain Unassigned PBMC. Membership in two',
    'or more sources is explicitly labeled Conflicting membership.',
    'Cluster labels are assigned only if a known source category has >=70% of',
    'sampled cells AND at least a 20 percentage-point lead over runner-up.',
    'Otherwise cluster_report_label is Mixed/unresolved. Labels are tentative.',
    'The gene heatmap uses signed PREPROCESSED GEO RNA values; each gene is',
    'standardized across new cluster means, NOT count normalization.',
    'Do not interpret negative entries as negative RNA molecules.',
    'No new NormalizeData, SCTransform, RPCA, DE, pseudobulk, or independent',
    'SingleR annotation was run.',
    'Full-cohort composition is based on original metadata of all 399,044 cells;',
    'UMAP/cluster purity is based only on the proportional 40,000-cell sample.',
    '25 visits are not 25 independent individuals. Severity comparisons need',
    'patient-level design and appropriate inferential methods.',
    'Multiomics Figure S2D is NOT reproduced here, and Figure S2A-C are',
    'exploratory analogous visualizations, NOT pixel-exact reproductions.'
  ),file.path(outdir,'3_methods_and_cautions.txt'))
  writeLines(capture.output(sessionInfo()),file.path(outdir,'3_session_info.txt'))
  zip_bin <- Sys.which('zip')
  if (nzchar(zip_bin)) {
    archive <- file.path(datadir,'3_results.zip')
    if (file.exists(archive)) unlink(archive)
    file_names <- list.files(outdir,full.names=TRUE,recursive=FALSE)
    file_names <- file_names[!dir.exists(file_names) &
               basename(file_names)!='3_pbmc_seurat_annotated.rds']
    relative_files <- file.path('3_results',basename(file_names))
    oldwd <- getwd(); setwd(datadir); on.exit(setwd(oldwd),add=TRUE)
    ret <- suppressWarnings(system2(zip_bin, args=c('-q','-9',shQuote(archive),
                                                   shQuote(relative_files))))
    if (ret != 0 || !file.exists(archive))
      logmsg('WARNING: ZIP could not be made; compress 3_results manually, excluding RDS')
    else logmsg('Created upload ZIP: ', archive)
  } else logmsg('WARNING: zip executable not installed; manually zip results excluding RDS')
  logmsg('SUCCESS | Step 3 complete. Upload 3_results.zip; keep annotated RDS locally.')
}

tryCatch(main(), error=function(e) {
  msg <- paste0('STEP 03 FAILED: ',conditionMessage(e))
  cat(msg,'\n',file=stderr())
  try({
    args <- commandArgs(trailingOnly=TRUE)
    where <- if (length(args)) args[1] else '.'
    od <- file.path(where,'3_results')
    dir.create(od,showWarnings=FALSE,recursive=TRUE)
    writeLines(msg,file.path(od,'3_error.txt'))
  },silent=TRUE)
  quit(save='no',status=1)
})
