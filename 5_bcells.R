#!/usr/bin/env Rscript
# Kim et al. (2023) GSE196198 | STEP 05: B-cell-specific exploratory Seurat atlas
# Target: initial RNA-only reconstruction toward supplementary Figure S7A/B.
# Assumptions: author-processed signed, integrated B-cell RNA values; NOT raw UMI counts.
# Run from the directory containing the GSE196198 files:
# env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE conda run --no-capture-output -n ngs_seurat Rscript --vanilla 5_bcells.R .
# Optional: GSE_BCELL_MAX_CELLS=20000 (default = 40000), or GSE_BCELL_REUSE=1 to only remake figures.

options(stringsAsFactors = FALSE, scipen = 99)
Sys.setenv(OMP_NUM_THREADS = '10', OPENBLAS_NUM_THREADS = '10',
           MKL_NUM_THREADS = '10', VECLIB_MAXIMUM_THREADS = '10')

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  datadir <- normalizePath(if (length(args)) args[1] else '.', mustWork = TRUE)
  setwd(datadir)
  outdir <- file.path(datadir, '5_results')
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(outdir, '5_bcells_run_log.txt')
  if (file.exists(log_path)) unlink(log_path)
  logmsg <- function(...) {
    msg <- paste0(..., collapse = '')
    cat(msg, '\n', sep = '')
    cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), msg, '\n', file = log_path, append = TRUE)
  }

  logmsg('STEP 05 | B-cell-specific Seurat processed-RNA atlas | Kim et al. GSE196198')
  logmsg('Data directory: ', datadir)
  if (!requireNamespace('Seurat', quietly = TRUE) ||
      !requireNamespace('SeuratObject', quietly = TRUE) ||
      !requireNamespace('Matrix', quietly = TRUE) ||
      !requireNamespace('ggplot2', quietly = TRUE)) {
    stop('Seurat, SeuratObject, Matrix and ggplot2 must be installed in the active isolated R environment')
  }
  suppressPackageStartupMessages(library(Seurat))
  logmsg('R: ', R.version.string, ' | Seurat: ', as.character(packageVersion('Seurat')))
  logmsg('INPUT LAYER = previously integrated signed RNA; never raw UMI counts')

  cap <- suppressWarnings(as.integer(Sys.getenv('GSE_BCELL_MAX_CELLS', '40000')))
  if (is.na(cap) || cap < 1000L) stop('GSE_BCELL_MAX_CELLS must be >= 1000')
  seed <- 196198L
  set.seed(seed)
  checkpoint <- file.path(outdir, '5_bcell_seurat_processed_sampled.rds')
  reuse <- identical(Sys.getenv('GSE_BCELL_REUSE', '0'), '1')

  if (reuse) {
    if (!file.exists(checkpoint)) stop('Cannot reuse because local checkpoint is missing: ', checkpoint)
    logmsg('Recovery mode: using existing saved checkpoint')
    obj <- readRDS(checkpoint)
    if (!inherits(obj, 'Seurat') || !'umap' %in% names(obj@reductions) ||
        !'analysis_cluster' %in% colnames(obj@meta.data)) {
      stop('The checkpoint does not contain completed B-cell UMAP and clusters')
    }
    if (!identical(as.integer(obj@misc$sample_cap), cap) ||
        !identical(as.integer(obj@misc$seed), seed)) stop('Checkpoint parameters differ')
  } else {
    metaf <- file.path(datadir, 'GSE196198_B_cell_meta.data.txt.gz')
    exprf <- file.path(datadir, 'GSE196198_B_cell_normalized_count_integrated.txt.gz')
    if (!file.exists(metaf) || !file.exists(exprf)) {
      stop('Missing author B-cell metadata or expression matrix: ', metaf, ' or ', exprf)
    }

    logmsg('Reading all B-cell metadata...')
    mcon <- gzfile(metaf, open = 'rt')
    meta <- tryCatch(read.table(mcon, header = TRUE, row.names = 1, sep = '',
                                quote = '"', comment.char = '', check.names = FALSE,
                                stringsAsFactors = FALSE), finally = close(mcon))
    required <- c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'percent.mt',
                  'id.match', 'severity', 'severity2')
    missing <- setdiff(required, colnames(meta))
    if (length(missing)) stop('B-cell metadata missing: ', paste(missing, collapse = ', '))
    if (anyDuplicated(rownames(meta)) || anyDuplicated(meta$id.match) ||
        anyNA(meta$id.match)) stop('B-cell metadata identifiers are not one-to-one')
    n_total <- nrow(meta)
    logmsg('B-cell records: ', n_total, ' | sample labels: ', length(unique(meta$orig.ident)))
    for (field in c('nCount_RNA', 'nFeature_RNA', 'percent.mt')) {
      meta[[field]] <- suppressWarnings(as.numeric(meta[[field]]))
      if (any(!is.finite(meta[[field]]))) stop('Invalid numeric QC field: ', field)
    }
    qc_pass <- meta$nFeature_RNA > 200 & meta$nFeature_RNA < 3500 &
               meta$nCount_RNA < 15000 & meta$percent.mt < 15
    logmsg('Supplied B cells meeting published QC thresholds: ', sum(qc_pass), '/', n_total)

    # The GEO export has a whitespace-delimited, exceptionally wide header.
    logmsg('Reading expression header and validating all B-cell column identifiers...')
    econ <- gzfile(exprf, 'rt')
    header <- readLines(econ, n = 1L, warn = FALSE)
    if (!length(header)) stop('Empty expression header')
    fmt <- if (grepl('\t', header, fixed = TRUE)) 'tab' else
           if (grepl(',', header, fixed = TRUE)) 'comma' else
           if (grepl('[[:space:]]', header)) 'whitespace' else 'unknown'
    if (fmt == 'unknown') stop('Cannot determine expression separator; inspect original file')
    parse_fields <- function(line) {
      line <- sub('\r$', '', line)
      if (fmt == 'tab') return(strsplit(line, '\t', fixed = TRUE)[[1]])
      if (fmt == 'comma') return(strsplit(line, ',', fixed = TRUE)[[1]])
      strsplit(trimws(line), '[[:space:]]+')[[1]]
    }
    stripq <- function(z) sub('"$', '', sub('^"', '', z))
    idnorm <- function(z) gsub('[^A-Za-z0-9_.]', '.', z)
    heads <- stripq(parse_fields(header))
    heads[1] <- sub('^\ufeff', '', heads[1])
    if (length(heads) == n_total + 1L) heads <- heads[-1L]
    if (length(heads) != n_total) stop('Expression header fields ', length(heads),
                                       ', expected ', n_total)
    meta_ids <- idnorm(rownames(meta))
    expr_ids <- idnorm(heads)
    if (anyDuplicated(meta_ids) || anyDuplicated(expr_ids)) {
      stop('A cell-ID collision was created by punctuation normalization')
    }
    ordered <- match(expr_ids, meta_ids)
    if (anyNA(ordered) || anyDuplicated(ordered)) {
      stop('Expression cell IDs do not map one-to-one to B-cell metadata. STOP: no positional matching')
    }
    meta <- meta[ordered, , drop = FALSE]
    qc_pass <- qc_pass[ordered]
    logmsg('Full matrix-header match: ', length(ordered), '/', n_total, '; separator = ', fmt)

    # Optionally cross-check canonical IDs against the Step 3 linkage audit.
    step3_audit <- file.path(datadir, '3_results', '3_subset_match_audit.csv')
    if (file.exists(step3_audit)) {
      previous <- read.csv(step3_audit)
      k <- which(previous$lineage == 'B cells')
      if (length(k) == 1L && previous$unique_matched_cells[k] != n_total) {
        stop('B-cell metadata size is inconsistent with successful Step 3 matching')
      }
    }

    eligible <- which(qc_pass)
    if (!length(eligible)) stop('No QC-eligible B cells')
    n_take <- min(cap, length(eligible))
    groups <- split(eligible, as.character(meta$orig.ident[eligible]))
    sizes <- lengths(groups)
    rawq <- sizes * n_take / sum(sizes)
    quotas <- floor(rawq)
    extra <- n_take - sum(quotas)
    if (extra > 0L) {
      idx <- order(rawq - quotas, decreasing = TRUE)[seq_len(extra)]
      quotas[idx] <- quotas[idx] + 1L
    }
    sampled <- sort(unlist(Map(function(index, quota) {
      if (quota == 0L) integer(0) else sample(index, quota, replace = FALSE)
    }, groups, quotas), use.names = FALSE))
    if (length(sampled) != n_take || anyDuplicated(sampled)) stop('Sampling integrity failure')
    sampling <- data.frame(sample_id = names(groups), original_B_cells = as.integer(sizes),
                           sampled_B_cells = as.integer(quotas),
                           fraction = as.numeric(quotas) / as.numeric(sizes))
    write.csv(sampling, file.path(outdir, '5_sampling_by_sample.csv'), row.names = FALSE)
    logmsg('Sample-balanced proportional sampling: ', n_take, '/', n_total,
           ' cells, seed ', seed)

    selmeta <- meta[sampled, , drop = FALSE]
    rownames(selmeta) <- heads[sampled]
    selmeta$sample_id <- as.character(selmeta$orig.ident)
    selmeta$donor_inferred <- sub('_.*$', '', sub('^Exp_', '', selmeta$sample_id))
    selmeta$visit_inferred <- ifelse(grepl('_', sub('^Exp_', '', selmeta$sample_id)),
                                    sub('^[^_]+_', '', sub('^Exp_', '', selmeta$sample_id)), '')
    selmeta$severity_broad <- sub('_.*$', '', as.character(selmeta$severity))
    selmeta$author_cluster <- if ('seurat_clusters' %in% names(selmeta))
                                  as.character(selmeta$seurat_clusters) else NA_character_

    first_gene <- readLines(econ, n = 1L, warn = FALSE)
    if (!length(first_gene)) stop('Expression matrix has no gene rows')
    pre <- parse_fields(first_gene)
    if (length(pre) != n_total + 1L) stop('First gene row has wrong field count')
    logmsg('First gene: ', stripq(pre[1]), ' | values: ', length(pre) - 1L)
    # We consumed the first gene row; reset only the connection and stream data once.
    close(econ)
    econ <- gzfile(exprf, 'rt')
    readLines(econ, n = 1L, warn = FALSE)

    n_features_expected <- 2000L  # previously verified author B-subset export
    logmsg('Streaming ', n_features_expected, ' gene rows; storing selected cells only')
    logmsg('Estimated dense double expression matrix = ',
           sprintf('%.2f', n_features_expected * n_take * 8 / 1024^3),
           ' GiB; Seurat requires additional working memory')
    emat <- matrix(NA_real_, nrow = n_features_expected, ncol = n_take)
    genes <- character(n_features_expected)
    seen <- 0L
    repeat {
      ln <- readLines(econ, n = 1L, warn = FALSE)
      if (!length(ln)) break
      seen <- seen + 1L
      if (seen > n_features_expected) stop('More than 2,000 B-cell genes found')
      fields <- parse_fields(ln)
      if (length(fields) != n_total + 1L) {
        stop('Gene row ', seen, ' has ', length(fields), ' fields; expected ', n_total + 1L)
      }
      genes[seen] <- stripq(fields[1])
      vals <- suppressWarnings(as.numeric(stripq(fields[sampled + 1L])))
      if (length(vals) != n_take || any(!is.finite(vals))) {
        stop('Nonfinite/non-numeric expression values for gene ', genes[seen])
      }
      emat[seen, ] <- vals
      if (seen %% 100L == 0L) logmsg('Expression genes read: ', seen, '/', n_features_expected)
    }
    close(econ)
    if (seen != n_features_expected || anyDuplicated(genes) || any(!nzchar(genes))) {
      stop('Unexpected B-cell gene count, repeated IDs or empty gene names')
    }
    rownames(emat) <- gsub('_', '-', genes, fixed = TRUE)
    if (anyDuplicated(rownames(emat))) stop('Gene ID collision after Seurat convention')
    colnames(emat) <- heads[sampled]
    if (!identical(colnames(emat), rownames(selmeta))) stop('Metadata cell order mismatch')
    logmsg('Expression range (sample): ', sprintf('%.3f', min(emat)), ' to ', sprintf('%.3f', max(emat)))

    assay <- SeuratObject::CreateAssay5Object(data = emat)
    rm(emat)
    invisible(gc(verbose = FALSE))
    obj <- SeuratObject::CreateSeuratObject(counts = assay, assay = 'RNA',
                                            project = 'GSE196198_Bcells_processed',
                                            meta.data = selmeta)
    rm(assay, meta, selmeta, groups, fields)
    invisible(gc(verbose = FALSE))
    DefaultAssay(obj) <- 'RNA'
    logmsg('Created data-only Seurat v5 object. NO invented counts / NormalizeData / SCTransform.')
    features <- rownames(obj)
    obj <- Seurat::ScaleData(obj, features = features, do.center = TRUE,
                             do.scale = TRUE, scale.max = 10, verbose = FALSE)
    dims <- min(30L, nrow(obj) - 1L, ncol(obj) - 1L)
    obj <- Seurat::RunPCA(obj, features = features, npcs = dims,
                          seed.use = seed, verbose = FALSE)
    use_dims <- seq_len(min(30L, ncol(Seurat::Embeddings(obj, 'pca'))))
    obj <- Seurat::FindNeighbors(obj, reduction = 'pca', dims = use_dims, verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = 0.5, algorithm = 1,
                                random.seed = seed, verbose = FALSE)
    obj$analysis_cluster <- as.character(SeuratObject::Idents(obj))
    obj <- Seurat::RunUMAP(obj, reduction = 'pca', dims = use_dims,
                           n.neighbors = 30, min.dist = 0.3,
                           seed.use = seed, verbose = FALSE)
    obj@misc$seed <- seed
    obj@misc$sample_cap <- cap
    obj@misc$input_total_bcells <- n_total
    obj@misc$method <- 'exploratory processed-RNA B-cell reclustering; not original raw count RPCA'
    logmsg('New B-cell clusters: ', length(unique(obj$analysis_cluster)))
    logmsg('Saving local checkpoint BEFORE plotting: ', checkpoint)
    saveRDS(obj, checkpoint, compress = 'gzip')
    logmsg('Checkpoint saved; if plotting fails, rerun with GSE_BCELL_REUSE=1')
  }

  # All charts below are descriptions of the processed-data analysis, not inferential tests.
  md <- obj@meta.data
  if (!all(c('sample_id', 'donor_inferred', 'severity', 'severity2',
             'severity_broad', 'analysis_cluster') %in% names(md))) {
    stop('Missing required metadata columns in B-cell checkpoint')
  }
  emb <- as.data.frame(Seurat::Embeddings(obj, reduction = 'umap'))
  emb$analysis_cluster <- md[rownames(emb), 'analysis_cluster']
  emb$severity_broad <- md[rownames(emb), 'severity_broad']
  emb$severity2 <- md[rownames(emb), 'severity2']
  emb$sample_id <- md[rownames(emb), 'sample_id']
  emb$donor_inferred <- md[rownames(emb), 'donor_inferred']
  xy <- colnames(emb)[seq_len(2)]
  savepng <- function(plt, name, width = 12, height = 8) {
    ggplot2::ggsave(file.path(outdir, name), plt, width = width, height = height,
                    units = 'in', dpi = 170, limitsize = FALSE)
    logmsg('Figure: ', name)
  }
  plotUMAP <- function(group, name, title) {
    plt <- ggplot2::ggplot(emb, ggplot2::aes(x = .data[[xy[1]]], y = .data[[xy[2]]], color = .data[[group]])) +
      ggplot2::geom_point(size = 0.12, alpha = 0.68, stroke = 0) +
      ggplot2::coord_equal() + ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(x = 'UMAP 1', y = 'UMAP 2', color = group, title = title,
                    subtitle = 'New Seurat clustering of sampled author-processed B-cell RNA') +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
    savepng(plt, name, 13, 8)
  }
  plotUMAP('analysis_cluster', '5_bcell_umap_clusters.png', 'B-cell subclusters')
  plotUMAP('severity_broad', '5_bcell_umap_severity.png', 'B-cell subclusters by broad clinical group')
  plotUMAP('severity2', '5_bcell_umap_severity2.png', 'B-cell subclusters by supplied severity2')
  plotUMAP('donor_inferred', '5_bcell_umap_donor.png', 'B-cell subclusters by inferred donor')
  if ('author_cluster' %in% colnames(md) && any(!is.na(md$author_cluster))) {
    emb$author_cluster <- md[rownames(emb), 'author_cluster']
    plotUMAP('author_cluster', '5_bcell_umap_author_clusters.png',
             'Author-supplied B-cell clusters (not new Seurat clustering)')
  }

  clusters <- sort(unique(as.character(md$analysis_cluster)))
  by_sample <- as.data.frame(table(sample_id = md$sample_id, cluster = md$analysis_cluster))
  by_sample <- by_sample[by_sample$Freq > 0, ]
  names(by_sample)[3] <- 'n_cells'
  totals <- as.data.frame(table(sample_id = md$sample_id))
  names(totals)[2] <- 'sample_total'
  by_sample <- merge(by_sample, totals, by = 'sample_id')
  by_sample$fraction <- by_sample$n_cells / by_sample$sample_total
  samp_meta <- md[!duplicated(md$sample_id), c('sample_id', 'donor_inferred', 'severity',
                                             'severity2', 'severity_broad'), drop = FALSE]
  by_sample <- merge(by_sample, samp_meta, by = 'sample_id', all.x = TRUE)
  write.csv(by_sample, file.path(outdir, '5_clusters_by_sample_descriptive.csv'), row.names = FALSE)
  write.csv(as.data.frame(table(cluster = md$analysis_cluster)),
            file.path(outdir, '5_cluster_sizes.csv'), row.names = FALSE)
  sev <- as.data.frame(table(severity2 = md$severity2))
  write.csv(sev, file.path(outdir, '5_severity2_cell_counts.csv'), row.names = FALSE)
  pfrac <- ggplot2::ggplot(by_sample, ggplot2::aes(x = sample_id, y = fraction, fill = cluster)) +
    ggplot2::geom_col(width = 0.93) + ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1)) +
    ggplot2::labs(title = 'Sample-wise fractions of new B-cell clusters',
                  subtitle = '40,000-cell proportional sample; descriptive only',
                  x = 'Collection (25 labels)', y = 'Fraction of sampled B cells')
  savepng(pfrac, '5_bcell_cluster_composition_by_sample.png', 16, 8)

  # Check all requested paper-relevant markers; nothing is asserted if absent.
  marker_sets <- list(
    'B-cell lineage' = c('MS4A1', 'CD19', 'CD79A', 'CD79B', 'CD22', 'CD74'),
    'Naive B candidates' = c('TCL1A', 'IGHD', 'IGHM', 'FCER2', 'IL4R'),
    'Memory B candidates' = c('CD27', 'TNFRSF13B', 'BANK1', 'AIM2', 'CD80'),
    'Atypical/activated B candidates' = c('FCRL5', 'ITGAX', 'TBX21', 'FCRL4', 'CD86'),
    'Plasma/secretory candidates' = c('JCHAIN', 'MZB1', 'XBP1', 'PRDM1', 'SDC1', 'IGHG1'),
    'Proliferation/interferon candidates' = c('MKI67', 'TOP2A', 'ISG15', 'IFIT3', 'MX1')
  )
  available <- data.frame(category = rep(names(marker_sets), lengths(marker_sets)),
                          gene = unlist(marker_sets, use.names = FALSE))
  available$present <- available$gene %in% rownames(obj)
  write.csv(available, file.path(outdir, '5_marker_availability.csv'), row.names = FALSE)
  present <- unique(available$gene[available$present])
  logmsg('Available tested marker genes: ', length(present), '/', nrow(available))
  if (length(present) < 6) stop('Too few expected markers in RNA matrix; verify feature names')

  # Access only the signed authors' processed RNA data layer. No DE is attempted.
  expr <- SeuratObject::LayerData(obj, assay = 'RNA', layer = 'data')
  if (!identical(colnames(expr), colnames(obj))) {
    expr <- expr[, colnames(obj), drop = FALSE]
  }
  meanmat <- vapply(clusters, function(cl) {
    cells <- which(md$analysis_cluster == cl)
    Matrix::rowMeans(expr[present, cells, drop = FALSE])
  }, numeric(length(present)))
  dimnames(meanmat) <- list(present, clusters)
  rawmean <- data.frame(gene = rep(rownames(meanmat), ncol(meanmat)),
                        cluster = rep(colnames(meanmat), each = nrow(meanmat)),
                        processed_mean = as.vector(meanmat))
  rawmean$category <- available$category[match(rawmean$gene, available$gene)]
  write.csv(rawmean, file.path(outdir, '5_marker_mean_by_cluster.csv'), row.names = FALSE)
  z <- t(scale(t(meanmat)))
  z[!is.finite(z)] <- 0
  # Important: element-wise replacement PRESERVES matrix dimensions/dimnames.
  # pmax(-3, pmin(3, z)) can turn a matrix into a plain vector in R.
  z[z < -3] <- -3
  z[z > 3] <- 3
  if (!is.matrix(z) || !identical(dim(z), dim(meanmat))) {
    stop('Marker heatmap matrix dimensions were lost; refusing to mislabel data')
  }
  dimnames(z) <- dimnames(meanmat)
  zlong <- data.frame(gene = rep(rownames(z), times = ncol(z)),
                      cluster = rep(colnames(z), each = nrow(z)),
                      z_mean = as.vector(z))
  zlong$gene <- factor(zlong$gene, levels = rev(present))
  zlong$cluster <- factor(zlong$cluster, levels = clusters)
  zlong$category <- available$category[match(as.character(zlong$gene), available$gene)]
  pheat <- ggplot2::ggplot(zlong, ggplot2::aes(x = cluster, y = gene, fill = z_mean)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = '#315a99', mid = 'white',
                                                         high = '#b93335', midpoint = 0,
                                                         limits = c(-3, 3)) +
    ggplot2::facet_grid(category ~ ., scales = 'free_y', space = 'free_y') +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 60, hjust = 1)) +
    ggplot2::labs(title = 'RNA marker patterns across new B-cell clusters',
                  subtitle = 'Per-gene Z-score of cluster means from signed integrated RNA; no raw-count DE',
                  x = 'New Seurat cluster', y = 'Gene', fill = 'Mean Z')
  savepng(pheat, '5_bcell_marker_heatmap.png', 16, 10)

  showmarkers <- intersect(c('MS4A1', 'CD79A', 'TCL1A', 'IGHD', 'CD27',
                             'FCRL5', 'ITGAX', 'JCHAIN', 'MZB1', 'XBP1',
                             'MKI67', 'ISG15'), rownames(obj))
  if (length(showmarkers)) {
    tryCatch({
      pfeat <- Seurat::FeaturePlot(obj, features = head(showmarkers, 12),
                                   reduction = 'umap', ncol = 4L, pt.size = 0.06,
                                   raster = FALSE, order = FALSE)
      savepng(pfeat, '5_bcell_marker_featureplots.png', 18, 12)
    }, error = function(e) {
      logmsg('WARNING: marker FeaturePlot failed but checkpoint and main results are saved: ',
             conditionMessage(e))
      writeLines(conditionMessage(e), file.path(outdir, '5_optional_featureplot_error.txt'))
    })
  }
  rm(expr)
  invisible(gc(verbose = FALSE))

  writeLines(c(
    'Kim et al. 2023 GSE196198 | Step 05',
    'Target: preliminary B-cell RNA-only panels related to S7A/B.',
    'Input: GSE196198_B_cell_normalized_count_integrated.txt.gz plus matching metadata.',
    paste('Input B-cell metadata records:', obj@misc$input_total_bcells),
    paste('Selected cells:', ncol(obj), ' | supplied genes:', nrow(obj)),
    'Input values contain negative noninteger entries; these are author-processed integrated RNA, NOT raw UMI counts.',
    'Created a Seurat v5 assay with a DATA layer, without inventing a counts layer.',
    'No NormalizeData, SCTransform, FindVariableFeatures, RPCA, FindAllMarkers or inferential DE.',
    'Original paper annotates B cells with combined RNA, antibody and BCR evidence; we have not reproduced that full procedure.',
    'Naive/memory/atypical/plasma programs are descriptive candidate marker signatures; no definitive B-cell subtypes assigned.',
    'Figures describe a proportional random sample by sample; not an exact reproduction of author UMAP coordinates or complete-cohort figures.',
    'Clinical severe categories are preserved exactly as severity2; no unsupported mapping to usual/delayed severe.',
    'Repeated visits belong to the same donors. Cell counts/stacked bars are descriptive and do not prove severity associations.',
    'Somatic hypermutation, BCR lineage relationships and severity-associated expression require later dedicated analysis.',
    'Saved local RDS checkpoint is deliberately excluded from 5_results.zip.'
  ), file.path(outdir, '5_methods_and_limitations.txt'))
  writeLines(capture.output(sessionInfo()), file.path(outdir, '5_session_info.txt'))

  # Remove stale failure report from a previous unsuccessful attempt.
  previous_error <- file.path(outdir, '5_error.txt')
  if (file.exists(previous_error)) unlink(previous_error)

  zipfile <- file.path(datadir, '5_results.zip')
  all_files <- list.files(outdir, full.names = TRUE, recursive = FALSE)
  to_share <- all_files[!dir.exists(all_files) & basename(all_files) != basename(checkpoint)]
  if (file.exists(zipfile)) unlink(zipfile)
  zip_exec <- Sys.which('zip')
  if (nzchar(zip_exec)) {
    status <- suppressWarnings(system2(zip_exec,
                  args = c('-q', '-9', shQuote(zipfile), shQuote(to_share))))
    if (is.numeric(status) && status == 0L && file.exists(zipfile)) {
      logmsg('Created upload ZIP: ', zipfile)
    } else {
      logmsg('WARNING: automatic ZIP failed. Manually ZIP 5_results excluding local RDS')
    }
  } else logmsg('WARNING: zip program missing; archive 5_results manually excluding RDS')
  logmsg('SUCCESS: Step 05 complete. Send 5_results.zip; retain checkpoint locally.')
}

tryCatch(main(), error = function(e) {
  msg <- paste0('STEP 05 FAILED: ', conditionMessage(e))
  cat(msg, '\n', file = stderr())
  try({
    args <- commandArgs(trailingOnly = TRUE)
    p <- file.path(if (length(args)) args[[1]] else '.', '5_results')
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
    writeLines(msg, file.path(p, '5_error.txt'))
  }, silent = TRUE)
  quit(save = 'no', status = 1L)
})
