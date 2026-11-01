#!/usr/bin/env Rscript
# GSE196198 | Step 08: donor-aware descriptive gene-expression contrasts
# Inputs: LOCAL Step 4/5 Seurat checkpoints + Step 6 clinical index.
# IMPORTANT: These are author-INTEGRATED, SIGNED values, NOT raw counts.
# Do not call FindMarkers(), edgeR, DESeq2, SCTransform, or log2FC here.
# All gene-level tests below are EXPLORATORY permutations of donor-mean
# processed expression, NOT count-based differential expression.
# Run: env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
#      conda run --no-capture-output -n ngs_seurat Rscript --vanilla 8_gene_expression.R .

options(stringsAsFactors=FALSE, scipen=99)
Sys.setenv(OMP_NUM_THREADS='10', OPENBLAS_NUM_THREADS='10',
           MKL_NUM_THREADS='10', VECLIB_MAXIMUM_THREADS='10')

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  root <- normalizePath(if (length(args)) args[1] else '.', mustWork=TRUE)
  setwd(root)
  out <- file.path(root, '8_results')
  dir.create(out, showWarnings=FALSE, recursive=TRUE)
  logfile <- file.path(out, '8_run_log.txt')
  if (file.exists(logfile)) unlink(logfile)
  logmsg <- function(...) {
    x <- paste0(..., collapse='')
    cat(x, '\n', sep='')
    cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), x, '\n',
        file=logfile, append=TRUE)
  }
  require_pkg <- c('Seurat', 'SeuratObject', 'Matrix', 'ggplot2')
  for (p in require_pkg)
    if (!requireNamespace(p, quietly=TRUE)) stop('Missing R package: ', p)
  logmsg('STEP 08 | Processed RNA: donor-level descriptive gene contrasts')
  logmsg('R: ', R.version.string, ' | Seurat: ', as.character(packageVersion('Seurat')))
  logmsg('No raw count layer, no log fold changes, no count-based DE.')
  min_cells <- suppressWarnings(as.integer(Sys.getenv('GSE_GENE_MIN_CELLS', '20')))
  if (is.na(min_cells) || min_cells < 5L) stop('GSE_GENE_MIN_CELLS must be >=5')
  logmsg('Minimum sampled cells per donor/collection/context: ', min_cells)

  clinfile <- file.path(root,'6_results','6_sample_donor_clinical_index.csv')
  if (!file.exists(clinfile)) stop('Need completed Step 6 index: ', clinfile)
  clinical <- read.csv(clinfile, check.names=FALSE)
  required <- c('sample_id','donor_inferred','disease_group','severity2','visit_code')
  if (!all(required %in% names(clinical))) stop('Clinical index columns missing')
  if (anyNA(clinical[,required]) || anyDuplicated(clinical$sample_id) ||
      anyDuplicated(clinical[c('donor_inferred','visit_code')])) {
    stop('Clinical index has missing/duplicate identities')
  }
  if (nrow(clinical)!=25L || length(unique(clinical$donor_inferred))!=11L)
    stop('Clinical index unexpectedly differs from earlier audit (25 samples, 11 donors)')
  donor_groups <- split(clinical$disease_group, clinical$donor_inferred)
  if (any(vapply(donor_groups,function(z) length(unique(z))!=1L, logical(1))))
    stop('A donor has changing disease-group labels')
  initial <- do.call(rbind, lapply(split(clinical, clinical$donor_inferred),
    function(z) z[order(as.numeric(z$visit_code)),,drop=FALSE][1L,,drop=FALSE]))
  rownames(initial) <- NULL
  if (nrow(initial)!=11L || anyDuplicated(initial$donor_inferred))
    stop('Failed to identify exactly one first recorded collection per donor')
  logmsg('First-recorded collection distribution: ',
         paste(capture.output(print(table(initial$disease_group))),collapse=' '))

  # Contexts are post-hoc RNA clusters discovered in Steps 4/5.
  # Selecting them using the same RNA limits biological independence.
  contexts <- list(
    T=list(
      T_all=list(clusters=NULL,desc='All sampled T cells'),
      T_cytotoxic_2_10=list(clusters=c('2','10'),desc='Cytotoxic-like T clusters'),
      T_naive_0_7=list(clusters=c('0','7'),desc='Naive-like T clusters'),
      T_interferon_16=list(clusters=c('16'),desc='Interferon-like T cluster')
    ),
    B=list(
      B_all=list(clusters=NULL,desc='All sampled B cells'),
      B_secretory_13=list(clusters=c('13'),desc='Secretory/plasmablast-like B cluster'),
      B_atypical_6=list(clusters=c('6'),desc='Atypical-like B cluster'),
      B_interferon_14=list(clusters=c('14'),desc='Interferon-like B cluster')
    )
  )
  marker_panels <- list(
    T=c('CD3D','CD4','CD8A','CCR7','TCF7','NKG7','GNLY','GZMB',
        'PRF1','IFIT1','ISG15','MX1','MKI67','CD38'),
    B=c('MS4A1','CD79A','TCL1A','IGHM','IGHD','JCHAIN','MZB1',
        'XBP1','PRDM1','FCRL5','ITGAX','IFIT3','ISG15','MX1')
  )
  all_contrasts <- list()
  all_coverage <- list()
  all_markers <- list()
  top_rows <- list()
  missing_markers <- list()
  calc_exact_p <- function(v, a, b) {
    # Rows = genes, columns = independent first-collection donors.
    na <- length(a); nb <- length(b)
    all_idx <- c(a,b); zz <- v[,all_idx,drop=FALSE]
    observed <- rowMeans(zz[,seq.int(na+1L, na+nb),drop=FALSE]) -
                rowMeans(zz[,seq_len(na),drop=FALSE])
    combinations <- combn(seq_len(na+nb), na)
    # Two-sided exact randomization over all allocations, no pseudoreplication.
    extremal <- integer(nrow(zz))
    for (j in seq_len(ncol(combinations))) {
      ai <- combinations[,j]
      bi <- setdiff(seq_len(na+nb), ai)
      delta <- rowMeans(zz[,bi,drop=FALSE]) - rowMeans(zz[,ai,drop=FALSE])
      extremal <- extremal + as.integer(abs(delta) >= abs(observed)-1e-10)
    }
    list(delta=observed, p=extremal / ncol(combinations),
         n_alloc=ncol(combinations))
  }

  for (lin in c('T','B')) {
    path <- if (lin=='T') file.path(root,'4_results','4_tcell_seurat_processed_sampled.rds')
            else file.path(root,'5_results','5_bcell_seurat_processed_sampled.rds')
    if (!file.exists(path)) stop('Local Seurat checkpoint missing: ',path,
      '\nIt was excluded from the uploaded ZIP. Do NOT rerun UMAP if file exists elsewhere.')
    logmsg('Reading ',lin,' checkpoint: ',path)
    obj <- readRDS(path)
    if (!inherits(obj,'Seurat')) stop('Not a Seurat object: ',path)
    if (!'RNA' %in% names(obj@assays)) stop('Missing RNA assay in ',lin)
    if (!'data' %in% SeuratObject::Layers(obj[['RNA']]))
      stop('Integrated data layer missing in ',lin,'; never substitute counts')
    md <- obj@meta.data
    req <- c('sample_id','donor_inferred','severity2','analysis_cluster')
    if (!all(req %in% names(md))) stop('Checkpoint metadata missing: ',lin)
    if (!setequal(unique(as.character(md$sample_id)), clinical$sample_id))
      stop(lin,' checkpoint sample inventory differs from Step 6')
    samples <- clinical$sample_id
    for (s in samples) {
      ix <- which(md$sample_id==s)
      if (!length(ix)) stop('Absent sample ',s)
      ci <- clinical[match(s,clinical$sample_id),]
      if (any(as.character(md$donor_inferred[ix])!=as.character(ci$donor_inferred)) ||
          any(as.character(md$severity2[ix])!=as.character(ci$severity2)))
        stop('Clinical mismatch in ',lin,' sample ',s)
    }
    logmsg(lin,': ',ncol(obj),' sampled cells and ',nrow(obj),' processed genes')
    mat <- SeuratObject::LayerData(obj,assay='RNA',layer='data')
    if (!setequal(colnames(mat),rownames(md))) stop('Expression/metadata cells differ in ',lin)
    if (!identical(colnames(mat),rownames(md))) md <- md[colnames(mat),,drop=FALSE]
    genes <- rownames(mat)
    if (anyNA(genes) || anyDuplicated(genes)) stop('Invalid gene names in ',lin)
    present_markers <- intersect(marker_panels[[lin]], genes)
    absent_markers <- setdiff(marker_panels[[lin]],genes)
    missing_markers[[lin]] <- data.frame(lineage=lin,gene=marker_panels[[lin]],
                              available=marker_panels[[lin]] %in% genes)
    logmsg(lin,': ',length(present_markers),'/',length(marker_panels[[lin]]),
           ' prespecified genes available')
    for (ctx in names(contexts[[lin]])) {
      spec <- contexts[[lin]][[ctx]]
      chosen <- if (is.null(spec$clusters)) rep(TRUE,nrow(md)) else
        as.character(md$analysis_cluster) %in% spec$clusters
      samplemeans <- matrix(NA_real_,nrow=length(genes),ncol=length(samples),
                            dimnames=list(genes,samples))
      ncs <- integer(length(samples))
      for (i in seq_along(samples)) {
        cell_idx <- which(chosen & md$sample_id==samples[i])
        ncs[i] <- length(cell_idx)
        if (ncs[i] >= min_cells) {
          x <- Matrix::rowMeans(mat[,cell_idx,drop=FALSE])
          samplemeans[,i] <- as.numeric(x)
        }
      }
      nflag <- ifelse(ncs >= min_cells, 'retained', 'below_min_cells')
      cov <- data.frame(lineage=lin,context=ctx,description=spec$desc,
                        sample_id=samples,donor=clinical$donor_inferred,
                        disease_group=clinical$disease_group,
                        severity2=clinical$severity2,visit_code=clinical$visit_code,
                        sampled_context_cells=ncs,min_cells=min_cells,
                        used_for_expression=nflag,stringsAsFactors=FALSE)
      all_coverage[[ctx]] <- cov
      logmsg(ctx,': ',sum(ncs>=min_cells),' / ',length(samples),
             ' collections with >=',min_cells,' cells')
      first_here <- initial$sample_id
      ix <- match(first_here, samples)
      eligible <- which(!is.na(samplemeans[1L,ix]))
      base <- initial[eligible,,drop=FALSE]
      if (!nrow(base)) logmsg('    No baseline donors with enough cells in ',ctx)
      mild <- which(base$disease_group=='Mild')
      sev <- which(base$disease_group %in% c('Severe.m','Severe.s'))
      n_mild <- length(mild); n_sev <- length(sev)
      # For reproducibility, no partial/masked gene values are allowed.
      baseline_matrix <- samplemeans[,ix[eligible],drop=FALSE]
      if (any(!is.finite(baseline_matrix))) stop('Incomplete baseline gene means for ',ctx)
      contrast <- data.frame(lineage=lin, context=ctx, gene=genes,
                             n_mild_donors=n_mild,n_severe_donors=n_sev,
                             mean_mild=NA_real_,mean_severe=NA_real_,
                             severe_minus_mild_processed_units=NA_real_,
                             median_mild=NA_real_,median_severe=NA_real_,
                             exploratory_exact_p=NA_real_,
                             BH_q_within_context=NA_real_,
                             permutation_allocations=NA_integer_,
                             contrast_status='insufficient_donors',
                             stringsAsFactors=FALSE)
      if (n_mild>=2L && n_sev>=3L) {
        contrast$mean_mild <- rowMeans(baseline_matrix[,mild,drop=FALSE])
        contrast$mean_severe <- rowMeans(baseline_matrix[,sev,drop=FALSE])
        contrast$median_mild <- apply(baseline_matrix[,mild,drop=FALSE],1,median)
        contrast$median_severe <- apply(baseline_matrix[,sev,drop=FALSE],1,median)
        ex <- calc_exact_p(baseline_matrix,mild,sev)
        contrast$severe_minus_mild_processed_units <- ex$delta
        contrast$exploratory_exact_p <- ex$p
        contrast$BH_q_within_context <- p.adjust(ex$p,method='BH')
        contrast$permutation_allocations <- ex$n_alloc
        contrast$contrast_status <- 'exploratory_donor_permutation'
        ord <- order(-abs(ex$delta), genes)
        top <- contrast[head(ord,20L),,drop=FALSE]
        top$ranking_note <- 'Largest absolute descriptive effect, NOT independently validated DE genes'
        top_rows[[ctx]] <- top
        logmsg('    Baseline comparison: Mild n=',n_mild,' vs pooled Severe n=',n_sev,
          '; genes=',length(genes),'; min p=',format(min(ex$p),digits=3),
          '; min BH q=',format(min(contrast$BH_q_within_context),digits=3))
      } else {
        logmsg('    SKIPPED severity contrast (Mild n=',n_mild,', Severe n=',n_sev,
               ' first-collection donors passing minimum-cell threshold)')
      }
      all_contrasts[[ctx]] <- contrast
      # Record all sample-level means for a SMALL prespecified marker panel.
      these <- present_markers
      if (length(these)) {
        markerlong <- data.frame(
          lineage=lin, context=ctx,
          gene=rep(these,times=length(samples)),
          sample_id=rep(samples,each=length(these)),
          donor=rep(clinical$donor_inferred,each=length(these)),
          disease_group=rep(clinical$disease_group,each=length(these)),
          visit_code=rep(clinical$visit_code,each=length(these)),
          n_cells=rep(ncs,each=length(these)),
          processed_mean=as.vector(samplemeans[these,,drop=FALSE]),
          stringsAsFactors=FALSE)
        all_markers[[ctx]] <- markerlong
      }
    }
    rm(obj, mat, md)
    invisible(gc(verbose=FALSE))
    logmsg('Released ',lin,' Seurat object before next lineage')
  }

  coverage <- do.call(rbind, all_coverage)
  results <- do.call(rbind, all_contrasts)
  markers <- do.call(rbind, all_markers)
  inventory <- do.call(rbind,missing_markers)
  rownames(coverage) <- NULL
  rownames(results) <- NULL
  rownames(markers) <- NULL
  rownames(inventory) <- NULL
  write.csv(coverage,file.path(out,'8_context_sample_coverage.csv'),row.names=FALSE)
  write.csv(inventory,file.path(out,'8_prespecified_gene_availability.csv'),row.names=FALSE)
  write.csv(markers,file.path(out,'8_marker_expression_by_sample.csv'),row.names=FALSE)
  prespecified <- results[results$gene %in% unique(unlist(marker_panels)),,drop=FALSE]
  write.csv(prespecified,file.path(out,'8_prespecified_gene_contrasts.csv'),row.names=FALSE)
  write.csv(do.call(rbind,top_rows),file.path(out,'8_top20_descriptive_genes_by_context.csv'),row.names=FALSE)
  con <- gzfile(file.path(out,'8_all_gene_exploratory_contrasts.csv.gz'),'wt')
  write.csv(results,con,row.names=FALSE)
  close(con)
  logmsg('Saved all gene contrasts, coverage, marker expression and targeted-gene inventory')

  # Visualization 1: unpooled first-collection donor means, all T and all B.
  group_levels <- c('Healthy','Asymptomatic','Mild','Severe.m','Severe.s')
  selected_genes <- list(T=c('GNLY','GZMB','IFIT1','ISG15'),
                         B=c('JCHAIN','MZB1','FCRL5','IFIT3'))
  f <- markers[markers$context %in% c('T_all','B_all') &
                 markers$sample_id %in% initial$sample_id,,drop=FALSE]
  chosen_marker <- (f$lineage=='T' & f$gene %in% selected_genes$T) |
                   (f$lineage=='B' & f$gene %in% selected_genes$B)
  f <- f[chosen_marker & is.finite(f$processed_mean),,drop=FALSE]
  if (nrow(f)) {
    f$disease_group <- factor(f$disease_group,levels=group_levels)
    f$label <- paste(f$lineage,f$gene,sep=' | ')
    f$label <- factor(f$label,levels=c(paste('T',selected_genes$T,sep=' | '),
                                     paste('B',selected_genes$B,sep=' | ')))
    p <- ggplot2::ggplot(f,ggplot2::aes(x=disease_group,y=processed_mean)) +
      ggplot2::geom_point(size=2,alpha=.85,position=ggplot2::position_jitter(width=.09,height=0)) +
      ggplot2::facet_wrap(~label,scales='free_y',ncol=4) +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=40,hjust=1)) +
      ggplot2::labs(title='First-recorded donor expression: processed RNA means',
                    subtitle='Each point is one donor (small groups; exploratory)',
                    x='Clinical group',y='Mean processed integrated value')
    ggplot2::ggsave(file.path(out,'8_first_collection_gene_expression.png'),p,
                    width=16,height=8,units='in',dpi=160)
    logmsg('Figure: 8_first_collection_gene_expression.png')
  }
  # Visualization 2: prespecified gene differences in whole T / B lineage.
  d <- prespecified[prespecified$context %in% c('T_all','B_all') &
                     prespecified$gene %in% unique(unlist(selected_genes)) &
                     is.finite(prespecified$severe_minus_mild_processed_units),,drop=FALSE]
  if (nrow(d)) {
    d$marker <- paste(d$lineage,d$gene,sep=' | ')
    p <- ggplot2::ggplot(d,ggplot2::aes(x=reorder(marker,severe_minus_mild_processed_units),
          y=severe_minus_mild_processed_units)) +
      ggplot2::geom_hline(yintercept=0,linetype='dashed',color='gray60') +
      ggplot2::geom_segment(ggplot2::aes(xend=marker,y=0,yend=severe_minus_mild_processed_units)) +
      ggplot2::geom_point(size=3) +
      ggplot2::coord_flip() + ggplot2::theme_classic(base_size=12) +
      ggplot2::labs(title='Mild versus severe: donor-mean processed expression',
         subtitle='Severe minus mild. NOT log2 fold change or count-based DE.',
         x='Measured marker',y='Difference in integrated-expression units')
    ggplot2::ggsave(file.path(out,'8_selected_gene_mean_differences.png'),p,
                    width=12,height=7,units='in',dpi=170)
    logmsg('Figure: 8_selected_gene_mean_differences.png')
  }
  # Visualization 3: actual within-patient visit codes (not elapsed days).
  z <- markers[markers$context %in% c('T_all','B_all') &
                 markers$disease_group %in% c('Mild','Severe.m','Severe.s') &
                 ((markers$context=='T_all' & markers$gene %in% c('IFIT1','GNLY')) |
                  (markers$context=='B_all' & markers$gene %in% c('JCHAIN','MZB1'))),,drop=FALSE]
  z <- z[is.finite(z$processed_mean),,drop=FALSE]
  if (nrow(z)) {
    z$gene_panel <- paste(z$lineage,z$gene,sep=' | ')
    z$visit_code <- as.numeric(z$visit_code)
    p <- ggplot2::ggplot(z,ggplot2::aes(x=visit_code,y=processed_mean,
                                    group=donor,color=donor)) +
      ggplot2::geom_line(linewidth=.7,alpha=.8) +
      ggplot2::geom_point(size=1.8) +
      ggplot2::facet_grid(gene_panel~disease_group,scales='free_y') +
      ggplot2::scale_x_continuous(breaks=c(1,2,4)) +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::labs(title='Within-donor processed gene expression over collection codes',
           subtitle='Visit codes 1/2/4 are not established days; descriptive only',
           x='Recorded collection code',y='Mean processed integrated value')
    ggplot2::ggsave(file.path(out,'8_within_donor_marker_trajectories.png'),p,
                    width=15,height=12,units='in',dpi=160)
    logmsg('Figure: 8_within_donor_marker_trajectories.png')
  }
  summary <- aggregate(sampled_context_cells ~ lineage+context+description,
                         coverage, function(x) c(min=min(x),max=max(x),total=sum(x)))
  write.csv(summary,file.path(out,'8_context_cell_count_summary.csv'),row.names=FALSE)
  writeLines(c(
    'GSE196198 | STEP 08 | Exploratory processed RNA gene-expression comparisons',
    'INPUT: Local Step 4 and Step 5 sampled Seurat checkpoints (NOT uploaded ZIPs).',
    'Layer RNA/data contains previously integrated signed author values, not raw UMI counts.',
    'One recorded first collection per inferred donor is used in clinical severity contrasts.',
    'Mild n=2; pooled Severe.m + Severe.s n=5 before per-context minimum-cell filtering.',
    paste('Minimum sampled cells per donor/collection/context:', min_cells),
    'Sample/donor IDs and severity labels cross-checked against Step 6 clinical index.',
    'For each gene, mean integrated value per donor per context; donors are independent units.',
    'Difference: severe minus mild mean of donor means in PROCESSED EXPRESSION UNITS.',
    'Difference is NOT log2FC, raw-count DE, cell-level statistical significance, or absolute RNA abundance.',
    'Exact donor-label permutation over all assignments; two-sided; BH q across measured genes WITHIN CONTEXT.',
    'Some small clusters have too few cells/donors and are intentionally not tested.',
    'No statistically confirmatory inference: tiny n, post-hoc cluster selection, confounding, integration effects.',
    'Multiple contexts tested, but q values are not jointly adjusted across contexts.',
    'Gene rankings highlight hypotheses only; the marker panel is filtered by availability.',
    'Visit codes are collection identifiers; they are not known time intervals.',
    'Next requirement for conventional DE: obtain raw/unintegrated UMI counts, validate sample IDs,',
    'aggregate donor-by-cell-state counts, and fit a donor-aware count model if adequately powered.',
    'Not an exact reproduction of the paper DE results.'
  ),file.path(out,'8_methods_and_limitations.txt'))
  writeLines(capture.output(sessionInfo()), file.path(out,'8_session_info.txt'))
  errorfile <- file.path(out,'8_error.txt')
  if(file.exists(errorfile)) unlink(errorfile)
  archive <- file.path(root,'8_results.zip')
  if(file.exists(archive)) unlink(archive)
  zz <- Sys.which('zip')
  if(!nzchar(zz)) stop('zip executable missing (results are saved in 8_results/)')
  paths <- list.files('8_results',full.names=TRUE,recursive=FALSE)
  paths <- paths[!dir.exists(paths)]
  status <- suppressWarnings(system2(zz,
                   c('-q','-9',shQuote(archive),shQuote(paths))))
  if(!identical(as.integer(status),0L) || !file.exists(archive))
    stop('Zip creation failed; review outputs in 8_results/')
  logmsg('Created upload ZIP: ',archive)
  logmsg('SUCCESS | Step 08 complete. Upload 8_results.zip.')
}

tryCatch(main(),error=function(e){
  msg <- paste('STEP 08 FAILED:',conditionMessage(e))
  cat(msg,'\n',file=stderr())
  args<-commandArgs(trailingOnly=TRUE)
  outdir<-file.path(if(length(args)) args[1] else '.', '8_results')
  dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
  writeLines(msg,file.path(outdir,'8_error.txt'))
  quit(save='no',status=1)
})
