#!/usr/bin/env Rscript
# GSE196198 | STEP 07 | Marker-program consistency and donor-aware comparisons
# INPUT: previously completed 4_results/, 5_results/, and 6_results/ CSV files.
# No raw counts, no integrated-matrix DE, no extra Seurat clustering.
# COMMAND:
#   env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
#   conda run --no-capture-output -n ngs_seurat \
#   Rscript --vanilla 7_validate_severity.R .

options(stringsAsFactors=FALSE, scipen=99)
set.seed(196198L) # reproducible visual jitter; statistical tests are exhaustive
Sys.setenv(OMP_NUM_THREADS='10', OPENBLAS_NUM_THREADS='10', MKL_NUM_THREADS='10')

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  root <- normalizePath(if (length(args)) args[[1]] else '.', mustWork=TRUE)
  setwd(root)
  out <- file.path(root, '7_results')
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  logpath <- file.path(out, '7_run_log.txt')
  if (file.exists(logpath)) unlink(logpath)
  logmsg <- function(...) {
    s <- paste0(..., collapse='')
    cat(s, '\n', sep='')
    cat(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), s, '\n',
        file=logpath, append=TRUE)
  }
  save_table <- function(x, name) {
    write.csv(x, file=file.path(out, name), row.names=FALSE, na='')
    logmsg('Table: ', name, ' (', nrow(x), ' rows)')
  }
  save_plot <- function(plot, name, width=12, height=7) {
    ggplot2::ggsave(filename=file.path(out,name), plot=plot,
                    width=width,height=height,units='in',dpi=170,limitsize=FALSE)
    logmsg('Figure: ',name)
  }
  logmsg('STEP 07 | GSE196198 | Marker evidence and donor-level severity comparisons')
  if (!requireNamespace('ggplot2', quietly=TRUE)) stop('ggplot2 is required')
  logmsg('R: ',R.version.string,' | ggplot2: ',as.character(packageVersion('ggplot2')))

  input <- function(step, filename) {
    p <- file.path(root,sprintf('%d_results',step), filename)
    if (!file.exists(p)) stop('Missing required input: ',p)
    d <- utils::read.csv(p, check.names=FALSE, stringsAsFactors=FALSE)
    if (!nrow(d)) stop('Input table is empty: ',p)
    d
  }
  require_columns <- function(d, cols, name) {
    absent <- setdiff(cols,names(d))
    if (length(absent)) stop(name,' is missing columns: ',paste(absent,collapse=', '))
  }
  
  # -------------------------------------------------------------
  # A. Check original cluster marker tables and metadata integrity
  # -------------------------------------------------------------
  read_lineage <- function(step, lineage) {
    means <- input(step, sprintf('%d_marker_mean_by_cluster.csv',step))
    available <- input(step,sprintf('%d_marker_availability.csv',step))
    sizes <- input(step,sprintf('%d_cluster_sizes.csv',step))
    require_columns(means,c('gene','cluster','processed_mean'),'marker means')
    require_columns(available,c('gene','present'),'marker availability')
    require_columns(sizes,c('cluster'),'cluster sizes')
    means$gene <- as.character(means$gene)
    means$cluster <- as.character(means$cluster)
    sizes$cluster <- as.character(sizes$cluster)
    available$gene <- as.character(available$gene)
    if (anyNA(means[c('gene','cluster','processed_mean')]) ||
        any(!is.finite(means$processed_mean)) ||
        anyDuplicated(means[c('gene','cluster')]) ||
        anyDuplicated(available$gene) || anyNA(available$present) ||
        anyDuplicated(sizes$cluster)) {
      stop('Invalid or duplicate ',lineage,' marker/cluster records')
    }
    if (!setequal(unique(means$gene),available$gene[available$present]))
      stop(lineage,' marker means do not correspond to measured marker genes')
    if (!setequal(unique(means$cluster),sizes$cluster))
      stop(lineage,' cluster identifiers differ between marker and size tables')
    ordered_clusters <- as.character(sort(as.integer(unique(means$cluster))))
    if (anyNA(as.integer(means$cluster))) stop('Non-numeric cluster labels in Step ',step)
    genes <- sort(unique(means$gene))
    expected_rows <- length(ordered_clusters)*length(genes)
    if (nrow(means)!=expected_rows)
      stop('Incomplete ',lineage,' marker-by-cluster matrix: ',nrow(means),' vs ',expected_rows)
    mat <- matrix(NA_real_,nrow=length(genes),ncol=length(ordered_clusters),
                  dimnames=list(genes,ordered_clusters))
    r <- match(means$gene, genes)
    c <- match(means$cluster,ordered_clusters)
    if (anyNA(r) || anyNA(c)) stop('Cannot align marker rows')
    mat[cbind(r,c)] <- means$processed_mean
    if (any(!is.finite(mat))) stop('Missing marker means after alignment')
    z <- t(scale(t(mat)))
    z[!is.finite(z)] <- 0
    if (!identical(dim(z),dim(mat))) stop('Marker Z-score matrix dimensions lost')
    dimnames(z) <- dimnames(mat)
    logmsg(lineage,': ',ncol(mat),' clusters; ',nrow(mat),' measured marker genes')
    list(raw=mat,z=z,available=available,clusters=ordered_clusters)
  }
  T <- read_lineage(4,'T')
  B <- read_lineage(5,'B')
  
  # Predefined for this script from Step 6 exploratory post-hoc candidate clusters;
  # signatures describe those same processed RNA markers, NOT independent validation.
  # Extra naive-like references are contextual and are NOT included in severity tests.
  specs <- list(
    T_Naive_context=list(lineage='T',clusters=c('0','7'),
         markers=c('CCR7','TCF7','LEF1','IL7R','LTB'),severity=FALSE),
    T_Cytotoxic_2_10=list(lineage='T',clusters=c('2','10'),
         markers=c('NKG7','GNLY','PRF1','GZMB','GZMK','CX3CR1'),severity=TRUE),
    T_Cycling_13=list(lineage='T',clusters='13',
         markers=c('MKI67','CD38'),severity=TRUE),
    T_Interferon_16=list(lineage='T',clusters='16',
         markers=c('ISG15','IFIT1','IFIT3','IFI6','MX1','STAT1'),severity=TRUE),
    B_Naive_context=list(lineage='B',clusters=c('2','10','11'),
         markers=c('TCL1A','IGHD','IGHM','FCER2','IL4R'),severity=FALSE),
    B_Secretory_13=list(lineage='B',clusters='13',
         markers=c('JCHAIN','MZB1','XBP1','PRDM1','IGHG1'),severity=TRUE),
    B_Atypical_6=list(lineage='B',clusters='6',
         markers=c('FCRL5','ITGAX','TBX21','CD86'),severity=TRUE),
    B_Interferon_14=list(lineage='B',clusters='14',
         markers=c('ISG15','IFIT3','MX1'),severity=TRUE)
  )
  candidate_names <- names(specs)[vapply(specs,function(x) isTRUE(x$severity),logical(1))]
  
  # Marker coverage and transparent gene-level support for each cluster.
  availability_rows <- list()
  marker_rows <- list()
  score_rows <- list()
  evidence_rows <- list()
  candidate_gene_rows <- list()
  k <- 0L
  for (name in names(specs)) {
    spec <- specs[[name]]
    lineage <- if (spec$lineage=='T') T else B
    missing <- setdiff(spec$markers,rownames(lineage$z))
    found <- intersect(spec$markers,rownames(lineage$z))
    if (length(found)<2L) stop(name,': fewer than 2 panel markers found')
    if (!all(spec$clusters %in% lineage$clusters))
      stop(name,': selected candidate cluster is missing')
    availability_rows[[name]] <- data.frame(
      program=name,lineage=spec$lineage,
      requested_markers=length(spec$markers),measured_markers=length(found),
      marker_coverage_pct=100*length(found)/length(spec$markers),
      present=paste(found,collapse=';'),missing=paste(missing,collapse=';'),
      expected_clusters=paste(spec$clusters,collapse='+'),
      included_in_severity=spec$severity)
    for (g in spec$markers) {
      marker_rows[[length(marker_rows)+1L]] <- data.frame(
        program=name,lineage=spec$lineage,gene=g,available=g %in% found,
        expected_cluster=paste(spec$clusters,collapse='+'))
    }
    values <- lineage$z[found,,drop=FALSE]
    scores <- colMeans(values)
    ranked <- setNames(rank(-scores,ties.method='min'),names(scores))
    for (cl in lineage$clusters) {
      v <- values[,cl]
      score_rows[[length(score_rows)+1L]] <- data.frame(
        program=name,lineage=spec$lineage,cluster=cl,
        is_original_candidate=cl %in% spec$clusters,
        score_mean_gene_z=unname(scores[cl]),
        score_rank_within_lineage=unname(ranked[cl]),
        markers_measured=length(found),markers_z_at_least_0_5=sum(v>=0.5),
        fraction_markers_z_at_least_0_5=mean(v>=0.5))
      if (cl %in% spec$clusters) {
        # Export the evidence for every gene, including weak/negative markers.
        for (g in found) {
          candidate_gene_rows[[length(candidate_gene_rows)+1L]] <- data.frame(
            program=name,lineage=spec$lineage,candidate_cluster=cl,gene=g,
            processed_mean=unname(lineage$raw[g,cl]),
            gene_z_relative_to_clusters=unname(lineage$z[g,cl]),
            supports_threshold_z_0_5=unname(lineage$z[g,cl])>=0.5)
        }
        # Leave-one-marker-out: does a candidate remain a top-three cluster?
        # With two-marker panels this sensitivity diagnostic is especially fragile.
        loo_rank <- vapply(found,function(drop_gene) {
          keep <- setdiff(found,drop_gene)
          other_scores <- colMeans(lineage$z[keep,,drop=FALSE])
          ranks <- setNames(rank(-other_scores,ties.method='min'),names(other_scores))
          unname(ranks[cl])
        },numeric(1))
        evidence_rows[[length(evidence_rows)+1L]] <- data.frame(
          program=name,lineage=spec$lineage,candidate_cluster=cl,
          mean_marker_z=unname(scores[cl]),
          marker_support=paste0(sum(v>=0.5),'/',length(v)),
          rank_of_17_clusters=unname(ranked[cl]),
          leave_one_gene_out_top3_fraction=mean(loo_rank<=3),
          marker_consistency=if (scores[cl]>=0.75 && sum(v>=0.5)>=ceiling(length(v)/2))
            'consistent with marker program' else 'mixed or weak marker evidence',
          interpretive_status='RNA-only, same-data marker corroboration; NOT independent subtype validation')
      }
    }
  }
  coverage <- do.call(rbind,availability_rows)
  scores <- do.call(rbind,score_rows)
  evidence <- do.call(rbind,evidence_rows)
  save_table(coverage,'7_marker_coverage.csv')
  save_table(do.call(rbind,marker_rows),'7_requested_marker_inventory.csv')
  save_table(scores,'7_all_clusters_marker_program_scores.csv')
  save_table(evidence,'7_candidate_cluster_marker_evidence.csv')
  save_table(do.call(rbind,candidate_gene_rows),'7_candidate_cluster_gene_level_evidence.csv')
  
  # Matrix-style heatmap for every cluster, with a cross marking the clusters
  # selected in the previous step. Gene standardization is across clusters.
  plot_signature <- function(source) {
    d <- scores[scores$lineage==source,]
    orders <- names(specs)[vapply(specs,function(s) s$lineage==source,logical(1))]
    d$program <- factor(d$program,levels=rev(orders))
    d$cluster <- factor(d$cluster,levels=as.character(sort(as.integer(unique(d$cluster)))))
    selected <- d[d$is_original_candidate,]
    ggplot2::ggplot(d,ggplot2::aes(x=cluster,y=program,fill=score_mean_gene_z)) +
      ggplot2::geom_tile(color='white',linewidth=.25) +
      ggplot2::geom_point(data=selected,ggplot2::aes(x=cluster,y=program),
                          inherit.aes=FALSE,shape=4,size=2.7,stroke=.8,color='black') +
      ggplot2::scale_fill_gradient2(low='#39629a',mid='white',high='#b63e40',
                                    midpoint=0,name='Mean gene Z') +
      ggplot2::labs(title=paste0(source,'-cell marker signature consistency'),
                    subtitle='Crosses: clusters selected before this step; scores standardized across clusters',
                    x='New Seurat cluster',y='Marker program') +
      ggplot2::theme_minimal(base_size=12) +
      ggplot2::theme(panel.grid=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_text(angle=0))
  }
  save_plot(plot_signature('T'),'7_T_marker_signature_heatmap.png',13,5.8)
  save_plot(plot_signature('B'),'7_B_marker_signature_heatmap.png',13,5.8)

  # -------------------------------------------------------------
  # B. Verify Step 6 proportions against Step 4/5 cluster count tables
  # -------------------------------------------------------------
  trajectories <- input(6,'6_sample_level_program_fractions.csv')
  require_columns(trajectories,
      c('sample_id','program','lineage','n_program','n_lineage_sampled','percentage',
        'donor_inferred','severity2','disease_group','visit_code'), 'Step 6 trajectories')
  if (anyNA(trajectories) || anyDuplicated(trajectories[c('sample_id','program')]))
    stop('Missing values or duplicate sample-program entries in Step 6')
  if (!setequal(unique(trajectories$program),candidate_names))
    stop('Step 6 program list differs from this Step 7 predefined candidate list')
  if (any(trajectories$n_program>trajectories$n_lineage_sampled) ||
      any(trajectories$n_program<0) ||
      any(trajectories$n_lineage_sampled<=0) ||
      any(abs(trajectories$percentage -
                100*trajectories$n_program/trajectories$n_lineage_sampled)>1e-7))
    stop('Step 6 cell proportions are inconsistent with their counts')
  for (step in c(4,5)) {
    src <- if (step==4) 'T' else 'B'
    d <- input(step,sprintf('%d_clusters_by_sample_descriptive.csv',step))
    require_columns(d,c('sample_id','cluster','n_cells','sample_total'),
                    paste0('Step ',step,' sample counts'))
    d$cluster <- as.character(d$cluster)
    if (anyDuplicated(d[c('sample_id','cluster')])) stop('Duplicate Step ',step,' sample-cluster counts')
    for (program in candidate_names) {
      s <- specs[[program]]
      if (s$lineage!=src) next
      x <- trajectories[trajectories$program==program,]
      for (i in seq_len(nrow(x))) {
        samp <- x$sample_id[i]
        sub <- d[d$sample_id==samp,]
        if (!nrow(sub)) stop('Step ',step,' lacks sample ',samp)
        n <- sum(sub$n_cells[sub$cluster %in% s$clusters])
        total <- unique(sub$sample_total)
        if (length(total)!=1L || n!=x$n_program[i] || total!=x$n_lineage_sampled[i])
          stop('Step 6 and Step ',step,' disagree for ',program,' / ',samp)
      }
    }
  }
  logmsg('Crosschecked all Step 6 candidate population counts against Steps 4 and 5.')

  # -------------------------------------------------------------
  # C. One initial collection per donor for cross-sectional comparisons
  # -------------------------------------------------------------
  clinical <- input(6,'6_sample_donor_clinical_index.csv')
  require_columns(clinical,c('sample_id','donor_inferred','severity2',
                             'disease_group','visit_code'),'clinical index')
  if (anyDuplicated(clinical$sample_id) ||
      anyDuplicated(clinical[c('donor_inferred','visit_code')]))
    stop('Duplicated donor/sample metadata')
  if (!setequal(unique(clinical$sample_id),unique(trajectories$sample_id)))
    stop('Step 6 clinical and abundance sample IDs differ')
  observed <- unique(trajectories[c('sample_id','donor_inferred','severity2',
                                  'disease_group','visit_code')])
  chk <- merge(observed,clinical,by='sample_id',suffixes=c('.t','.c'),sort=FALSE)
  if (nrow(chk)!=nrow(clinical)) stop('Unexpected clinical index row count')
  for (nm in c('donor_inferred','severity2','disease_group','visit_code'))
    if (!identical(as.character(chk[[paste0(nm,'.t')]]),
                   as.character(chk[[paste0(nm,'.c')]])))
      stop('Clinical labels disagree: ',nm)
  donor_groups <- split(clinical$disease_group,clinical$donor_inferred)
  if (any(vapply(donor_groups,function(x) length(unique(x))!=1L,logical(1))))
    stop('Donor changes disease-group label across collections')
  initial <- do.call(rbind,lapply(split(clinical,clinical$donor_inferred),function(x) {
    x <- x[order(x$visit_code),,drop=FALSE]
    x[1L,,drop=FALSE]
  }))
  rownames(initial) <- NULL
  if (anyDuplicated(initial$donor_inferred)) stop('Duplicate first collection donor')
  group_counts <- as.data.frame(table(group=initial$disease_group),stringsAsFactors=FALSE)
  group_counts <- group_counts[group_counts$Freq>0,]
  names(group_counts)[2] <- 'n_donors'
  save_table(group_counts,'7_baseline_donors_by_clinical_group.csv')
  first <- trajectories[trajectories$sample_id %in% initial$sample_id,,drop=FALSE]
  if (nrow(first)!=nrow(initial)*length(candidate_names))
    stop('Missing initial donor-program rows')
  save_table(first,'7_first_collection_donor_program_fractions.csv')
  
  group_levels <- c('Healthy','Asymptomatic','Mild','Severe.m','Severe.s')
  first$disease_group <- factor(first$disease_group,levels=group_levels)
  first$program <- factor(first$program,levels=c('T_Cytotoxic_2_10','T_Cycling_13',
                 'T_Interferon_16','B_Secretory_13','B_Atypical_6','B_Interferon_14'))
  p_baseline <- ggplot2::ggplot(first,ggplot2::aes(x=disease_group,y=percentage)) +
    ggplot2::geom_point(ggplot2::aes(color=donor_inferred),size=2.6,
                       position=ggplot2::position_jitter(width=.06,height=0)) +
    ggplot2::facet_wrap(~program,scales='free_y',ncol=3) +
    ggplot2::labs(title='First recorded collection: donor-level immune cluster fractions',
                  subtitle='One point per donor; original clinical groups retained',
                  x='Author clinical code',
                  y='Percentage of sampled B cells or T cells within lineage',color='Donor') +
    ggplot2::theme_classic(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1))
  save_plot(p_baseline,'7_baseline_program_fractions_by_severity.png',17,9)

  group_rows <- lapply(split(first,interaction(first$program,first$disease_group,drop=TRUE)),
    function(d) {
      v <- d$percentage
      data.frame(program=as.character(d$program[1]),lineage=d$lineage[1],
        disease_group=as.character(d$disease_group[1]),n_donors=length(v),
        mean_percent=mean(v),median_percent=median(v),min_percent=min(v),
        max_percent=max(v),sd_percent=if(length(v)>1L) stats::sd(v) else NA_real_)
    })
  group_summary <- do.call(rbind,group_rows)
  rownames(group_summary) <- NULL
  save_table(group_summary,'7_first_collection_group_descriptive.csv')

  # -------------------------------------------------------------
  # D. EXACT donor-label permutation comparisons: only first collections
  # -------------------------------------------------------------
  # With 2 mild vs 5 severe donors: only choose(7,2)=21 assignments.
  # These p-values cannot establish clinical associations: post-hoc cluster
  # selection, extremely low power, unmeasured confounding, multiple tests.
  exact_perm <- function(x,y) {
    x <- as.numeric(x); y <- as.numeric(y)
    n1 <- length(x); n2 <- length(y)
    if (n1<2L || n2<2L) return(c(delta=NA_real_,p=NA_real_,nperm=0))
    values <- c(x,y)
    if (choose(n1+n2,n1)>50000) stop('Exact permutation set unexpectedly large')
    combos <- utils::combn(n1+n2,n1)
    effects <- apply(combos,2L,function(idx) {
      mean(values[-idx])-mean(values[idx])
    })
    obs <- mean(y)-mean(x)
    # Two-sided conservative discrete exact p-value, not a t-test on cells.
    p <- mean(abs(effects)>=abs(obs)-1e-10)
    c(delta=obs,p=p,nperm=ncol(combos))
  }
  tests <- list(
    list(label='Mild vs Severe (pooled)',A='Mild',B=c('Severe.m','Severe.s')),
    list(label='Severe.m vs Severe.s',A='Severe.m',B='Severe.s')
  )
  test_rows <- list()
  for (program in candidate_names) {
    rows <- first[as.character(first$program)==program,]
    for (comparison in tests) {
      x <- rows$percentage[rows$disease_group %in% comparison$A]
      y <- rows$percentage[rows$disease_group %in% comparison$B]
      test <- exact_perm(x,y)
      test_rows[[length(test_rows)+1L]] <- data.frame(
        program=program,lineage=as.character(rows$lineage[1]),
        contrast=comparison$label,
        group_A=paste(comparison$A,collapse='+'),group_B=paste(comparison$B,collapse='+'),
        n_donors_A=length(x),n_donors_B=length(y),
        mean_percent_A=mean(x),mean_percent_B=mean(y),
        delta_percentage_points_B_minus_A=unname(test['delta']),
        exact_two_sided_permutation_p=unname(test['p']),
        number_possible_assignments=as.integer(test['nperm']))
    }
  }
  test_table <- do.call(rbind,test_rows)
  test_table$BH_q_across_12_exploratory_tests <- stats::p.adjust(
    test_table$exact_two_sided_permutation_p,method='BH')
  save_table(test_table,'7_exploratory_exact_donor_permutation_tests.csv')
  logmsg('Exploratory exact permutations computed (discrete, tiny n; NOT confirmatory).')
  
  effect_plot <- test_table
  effect_plot$program <- factor(effect_plot$program,levels=rev(candidate_names))
  pe <- ggplot2::ggplot(effect_plot,
    ggplot2::aes(x=delta_percentage_points_B_minus_A,y=program)) +
    ggplot2::geom_vline(xintercept=0,linetype=2,color='grey50') +
    ggplot2::geom_point(size=3,color='#264f81') +
    ggplot2::facet_wrap(~contrast,scales='free_x',ncol=2) +
    ggplot2::labs(title='Exploratory first-collection group differences',
      subtitle='B minus A group mean; unweighted donor percentages, very small sample sizes',
      x='Mean difference (percentage points)',y='Candidate population') +
    ggplot2::theme_classic(base_size=12)
  save_plot(pe,'7_first_collection_mean_differences.png',15,6)

  # -------------------------------------------------------------
  # E. Consistency of first-to-last within-donor changes (from Step 6)
  # -------------------------------------------------------------
  changes <- input(6,'6_within_donor_first_to_last_change.csv')
  require_columns(changes,c('donor_inferred','disease_group','program','lineage',
                           'first_visit','last_visit','first_percent','last_percent',
                           'delta_percentage_points'),'Step 6 paired donor table')
  if (anyDuplicated(changes[c('donor_inferred','program')]))
    stop('Duplicate donor/program pairs in Step 6 changes')
  if (any(abs(changes$delta_percentage_points-
                  (changes$last_percent-changes$first_percent))>1e-7))
    stop('Step 6 paired percentage changes do not reconcile')
  change_summary <- do.call(rbind,lapply(
    split(changes,interaction(changes$program,changes$disease_group,drop=TRUE)),
    function(d) data.frame(program=d$program[1],lineage=d$lineage[1],
      disease_group=d$disease_group[1],n_donors=nrow(d),
      mean_delta_pp=mean(d$delta_percentage_points),
      median_delta_pp=median(d$delta_percentage_points),
      n_increased=sum(d$delta_percentage_points>0),
      n_decreased=sum(d$delta_percentage_points<0),
      n_unchanged=sum(d$delta_percentage_points==0))))
  rownames(change_summary) <- NULL
  save_table(change_summary,'7_donor_change_direction_summary.csv')
  changes$disease_group <- factor(changes$disease_group,
                                  levels=c('Mild','Severe.m','Severe.s'))
  changes$program <- factor(changes$program,levels=rev(candidate_names))
  pchange <- ggplot2::ggplot(changes,
         ggplot2::aes(x=delta_percentage_points,y=program)) +
     ggplot2::geom_vline(xintercept=0,linetype=2,color='grey55') +
     ggplot2::geom_point(ggplot2::aes(color=donor_inferred),size=2.5,
       position=ggplot2::position_jitter(height=.10,width=0)) +
     ggplot2::facet_wrap(~disease_group,ncol=3) +
     ggplot2::labs(title='Individual donors: change from first to last collection',
                   subtitle='Collection codes 1 to 4 are NOT verified elapsed days',
                   x='Last minus first fraction (percentage points)',y='Candidate program',
                   color='Donor') + ggplot2::theme_classic(base_size=11)
  save_plot(pchange,'7_within_donor_changes_by_group.png',17,7)
  
  # -------------------------------------------------------------
  # F. Scientific boundaries and sharing
  # -------------------------------------------------------------
  cautions <- c(
    'GSE196198 | STEP 07 | Marker evidence and exploratory donor-aware severity comparisons',
    'Inputs: completed Step 4/5 marker mean and sample-composition CSVs; Step 6 program/clinical CSVs.',
    paste('Donors:',nrow(initial),'| Initial donor counts:',
      paste(paste(group_counts$group,group_counts$n_donors,sep='='),collapse='; ')),
    'New T/B clusters are from proportionate 40k-cell samples per lineage, not the complete cohort.',
    'Marker program Z-scores are computed FROM THE SAME processed RNA matrix that was used to select clusters.',
    'Marker consistency uses an arbitrary descriptive rule: mean gene z >=0.75 and at least half genes z >=0.5.',
    'Leave-one-marker-out top-three fraction is a signature sensitivity check, NOT validation on new data.',
    'Therefore marker agreement is corroboration of an exploratory label, NOT independent subtype validation.',
    'Scores represent across-cluster standardized gene means and cannot be treated as single-cell signature scores.',
    'One first-recorded collection per donor is used for all cross-sectional severity contrasts.',
    'Mild versus pooled Severe is distinct from Severe.m versus Severe.s; source labels are preserved.',
    'The clinical meaning of Severe.m versus Severe.s has NOT been independently established.',
    'Permutation comparisons use DONORS as the exchangeable unit and unweighted sample percentages.',
    'Even if p<0.05, any test here is exploratory due to post-hoc clusters and donor sample sizes of 2 to 5.',
    'Exact p values are discrete and very low-powered; BH q values are across 12 exploratory tests.',
    'Disease severity, patient, treatment, batch and collection code may be confounded.',
    'Within-donor collection identifiers 1/2/4 are not established elapsed days.',
    'T and B percentages use DIFFERENT denominators (sampled cells of their own lineage).',
    'No raw UMI counts, valid pseudo-bulk DE, mass-normalized gene DE, CITE-seq or receptor validation is performed.',
    'In particular, this stage is NOT a full reproduction of the original paper subtype or severity figures.',
    'The processed RNA marker panel has missing genes and author integration may reduce biological contrasts.'
  )
  writeLines(cautions,file.path(out,'7_methods_limitations.txt'))
  writeLines(capture.output(sessionInfo()),file.path(out,'7_session_info.txt'))
  olderr <- file.path(out,'7_error.txt'); if (file.exists(olderr)) unlink(olderr)
  zip <- Sys.which('zip')
  if (nzchar(zip)) {
    dest <- file.path(root,'7_results.zip')
    if (file.exists(dest)) unlink(dest)
    files <- list.files('7_results',full.names=TRUE,recursive=FALSE)
    files <- files[!dir.exists(files)]
    status <- suppressWarnings(system2(zip,c('-q','-9',shQuote(dest),shQuote(files))))
    if (identical(as.integer(status),0L) && file.exists(dest))
      logmsg('Created archive: ',dest)
    else logmsg('WARNING: Could not build ZIP; files remain in 7_results/')
  } else logmsg('WARNING: zip unavailable, outputs remain in 7_results/')
  logmsg('SUCCESS: Step 07 complete. Send 7_results.zip for review.')
}

tryCatch(main(),error=function(e) {
  s <- paste0('STEP 07 FAILED: ',conditionMessage(e))
  cat(s,'\n',file=stderr())
  a <- commandArgs(trailingOnly=TRUE)
  dd <- if(length(a)) a[[1]] else '.'
  od <- file.path(dd,'7_results')
  dir.create(od,recursive=TRUE,showWarnings=FALSE)
  writeLines(c(s,capture.output(sys.calls())),file.path(od,'7_error.txt'))
  quit(save='no',status=1)
})
