#!/usr/bin/env Rscript
# STEP 09 | Kim et al. 2023, GSE196198 | TCR/BCR repertoire and RNA-cluster linkage
# Paper: Supplementary S1(A-D), selected receptor concepts from S5 and S7.
# Input: original receptor workbooks, PBMC/T/B metadata, local Step4/5 RDS files.
# R >=4.0; needs readxl and ggplot2. Processed RNA is not normalized here.
# Run: Rscript --vanilla 9_receptor_repertoire.R .
# To replot from saved local checkpoint: GSE_RECEPTOR_REUSE=1 Rscript ...

options(stringsAsFactors=FALSE, scipen=99)
Sys.setenv(OMP_NUM_THREADS='10', OPENBLAS_NUM_THREADS='10',
           MKL_NUM_THREADS='10', VECLIB_MAXIMUM_THREADS='10')

run <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  root <- normalizePath(if (length(args)) args[[1]] else '.', mustWork=TRUE)
  setwd(root)
  od <- file.path(root,'9_results')
  dir.create(od, recursive=TRUE, showWarnings=FALSE)
  logfile <- file.path(od,'9_run_log.txt')
  cat('', file=logfile)
  say <- function(...) {
    line <- paste0(...,collapse='')
    cat(line,'\n',sep='')
    cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),line,'\n',file=logfile,append=TRUE)
  }
  say('STEP 09 | Kim et al. (2023) | paired TCR/BCR, longitudinal clone tracking')
  if (!requireNamespace('ggplot2',quietly=TRUE)) stop('ggplot2 is required')
  if (!requireNamespace('readxl',quietly=TRUE)) {
    stop("Package 'readxl' is missing. Install with: conda install -n ngs_seurat -c conda-forge r-readxl")
  }
  checkpoint <- file.path(od,'9_receptor_checkpoint.rds')
  reuse <- identical(Sys.getenv('GSE_RECEPTOR_REUSE','0'),'1')
  save_csv <- function(x,name) utils::write.csv(x,file.path(od,name),row.names=FALSE,na='')
  save_gz <- function(x,name) {
    con <- gzfile(file.path(od,name),'wt')
    on.exit(close(con))
    utils::write.csv(x,con,row.names=FALSE,na='')
  }
  metadata <- function(type) {
    path <- file.path(root,paste0('GSE196198_',type,'_meta.data.txt.gz'))
    if (!file.exists(path)) stop('Missing metadata: ',path)
    con <- gzfile(path,'rt')
    on.exit(close(con))
    x <- utils::read.table(con,header=TRUE,row.names=1,sep='',quote='"',
                          comment.char='',check.names=FALSE,stringsAsFactors=FALSE)
    needed <- c('id.match','orig.ident','severity','severity2')
    if (!all(needed %in% names(x))) stop('Metadata missing: ',paste(setdiff(needed,names(x)),collapse=', '))
    if (anyNA(x$id.match) || anyDuplicated(x$id.match)) stop('Missing or duplicate canonical IDs: ',type)
    x
  }
  parse_sample <- function(sample_id) {
    label <- sub('^Exp_','',sample_id)
    donor <- sub('_.*$','',label)
    visit <- ifelse(grepl('_',label),sub('^[^_]+_','',label),'1')
    if (any(!grepl('^[HCMS][0-9]+$',donor)) || any(!grepl('^[0-9]+$',visit))) {
      stop('Unexpected sample identifiers. Inspect PBMC orig.ident before assigning visits.')
    }
    data.frame(donor=donor,visit_index=as.integer(visit))
  }
  truth <- function(x,field,receptor) {
    val <- tolower(trimws(as.character(x)))
    bad <- !val %in% c('true','false','1','0','', 'na','nan','none')
    if (any(bad)) stop(receptor,': unrecognized ',field,' flag: ',head(unique(val[bad]),1))
    val %in% c('true','1')
  }
  required_receptor <- c('Sample ID','barcode','contig_id','chain','v_gene','j_gene','cdr3','cdr3_nt')
  
  receptor <- function(label, suffix, pbmc, lineage_ids) {
    filename <- file.path(root,paste0('GSE196198_',label,'_seq_data.xlsx'))
    if (!file.exists(filename)) stop('Missing receptor workbook: ',filename)
    sheets <- readxl::excel_sheets(filename)
    if (!label %in% sheets) stop('Sheet ',label,' not found in ',basename(filename),'. Sheets: ',paste(sheets,collapse=', '))
    say('Reading ',label,' contigs from sheet ',label,' (may take several minutes)...')
    d <- as.data.frame(readxl::read_excel(filename,sheet=label,col_types='text',.name_repair='minimal'))
    if (!all(required_receptor %in% names(d))) {
      stop(label,' workbook missing columns: ',paste(setdiff(required_receptor,names(d)),collapse=', '))
    }
    for (field in required_receptor) {
      d[[field]] <- trimws(as.character(d[[field]]))
      d[[field]][is.na(d[[field]])] <- ''
    }
    n_orig <- nrow(d)
    if (!n_orig) stop(label,' workbook is empty')
    if (any(d$contig_id=='') || anyDuplicated(paste(d[['Sample ID']],d$contig_id,sep='|'))) {
      stop(label,' has empty/duplicate contig_id within a sample')
    }
    pattern <- if (suffix=='T') '^([HC][0-9]+_T|[MS][0-9]+_[0-9]+T)$' else
                              '^([HC][0-9]+_B|[MS][0-9]+_[0-9]+B)$'
    if (any(!grepl(pattern,d[['Sample ID']])) || any(!grepl('^[ACGT]+-[0-9]+$',d$barcode))) {
      stop(label,' unexpected Sample ID or barcode. Example: ',paste(head(unique(d[['Sample ID']]),5),collapse=', '))
    }
    sample_lab <- sub(paste0('_?',suffix,'$'),'',d[['Sample ID']])
    d$sample_id <- paste0('Exp_',sample_lab)
    d$id.match <- paste0(sample_lab,'_',d$barcode)
    mi <- match(d$id.match,pbmc$id.match)
    d$rna_match <- !is.na(mi)
    if (any(d$rna_match & d$sample_id!=pbmc$orig.ident[mi],na.rm=TRUE)) {
      stop(label,' receptor sample and PBMC identity disagree')
    }
    if (any(!d$sample_id %in% pbmc$orig.ident)) stop(label,' contains unknown sample labels')
    flags <- c('is_cell','high_confidence','full_length','productive')
    present_flags <- intersect(flags,names(d))
    if (label=='BCR' && length(present_flags)!=length(flags)) {
      stop('BCR workbook lacks required contig QC flags: ',paste(setdiff(flags,present_flags),collapse=','))
    }
    good_flag <- rep(TRUE,nrow(d))
    flags_report <- data.frame(receptor=character(),flag=character(),present=logical(),
                               passing=integer(),total=integer())
    for (field in flags) {
      if (field %in% present_flags) {
        ok <- truth(d[[field]],field,label)
        good_flag <- good_flag & ok
        flags_report <- rbind(flags_report,data.frame(receptor=label,flag=field,
                                      present=TRUE,passing=sum(ok),total=nrow(d)))
      } else {
        flags_report <- rbind(flags_report,data.frame(receptor=label,flag=field,
                                      present=FALSE,passing=NA_integer_,total=nrow(d)))
      }
    }
    allowed <- if(label=='TCR') c('TRA','TRB','TRG','TRD') else c('IGH','IGK','IGL')
    allowed_chain <- d$chain %in% allowed
    valid_vj <- nzchar(d$v_gene) & nzchar(d$j_gene) &
                grepl('^[A-Z0-9*/.-]+$',d$v_gene) & grepl('^[A-Z0-9*/.-]+$',d$j_gene)
    for (ch in allowed) {
      idx <- which(d$chain==ch)
      if (length(idx)) valid_vj[idx] <- valid_vj[idx] &
                            startsWith(d$v_gene[idx],paste0(ch,'V')) &
                            startsWith(d$j_gene[idx],paste0(ch,'J'))
    }
    valid_seq <- grepl('^[ACGT]+$',d$cdr3_nt) &
                 grepl('^[ACDEFGHIKLMNPQRSTVWY]+$',d$cdr3) &
                 nchar(d$cdr3_nt)==3L*nchar(d$cdr3)
    d$eligible <- good_flag & allowed_chain & valid_vj & valid_seq
    # With a flag-passing but invalid contig, another receptor chain may be hidden.
    unresolved <- unique(d$id.match[good_flag & !d$eligible])
    all_cells <- unique(d$id.match)
    # Exclude duplicated identical contigs from chain multiplicity, but never collapse
    # biologically distinct TRA/TRB/IGH/IGK/IGL sequences into one.
    u <- d[d$eligible,c('id.match','chain','v_gene','j_gene','cdr3_nt','cdr3')]
    u <- unique(u)
    if (!nrow(u)) stop(label,' has no eligible receptor chains')
    tab <- table(factor(u$id.match,levels=all_cells),factor(u$chain,levels=allowed))
    counts <- matrix(as.integer(tab),nrow=nrow(tab),dimnames=dimnames(tab))
    n_total <- rowSums(counts)
    status <- ifelse(n_total==0,'no_eligible_chain',ifelse(n_total==1,'single_chain','multiple_chains'))
    subtype <- rep('unresolved',length(all_cells))
    choose <- list()
    definitions <- if(label=='TCR') list(c('TCR_ab','TRA','TRB'),c('TCR_gd','TRG','TRD')) else
                                    list(c('BCR_IGK','IGH','IGK'),c('BCR_IGL','IGH','IGL'))
    for (def in definitions) {
      yes <- n_total==2 & counts[,def[2]]==1 & counts[,def[3]]==1
      status[yes] <- 'single_pair'; subtype[yes] <- def[1]
      choose[[def[1]]] <- all_cells[yes]
    }
    status[all_cells %in% unresolved] <- 'unresolved_contig'
    subtype[status!='single_pair'] <- 'unresolved'
    cell <- data.frame(id.match=all_cells, pair_status=status, receptor_subtype=subtype,
                       eligible_chains=n_total,stringsAsFactors=FALSE)
    cell$rna_match <- cell$id.match %in% pbmc$id.match
    cell$in_author_subset <- cell$id.match %in% lineage_ids
    cell$clone_key <- NA_character_
    cell$v1 <- cell$j1 <- cell$v2 <- cell$j2 <- NA_character_
    cell$cdr3_1 <- cell$cdr3_2 <- NA_character_
    # Match one unique first and second chain (not just matching a CDR3 sequence).
    for (def in definitions) {
      these <- cell$id.match[cell$pair_status=='single_pair' & cell$receptor_subtype==def[1]]
      if (!length(these)) next
      a <- u[u$id.match %in% these & u$chain==def[2],]
      b <- u[u$id.match %in% these & u$chain==def[3],]
      ai <- match(these,a$id.match); bi <- match(these,b$id.match)
      if (anyNA(ai) || anyNA(bi) || anyDuplicated(a$id.match) || anyDuplicated(b$id.match)) {
        stop(label,' chain pairing mismatch: ',def[1])
      }
      idx <- match(these,cell$id.match)
      cell$v1[idx] <- a$v_gene[ai]; cell$j1[idx] <- a$j_gene[ai]
      cell$v2[idx] <- b$v_gene[bi]; cell$j2[idx] <- b$j_gene[bi]
      cell$cdr3_1[idx] <- a$cdr3_nt[ai]; cell$cdr3_2[idx] <- b$cdr3_nt[bi]
      cell$clone_key[idx] <- paste(def[1],
        paste0(def[2],':',a$v_gene[ai],':',a$j_gene[ai],':',a$cdr3_nt[ai]),
        paste0(def[3],':',b$v_gene[bi],':',b$j_gene[bi],':',b$cdr3_nt[bi]),sep='|')
    }
    j <- match(cell$id.match,pbmc$id.match)
    cell$sample_id <- pbmc$orig.ident[j]
    cell$severity <- pbmc$severity[j]
    cell$severity2 <- pbmc$severity2[j]
    # Receptor records not matching RNA are logged, but never silently assigned.
    cell$sample_id[is.na(j)] <- d$sample_id[match(cell$id.match,d$id.match)][is.na(j)]
    names2 <- parse_sample(cell$sample_id)
    cell$donor <- names2$donor;cell$visit_index <- names2$visit_index
    say(label,': ',n_orig,' contigs; ',length(all_cells),' cells; ',
        sum(cell$rna_match),' matched PBMC; ',sum(cell$pair_status=='single_pair'),' uniquely paired')
    qc <- data.frame(receptor=label,contigs=n_orig,contig_cells=length(all_cells),
                     rna_matched_contig_cells=sum(cell$rna_match),
                     eligible_contigs=sum(d$eligible),bad_sequence=sum(!valid_seq),
                     invalid_gene=sum(!valid_vj),unresolved_flag_passing_cells=length(unresolved),
                     single_pair_cells=sum(cell$pair_status=='single_pair'),
                     subset_valid_pairs=sum(cell$pair_status=='single_pair' & cell$in_author_subset & cell$rna_match))
    list(cells=cell,qc=qc,flags=flags_report,invalid=unresolved,
         eligible_ids=unique(d$id.match[d$eligible]))
  }
  
  if (!reuse) {
    say('Reading original metadata; linking by canonical id.match (sample + barcode)...')
    pbmc <- metadata('pbmc')
    tmeta <- metadata('T_cell');bmeta <- metadata('B_cell')
    for (subset in list(tmeta,bmeta)) {
      ii <- match(subset$id.match,pbmc$id.match)
      if (anyNA(ii) || !identical(as.character(subset$orig.ident),as.character(pbmc$orig.ident[ii]))) {
        stop('Author subset identity mapping failed; do not use row-name matching')
      }
    }
    t <- receptor('TCR','T',pbmc,tmeta$id.match)
    b <- receptor('BCR','B',pbmc,bmeta$id.match)
    qc <- rbind(t$qc,b$qc)
    flags_report <- rbind(t$flags,b$flags)
    save_csv(qc,'9_receptor_qc.csv')
    save_csv(flags_report,'9_contig_flag_audit.csv')
    # Do not count the same RNA barcode as an independently pure T and B cell.
    for (lbl in c('TCR','BCR')) {
      obj <- if(lbl=='TCR') t else b
      other <- if(lbl=='TCR') b else t
      z <- obj$cells
      ok <- z$pair_status=='single_pair' & z$rna_match & z$in_author_subset &
            !(z$id.match %in% other$eligible_ids)
      if (lbl=='TCR') ok <- ok & z$receptor_subtype=='TCR_ab'
      z$primary_pair <- ok
      z$clone_id <- NA_character_
      keys <- paste(z$donor[ok],z$clone_key[ok],sep='||')
      uniq <- sort(unique(keys))
      z$clone_id[ok] <- sprintf('%s_%07d',lbl,match(keys,uniq))
      if (lbl=='TCR') t$cells <- z else b$cells <- z
      say(lbl,': primary paired RNA-lineage cells (excluding joint chain evidence) = ',sum(ok))
    }
    # Count a clonotype separately for each collection; clone ID is donor-scoped.
    make_clones <- function(z,label) {
      keep <- z[z$primary_pair,c('id.match','sample_id','donor','visit_index','severity','severity2',
                               'clone_id','clone_key','v1','j1','v2','j2','cdr3_1','cdr3_2')]
      if (!nrow(keep)) stop('No primary ',label,' pairs survived QC')
      size <- aggregate(rep(1L,nrow(keep)),keep[c('sample_id','clone_id')],sum)
      names(size)[3] <- 'clone_cells';size$receptor <- label
      n_by_sample <- aggregate(size$clone_cells,list(sample_id=size$sample_id),sum)
      names(n_by_sample)[2] <- 'total_paired_cells'
      size <- merge(size,n_by_sample,by='sample_id',sort=FALSE)
      size$fraction <- size$clone_cells/size$total_paired_cells
      names2 <- parse_sample(size$sample_id)
      size$donor <- names2$donor;size$visit_index <- names2$visit_index
      size$size_bin <- as.character(cut(size$clone_cells,c(0,1,5,20,100,Inf),
                    labels=c('1','2-5','6-20','21-100','>100'),right=TRUE))
      catalog <- keep[!duplicated(keep$clone_id),c('clone_id','donor','clone_key','v1','j1','v2','j2','cdr3_1','cdr3_2')]
      list(cells=keep,sizes=size,catalog=catalog)
    }
    tc <- make_clones(t$cells,'TCR');bc <- make_clones(b$cells,'BCR')
    state <- list(TCR=tc,BCR=bc,t_cells=t$cells,b_cells=b$cells,
                  qc=qc,flags=flags_report,pbmc_samples=unique(pbmc[c('orig.ident','severity','severity2')]))
    say('Saving receptor checkpoint BEFORE figures...')
    saveRDS(state,checkpoint,compress='gzip')
    say('Checkpoint saved: ',checkpoint)
    rm(pbmc,tmeta,bmeta,t,b,tc,bc);invisible(gc(FALSE))
  } else {
    if (!file.exists(checkpoint)) stop('No Step 9 checkpoint to reuse')
    say('Reusing existing Step 9 receptor checkpoint...')
    state <- readRDS(checkpoint)
  }

  # Clone-size fractions within each sample; no cell-level pseudoreplication tests.
  sizes <- rbind(state$TCR$sizes,state$BCR$sizes)
  samples <- state$pbmc_samples
  names(samples)[1] <- 'sample_id'
  samples <- unique(samples)
  samples <- samples[order(samples$sample_id),]
  sample_qc <- do.call(rbind,lapply(c('TCR','BCR'),function(lbl) {
    dat <- state[[lbl]]$sizes
    do.call(rbind,lapply(seq_len(nrow(samples)),function(i) {
      si <- samples$sample_id[i];d <- dat[dat$sample_id==si,]
      denom <- sum(d$clone_cells)
      data.frame(receptor=lbl,sample_id=si,severity=samples$severity[i],
        severity2=samples$severity2[i],
        donor=parse_sample(si)$donor,visit_index=parse_sample(si)$visit_index,
        paired_cells=denom,unique_clones=nrow(d),
        multicell_clone_cells=sum(d$clone_cells[d$clone_cells>=2]),
        fraction_multicell=if(denom>0) sum(d$clone_cells[d$clone_cells>=2])/denom else NA_real_,
        largest_clone=if(nrow(d)) max(d$clone_cells) else 0L,
        largest_clone_fraction=if(denom>0) max(d$clone_cells)/denom else NA_real_,
        simpson_concentration=if(denom>0) sum(d$fraction^2) else NA_real_)
    }))
  }))
  save_csv(sample_qc,'9_sample_repertoire_summary.csv')
  save_gz(sizes,'9_clone_sizes_by_sample.csv.gz')
  save_gz(rbind(data.frame(receptor='TCR',state$TCR$catalog),
                data.frame(receptor='BCR',state$BCR$catalog)), '9_clone_catalog.csv.gz')
  gg <- ggplot2::ggplot
  theme <- ggplot2::theme_classic
  pngout <- function(p,name,w=13,h=7) {
    ggplot2::ggsave(file.path(od,name),p,width=w,height=h,units='in',dpi=170,limitsize=FALSE)
    say('Figure: ',name)
  }
  # Main plot uses CELL-WEIGHTED proportions, not clone counts.
  class_cells <- aggregate(sizes$clone_cells,sizes[c('sample_id','receptor','size_bin')],sum)
  names(class_cells)[4] <- 'n_cells'
  den <- aggregate(class_cells$n_cells,class_cells[c('sample_id','receptor')],sum)
  names(den)[3] <- 'total'
  class_cells <- merge(class_cells,den,by=c('sample_id','receptor'))
  class_cells$percentage <- 100*class_cells$n_cells/class_cells$total
  class_cells$size_bin <- factor(class_cells$size_bin,levels=c('1','2-5','6-20','21-100','>100'))
  save_csv(class_cells,'9_clone_size_bins_by_sample.csv')
  p <- gg(class_cells,ggplot2::aes(x=sample_id,y=percentage,fill=size_bin)) +
       ggplot2::geom_col() + ggplot2::facet_wrap(~receptor,ncol=1) + theme(base_size=11)+
       ggplot2::theme(axis.text.x=ggplot2::element_text(angle=70,hjust=1)) +
       ggplot2::labs(title='Captured cells by paired clonotype size | S1A-related',
                     subtitle='Denominator = successfully paired, RNA-lineage-matched receptor cells',
                     x='Collection',y='Percent of paired cells',fill='Clone size')
  pngout(p,'9_S1A_clone_size_distribution.png',15,9)
  
  # Figure S1B-C concept: V and J usage across TCR beta chains by sample.
  vj <- state$TCR$cells
  count_vj <- function(field) {
    x <- aggregate(rep(1L,nrow(vj)),vj[c('sample_id',field)],sum)
    names(x)[3] <- 'n_cells'
    names(x)[2] <- 'gene'
    x$gene_type <- field
    x
  }
  gene_usage <- rbind(count_vj('v2'),count_vj('j2'))
  gene_usage <- merge(gene_usage,aggregate(gene_usage$n_cells,
                          gene_usage[c('sample_id','gene_type')],sum),by=c('sample_id','gene_type'))
  names(gene_usage)[names(gene_usage)=='x'] <- 'all_paired'
  gene_usage$percent <- 100*gene_usage$n_cells/gene_usage$all_paired
  save_csv(gene_usage,'9_S1_BC_TCR_beta_VJ_gene_usage.csv')
  for (field in c('v2','j2')) {
    xx <- gene_usage[gene_usage$gene_type==field,]
    totals <- sort(tapply(xx$n_cells,xx$gene,sum),decreasing=TRUE)
    xx <- xx[xx$gene %in% head(names(totals),18),]
    xx$gene <- factor(xx$gene,levels=rev(head(names(totals),18)))
    pp <- gg(xx,ggplot2::aes(x=sample_id,y=gene,fill=percent))+
      ggplot2::geom_tile()+ggplot2::scale_fill_gradient(low='white',high='#235b91')+
      theme(base_size=11)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=70,hjust=1))+
      ggplot2::labs(title=paste('TCR beta',if(field=='v2') 'V' else 'J','gene usage | S1B/C-related'),
        subtitle='Top 18 observed genes; fraction of paired TCR-alpha/beta cells in each collection',
        x='Collection',y='Gene',fill='% paired cells')
    pngout(pp,paste0('9_S1_',if(field=='v2')'B_TRBV' else 'C_TRBJ','_usage.png'),15,7)
  }
  
  # S1D: rank top 13 clonotypes at donor's FIRST sample, hold identities fixed
  # and track their captured relative frequencies over that donor's later samples.
  make_tracking <- function(lbl,k=13L) {
    dat <- state[[lbl]]$sizes
    if (!nrow(dat)) return(data.frame())
    selected <- do.call(rbind,lapply(split(dat,dat$donor),function(s) {
      earliest <- min(s$visit_index)
      at_first <- s[s$visit_index==earliest,,drop=FALSE]
      at_first <- at_first[order(-at_first$clone_cells,at_first$clone_id),,drop=FALSE]
      if (!nrow(at_first)) return(NULL)
      at_first <- head(at_first,k)
      at_first$initial_rank <- seq_len(nrow(at_first))
      at_first[c('donor','clone_id','initial_rank')]
    }))
    if (is.null(selected) || !nrow(selected)) return(data.frame())
    sampled <- sample_qc[sample_qc$receptor==lbl,c('sample_id','donor','visit_index','paired_cells')]
    # Each donor's selected clones are evaluated at every available collection,
    # including genuine zeros (unobserved sequences).
    expected <- merge(selected,sampled,by='donor',sort=FALSE)
    lookup <- paste(dat$sample_id,dat$clone_id,sep='|')
    expected$clone_cells <- dat$clone_cells[match(paste(expected$sample_id,expected$clone_id,sep='|'),lookup)]
    expected$clone_cells[is.na(expected$clone_cells)] <- 0L
    expected$frequency_pct <- ifelse(expected$paired_cells>0,
                             100*expected$clone_cells/expected$paired_cells,NA_real_)
    expected$receptor <- lbl
    expected
  }
  track <- rbind(make_tracking('TCR'),make_tracking('BCR'))
  save_csv(track,'9_S1D_top13_donor_clonotype_trajectories.csv')
  # Limit plotting to longitudinal donors; healthy and single-visit donors stay in CSV.
  nvis <- tapply(track$sample_id[track$receptor=='TCR'],track$donor[track$receptor=='TCR'],function(x) length(unique(x)))
  longitudinal <- names(nvis)[nvis>1]
  tr <- track[track$receptor=='TCR' & track$donor %in% longitudinal,]
  if(nrow(tr)) {
    tr$clone_rank <- factor(paste0('Rank ',tr$initial_rank),levels=paste0('Rank ',1:13))
    pp <- gg(tr,ggplot2::aes(x=visit_index,y=frequency_pct,group=clone_id,color=clone_rank))+
        ggplot2::geom_line(linewidth=0.45,alpha=0.65)+ggplot2::geom_point(size=0.75)+
        ggplot2::facet_wrap(~donor,scales='free_y')+theme(base_size=10)+
        ggplot2::scale_x_continuous(breaks=sort(unique(tr$visit_index)))+
        ggplot2::labs(title='TCR top 13 early clonotypes tracked over visits | S1D-related',
          subtitle='Clone IDs fixed at each donor first recorded collection (not antigen specificity)',
          x='Ordinal collection number (not verified weeks)',y='% of paired TCR cells',color='Early rank')
    pngout(pp,'9_S1D_TCR_top13_longitudinal.png',16,10)
  }
  
  # Gene-pair frequencies are descriptive approximations of S5C, not a Circos plot.
  tcells <- state$TCR$cells
  tcells$v_j_pair <- paste(tcells$v2,tcells$j2,sep=' / ')
  vu <- aggregate(rep(1L,nrow(tcells)),tcells[c('sample_id','donor','visit_index','v_j_pair')],sum)
  names(vu)[5] <- 'n_paired_cells'
  save_csv(vu,'9_S5C_TCR_beta_VJ_pairs_by_collection.csv')
  
  # Link to sampled Seurat-cluster tables ONLY for retained sampled cells.
  # Clones are defined from FULL receptor workbook, not from a 40k RNA sample.
  for (lbl in c('TCR','BCR')) {
    step <- if(lbl=='TCR') '4' else '5'
    pth <- if(lbl=='TCR') file.path(root,'4_results','4_tcell_seurat_processed_sampled.rds') else
                            file.path(root,'5_results','5_bcell_seurat_processed_sampled.rds')
    if(!file.exists(pth)) {
      say('WARNING: ',lbl,' sampled Seurat checkpoint absent. Cannot make receptor/RNA cluster linkage: ',pth)
      next
    }
    say('Linking ',lbl,' paired receptors to Step ',step,' sampled Seurat clusters...')
    obj <- readRDS(pth)
    md <- obj@meta.data
    if (!all(c('id.match','analysis_cluster') %in% names(md))) {
      # T-cell Step4 may use analysis_clusters; name it without guesswork.
      if (all(c('id.match','analysis_clusters') %in% names(md))) {
        md$analysis_cluster <- md$analysis_clusters
      } else stop('Sampled ',lbl,' Seurat object lacks id.match or cluster labels')
    }
    if (anyNA(md$id.match) || anyDuplicated(md$id.match)) stop('Duplicate sampled Seurat canonical IDs')
    matched <- state[[lbl]]$cells
    matched <- matched[matched$primary_pair,c('id.match','sample_id','clone_id')]
    if (!nrow(matched)) next
    matched$sample_clone_cells <- state[[lbl]]$sizes$clone_cells[
      match(paste(matched$sample_id,matched$clone_id,sep='|'),
            paste(state[[lbl]]$sizes$sample_id,state[[lbl]]$sizes$clone_id,sep='|'))]
    joined <- merge(data.frame(id.match=as.character(md$id.match),
                  cluster=as.character(md$analysis_cluster),stringsAsFactors=FALSE),matched,
                  by='id.match',all=FALSE)
    joined$clone_category <- cut(joined$sample_clone_cells,c(0,1,5,Inf),
                                 labels=c('singleton','2-5 cells','6+ cells'))
    cluster_table <- as.data.frame(table(cluster=joined$cluster,clone_category=joined$clone_category))
    cluster_table <- cluster_table[cluster_table$Freq>0,]
    names(cluster_table)[3] <- 'paired_cells'
    denom <- aggregate(cluster_table$paired_cells,list(cluster=cluster_table$cluster),sum)
    names(denom)[2] <- 'cluster_paired_total'
    cluster_table <- merge(cluster_table,denom,by='cluster')
    cluster_table$percent <- 100*cluster_table$paired_cells/cluster_table$cluster_paired_total
    save_csv(cluster_table,paste0('9_',lbl,'_clone_sizes_by_sampled_RNA_cluster.csv'))
    say(lbl,': ',nrow(joined),' primary paired receptor cells among sampled RNA cells')
    if(nrow(cluster_table)) {
      pp <- gg(cluster_table,ggplot2::aes(x=cluster,y=percent,fill=clone_category))+
        ggplot2::geom_col()+theme(base_size=11)+
        ggplot2::labs(title=paste(lbl,'paired-clone sizes across sampled RNA clusters'),
          subtitle='Full-workbook clone sizes; denominator = sampled RNA cells with primary paired receptor',
          x='De novo sampled RNA cluster',y='Percent of paired cells',fill='Clone size')
      pngout(pp,paste0('9_',lbl,'_clone_expansion_by_RNA_cluster.png'),13,6)
    }
    rm(obj,md,joined);invisible(gc(FALSE))
  }
  save_csv(state$qc,'9_receptor_qc.csv')
  save_csv(state$flags,'9_contig_flag_audit.csv')
  writeLines(c(
    'Source: Kim et al. 2023, Tracking antigen-specific TCR clonotypes in SARS-CoV-2 infection reveals distinct severity trajectories.',
    'GEO accession: GSE196198; user-provided supplementary figure captions and methods were consulted.',
    'S1A: clone size distribution by collection; S1B-C: TCR beta V/J usage; S1D: top 13 early clones tracked by donor.',
    'S5C: donor-collection TCR beta V-J gene-pair frequencies are supplied as tables only; NOT the authors Circos graph.',
    'S7G: BCR somatic hypermutation lineage reconstruction NOT reproduced here. Requires germline alignment and clone-family inference.',
    'TCR primary: exactly one distinct TRA and one TRB chain, including nucleotide CDR3 and V/J identity, with RNA T-subset membership.',
    'BCR primary: exactly one IGH and IGK/IGL chain, with supplied QC flags and RNA B-subset membership.',
    'Gamma/delta and ambiguous/multiple-chain cells excluded from primary TCR-alpha/beta analysis; joint TCR/BCR chain evidence excluded.',
    'Unique donor-scoped clone identity across visits uses exact chain type, V gene, J gene and CDR3 nucleotide sequence.',
    'TCR productivity is unverified if absent from the TCR workbook; do not infer experimental antigen specificity.',
    'Basic CDR3 nucleotide/amino-acid length and alphabet checks, not full contig sequence verification.',
    'All clone sizes describe captured cells, not absolute in-vivo clonal abundance.',
    'Donor and visit indices inferred from sample names; ordinal visits are not biological days/weeks.',
    'Denominator for frequency: primary paired cells in the same collection; not all PBMC cells.',
    'No donor-group inferential tests, DE p-values or antigen-specificity claims in this exploratory receptor stage.',
    'The Seurat/RNA link is limited to the existing 40,000-cell sampled T/B objects, unlike the full receptor-workbook summary.',
    'Missing or unpaired receptor cells should not be equated with absent antigen-specific immunity.'
  ),file.path(od,'9_methods_and_paper_alignment.txt'))
  writeLines(capture.output(sessionInfo()),file.path(od,'9_session_info.txt'))
  err <- file.path(od,'9_error.txt'); if(file.exists(err)) unlink(err)
  zipfile <- file.path(root,'9_results.zip')
  share <- list.files(od,full.names=TRUE)
  share <- share[!grepl('9_receptor_checkpoint\\.rds$',share)]
  zip <- Sys.which('zip')
  if(!nzchar(zip)) stop('System zip not available; install zip')
  if(file.exists(zipfile)) unlink(zipfile)
  status <- suppressWarnings(system2(zip,c('-q','-9',shQuote(zipfile),shQuote(share))))
  if(!identical(as.integer(status),0L) || !file.exists(zipfile)) stop('Could not create 9_results.zip')
  say('Created upload ZIP: ',zipfile)
  say('SUCCESS | Step 09 completed; keep 9_receptor_checkpoint.rds on your computer')
}

tryCatch(run(), error=function(e) {
  msg <- paste0('STEP 09 FAILED: ',conditionMessage(e))
  cat(msg,'\n',file=stderr())
  args <- commandArgs(trailingOnly=TRUE)
  base <- if(length(args)) args[[1]] else '.'
  dir.create(file.path(base,'9_results'),recursive=TRUE,showWarnings=FALSE)
  cat(msg,'\n',file=file.path(base,'9_results','9_error.txt'))
  quit(save='no',status=1)
})
