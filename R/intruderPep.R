# library(EnsDb.Hsapiens.v86)
# library(Biostrings)
# library(ensembldb)


intruderPep <- function(tx.df,blat.gr,
                        file.prefix,out.dir,
                        use_hash_table = FALSE,
                        diamond.df = diamond.df,

                        # merged_tx_only = TRUE,
                        polyA_threshold = 15,
                        min_length_novel = 2,
                        kmer_size = 5,
                        blat.filter = FALSE){
   # defaults: use_hash_table = FALSE; polyA_threshold = 15; min_length_novel = 2; kmer_size = 5;

   sid <- file.prefix
   np.dir <- file.path(out.dir,'prot')
   neopep_fa_file <- file.path(np.dir,paste(file.prefix,'neopeptides.fa',sep='_'))
   neopep_csv_file <- file.path(np.dir,paste(file.prefix,'neopeptides.csv',sep='_'))
   orf_csv_file <- file.path(np.dir,paste(file.prefix,'_orfs.csv',sep=''))

   tx_log_filename <- file.path(np.dir,paste(file.prefix,'TX_to_NP_log.csv',sep='_'))

   if (!dir.exists(np.dir)){dir.create(np.dir)}
   if (kmer_size + min_length_novel >= 9){
      message('WARNING: setting kmer_size to greater than 9-min_length_novel may cause overestimation of neopeptide sequence.')
   }

   event_list <- sort(unique(tx.df$DNA_eventId))
   tx_log <- character(nrow(tx.df))

   rm_nested = TRUE

   orf.df <- data.frame(stringsAsFactors = FALSE)
   reads.df <- data.frame(stringsAsFactors = FALSE)
   np.ss <- AAStringSet()
   np.df <- data.frame(stringsAsFactors = FALSE)




   # Begin Event Loop ----

   tic()
   for (event in event_list){

      e_tx.gr <- blat.gr[grepl(paste('J',event,'',sep = '_'),blat.gr$query)]
      e_tx.df <- tx.df[tx.df$DNA_eventId == event,]
      tx_idx <- which(tx.df$DNA_eventId == event)

      gene_list <- unique(x = unlist(strsplit(
         c(e_tx.df$geneA[e_tx.df$geneA != 'no_gene'],e_tx.df$geneB[e_tx.df$geneB != 'no_gene']),
         ';')))

      message(paste(
         'Event:',event,'\tnTX:',nrow(e_tx.df),
         '\tRPs:',unique(e_tx.df$RNA_RPs),
         '\tChrA',unique(e_tx.df$chrA),
         '\tChrB',unique(e_tx.df$chrB),
         '\tGene(s):',paste(gene_list,collapse=', ')
      ))


      # Read protein reference fasta ----

      e_prot_ref_fa_filename <- file.path(np.dir,'ref',paste('event',event,'prot_ref.fa',sep='_'))

      if (!file.exists(e_prot_ref_fa_filename)){
         gene_prot <- AAStringSet()
      } else {
         gene_prot <- readAAStringSet(e_prot_ref_fa_filename)
      }

      check_for_gene_prot <- length(gene_prot) > 0

      if (check_for_gene_prot & use_hash_table){
         hash_file <- file.path(np.dir,'ref',paste('event',event,'prot_ref_hash.rds',sep='_'))
         hash.tab <- readRDS(file = file.path(hash_file))
      }

      if (!check_for_gene_prot & length(e_tx.gr) > 0){
         suppressWarnings(cds.grl <- cdsBy(x = edb.ucsc,by = 'gene',
                                           filter = GRangesFilter(value = c(e_tx.gr,invertStrand(e_tx.gr)))))
         if (length(cds.grl) > 0){
            cds.gr <- subsetByOverlaps(x = unlist(cds.grl),ranges = e_tx.gr,ignore.strand = TRUE)
         } else {
            cds.gr <- GRanges()
         }

         filter_list <- AnnotationFilterList(TxBiotypeFilter('protein_coding'),
                                             GRangesFilter(value = c(e_tx.gr,invertStrand(e_tx.gr))))

         exons.grl <- exonsBy(x = edb.ucsc,by = 'tx',
                              filter = filter_list)
         if (length(exons.grl) > 0){
            exons.gr <- subsetByOverlaps(x = unlist(exons.grl),ranges = e_tx.gr,ignore.strand = TRUE)
         } else {
            exons.gr <- GRanges()
         }

         suppressWarnings(utr5.grl <- fiveUTRsByTranscript(x = edb.ucsc,
                                                           filter = filter_list))
         if (length(utr5.grl) > 0){
            utr5.gr <- subsetByOverlaps(x = unlist(utr5.grl),ranges = e_tx.gr,ignore.strand = TRUE)
         } else {
            utr5.gr <- GRanges()
         }

         suppressWarnings(utr3.grl <- threeUTRsByTranscript(x = edb.ucsc,
                                                            filter = filter_list))
         if (length(utr3.grl) > 0){
            utr3.gr <- subsetByOverlaps(x = unlist(utr3.grl),ranges = e_tx.gr,ignore.strand = TRUE)
         } else {
            utr3.gr <- GRanges()
         }

         suppressWarnings(promoter.grl <- promoters(x = edb.ucsc,
                                                    filter = filter_list))

         if (length(promoter.grl) > 0){
            promoter.gr <- subsetByOverlaps(x = promoter.grl,ranges = e_tx.gr,ignore.strand = TRUE)
         } else {
            promoter.gr <- GRanges()
         }


         if (length(cds.gr) > 0){
            msg = 'CDS hit; no reference protein.'
            message(paste('Event',event,'has no reference protein but hits a CDS!'))
            tx_log[tx_idx] <- rep(msg,length(tx_idx))
         } else if (length(exons.gr) > 0) {
            msg = 'Exon hit; no CDS.'
            message(paste('Event',event,'has no reference protein but hits an exon!'))
            tx_log[tx_idx] <- rep(msg,length(tx_idx))
         } else if (length(utr5.gr) > 0){
            msg = '5pUTR hit; no CDS.'
            message(paste('Event',event,'has no reference protein but hits a 5pUTR!'))
            tx_log[tx_idx] <- rep(msg,length(tx_idx))
         } else if (length(utr3.gr) > 0){
            msg = '3pUTR hit; no CDS.'
            message(paste('Event',event,'has no reference protein but hits a 3pUTR!'))
            tx_log[tx_idx] <- rep(msg,length(tx_idx))
         } else if (length(promoter.gr) > 0){
            msg = 'promoter hit; no CDS.'
            message(paste('Event',event,'has no reference protein but hits a promoter!'))
            tx_log[tx_idx] <- rep(msg,length(tx_idx))
         } else {
            msg = 'No functional hits'
            message(paste('Event',event,'hits no protein-impacting regions.'))
            tx_log[tx_idx] <- rep(msg,length(tx_idx))
         }
         if (nrow(e_tx.df) >= quantile(table(tx.df$DNA_eventId),probs = .25)){
            captainsLog(sid = sid,eid = event,tid = NA,msg = msg)
         }

      } else if (length(e_tx.gr) == 0){
         message(paste('Event',event,'has no BLAT hits.'))
         tx_log[tx_idx] <- rep('No BLAT hits',length(tx_idx))

      } else {

         ## Transcript Checkpoint ----

         tx_list <- 1:nrow(e_tx.df)

         for (i_tx in tx_list){

            tx_seq = as.character(e_tx.df[i_tx,'tx'])
            tx_name <- as.character(e_tx.df[i_tx,'tx_num'])
            tx_idx <- which(tx.df$tx_num == tx_name)
            tx.gr <- blat.gr[blat.gr$query == tx_name,]

            i_tx_diamond.df <- diamond.df[diamond.df$query == tx_name,]

            input_diamond.df <- data.frame()
            filter_gene_prot <- FALSE
            if (filter_gene_prot){
               for (g in unique(i_tx_diamond.df$gene_name)){
                  g_tx_diamond.df <- i_tx_diamond.df[which(i_tx_diamond.df$gene_name == g),]

                  if (nrow(g_tx_diamond.df) > 1 & length(unique(g_tx_diamond.df$gene_name)) == 1){
                     g_tx_diamond.df <- g_tx_diamond.df[g_tx_diamond.df$bit_score == max(g_tx_diamond.df$bit_score),]
                  }
                  if (nrow(g_tx_diamond.df) > 1 & any(duplicated(gene_prot[names(gene_prot) %in% g_tx_diamond.df$gp_match]))){
                     g_tx_diamond.df <- g_tx_diamond.df[g_tx_diamond.df$gp_match %in% names(gene_prot)[!duplicated(gene_prot[names(gene_prot) %in% g_tx_diamond.df$gp_match])],]
                  }
                  if (nrow(input_diamond.df) > 0){input_diamond.df <- c(input_diamond.df,g_tx_diamond.df)} else {input_diamond.df <- g_tx_diamond.df}
               }
               tx_diamond_match <- input_diamond.df$gp_match
            } else {
               tx_diamond_match <- i_tx_diamond.df$gp_match
            }
            #tx_diamond_match <- input_diamond.df$gp_match
            message(paste('Checking transcript',e_tx.df[i_tx,'tx_num'],'...'))
            check_no_gaps <- !grepl('N',tx_seq)
            #check_no_poly_A <- !grepl(paste(rep('A',polyA_threshold),collapse=''),tx_seq) & !grepl(paste(rep('T',polyA_threshold),collapse=''),tx_seq)
            check_blat <- length(tx.gr) > 0



            if (check_no_gaps &&  check_blat){

               suppressWarnings(fivePrimeUTR.gr <- get5pUTRbyOverlap(tx.gr))
               suppressMessages(breakpoint.gr <- getBreakpointByOverlap(tx.gr))



               plot_indiv_tx <- FALSE
               if (plot_indiv_tx){
                  #
                  # plot region A
                  A.gr <- GRanges(seqnames = e_tx.df[i_tx,'chrA'],ranges = IRanges(start = e_tx.df[i_tx,'posAmin'],end = e_tx.df[i_tx,'posAmax']))
                  A.gviz <- suppressWarnings(plotJunctionTX(tx.gr,pos.gr = A.gr,highlight = TRUE))
                  #View(A.gviz$gvizObj)


                  # plot region B
                  B.gr <- GRanges(seqnames = e_tx.df[i_tx,'chrB'],ranges = IRanges(start = e_tx.df[i_tx,'posBmin'],end = e_tx.df[i_tx,'posBmax']))
                  B.gviz <- plotJunctionTX(tx.gr = tx.gr,pos.gr = B.gr,highlight = TRUE)

               } # close plot case

               if (length(gene_prot) > 0){

                  # Call matchORF ----

                  matchORF_output <- matchORF(tx_seq = tx_seq,
                                              tx_name = tx_name,
                                              #gene_prot = gene_prot[names(gene_prot) %in% tx_diamond_match],
                                              sid = sid,
                                              diamond.df = diamond.df[diamond.df$query == tx_name,],
                                              kmer_size = kmer_size,
                                              #fivePrimeUTR.gr = fivePrimeUTR.gr,
                                              #breakpoint.gr = breakpoint.gr,
                                              min_length_novel = min_length_novel)
                  if (length(matchORF_output[['log']])>1){
                     why()
                  }
                  if (any(grepl(',',matchORF_output[['log']]))){
                     why()
                  }
                  if (length(matchORF_output[['log']]) > 0 && matchORF_output[['log']] != ''){
                     tx_log[tx_idx] <- matchORF_output[['log']]
                  } else if (length(matchORF_output[['log']]) > length(tx_idx)){
                     tx_log[tx_idx] <- 'ERROR: Catch case of multiple log output'
                  } else {
                     tx_log[tx_idx] <- 'matchORF returned no tx log info'
                  }
                  iORF.df <- matchORF_output[['orf.df']]

                  if(nrow(iORF.df) == 0){
                     tx_log[tx_idx] <- ifelse(tx_log[tx_idx] == '','No matchORF output.',tx_log[tx_idx])
                  } else if  (nrow(iORF.df) >= 1){
                     # iORF.df$gene_5p <- iORF.df$gene
                     if (!'gene_3p' %in% names(iORF.df)){
                        iORF.df$gene_3p <- findThreePrimeUTR(iORF.df = iORF.df, tx.gr = tx.gr)
                     }


                     iORF.df$event <- event
                     iORF.df$tx_name <- tx_name
                     iORF.df$tx_conf <- e_tx.df[i_tx,'tx_conf']
                     iORF.df$rps <- e_tx.df[i_tx,'RNA_RPs']
                     iORF.df$gene_prediction <- e_tx.df[i_tx,'gene_prediction']
                     iORF.df$n_reads <- e_tx.df[i_tx,'n_reads']

                     # iORF.df$seq_list <- e_tx.df[i_tx,'seq_list']
                     #                        iORF.df$novel_seq <- novel_seq

                     if (nrow(orf.df) == 0) {
                        orf.df <- as.data.frame(iORF.df)
                     } else {
                        orf.df <- rbind(orf.df,iORF.df)
                     } #
                  } #if more than 1 ORF
               } else if (length(gene_prot) == 0){
                  message('gene_prot is empty')
                  tx_log[tx_idx] <- 'gene_prot is empty'
               } # if |gene_prot| > 0
            } else if (!check_no_gaps){ # Close transcript checkpoint
               message('Transcript has Ns in sequence')
               tx_log[tx_idx] <- 'Transcript has gaps'
            } else if (!check_no_poly_A){
               message('Poly-A segment found.')
               tx_log[tx_idx] <- 'polyA'
            } else if (!check_blat){
               message('Failed BLAT QC.')
               tx_log[tx_idx] <- 'failed blat QC'
            } else {
               message('wut')
               tx_log[tx_idx] <- 'case I forgot to make log for'
            }
            # loop over each TX
         } # close loop: nrow(e_tx.df) > 0
      } # close event-gene-overlap checkpoint
   } # close event loop
   dt <- toc()

   ### Formatting stringset output ----
   if (nrow(orf.df)>0){


      fs_filter <- !orf.df$mid_exon_break_flag | orf.df$breakpoint
      init_filter <- !orf.df$late_start_seq_flag | orf.df$fivepUTR

      # merged_filter <- grepl('merged',orf.df$name)

      orf.df$ex_crit <- character(nrow(orf.df))

      for (i in 1:nrow(orf.df)){
         if (!orf.df[i,'mid_exon_break_flag'] & !orf.df[i,'late_start_seq_flag']){
            orf.df[i,'ex_crit'] <- 'PASS'
         } else if ( (orf.df[i,'breakpoint'] & orf.df[i,'mid_exon_break_flag']) &
                     (orf.df[i,'fivepUTR']   & orf.df[i,'late_start_seq_flag']) ){
            orf.df[i,'ex_crit'] <- 'PASS w bp + 5pUTR exceptions'
         } else if (orf.df[i,'mid_exon_break_flag'] & orf.df[i,'breakpoint'] ){
            orf.df[i,'ex_crit'] <- 'PASS w breakpoint exception'
         } else if (orf.df[i,'late_start_seq_flag'] & orf.df[i,'fivepUTR'] ){
            orf.df[i,'ex_crit'] <- 'PASS w 5pUTR exception'
         } else if (orf.df[i,'late_start_seq_flag'] & !orf.df[i,'fivepUTR']){
            orf.df[i,'ex_crit'] <- 'late start + no 5pUTR explanation'
         } else if (orf.df[i,'mid_exon_break_flag'] & !orf.df[i,'breakpoint']){
            orf.df[i,'ex_crit'] <- 'kmers stop short of exon end; no breakpoint detected'
         } else {
            orf.df[i,'ex_crit'] <- 'unknown'
         }
      }
      # sum(grepl('PASS',orf.df[fs_filter & merged_filter & init_filter,'ex_crit']))
      write.table(x = orf.df,file = orf_csv_file,quote = FALSE,col.names=TRUE,sep=',',row.names = FALSE)



      filtered_orf.df <- orf.df[fs_filter & init_filter,]
      filtered_orf.df$nObvs <- integer(nrow(filtered_orf.df))
      filtered_orf.df$sumConf <- numeric(nrow(filtered_orf.df))

      if (nrow(filtered_orf.df) > 0 && anyDuplicated(filtered_orf.df$neopep)){
         for (i_orf in 1:nrow(filtered_orf.df)){
            filtered_orf.df[i_orf,'nObvs'] <- sum(filtered_orf.df$neopep == filtered_orf.df[i_orf,'neopep'])
            filtered_orf.df[i_orf,'sumConf'] <- sum(filtered_orf.df$tx_conf[filtered_orf.df$neopep == filtered_orf.df[i_orf,'neopep']])
         }
         filtered_orf.df <- filtered_orf.df[!duplicated(filtered_orf.df$neopep),]
      }

   } else {filtered_orf.df <- orf.df}
   #np_duplicate_filter <- duplicated(filtered_orf.df$neopep)

   if (nrow(filtered_orf.df) > 0){
      neopep_names <- character(length = nrow(filtered_orf.df))
      for (i in 1:nrow(filtered_orf.df)){
         source_prot <- sapply(strsplit(unlist(strsplit(filtered_orf.df$source_prot[i],split = '-')),split = '_'),FUN = '[',1)
         if (length(source_prot) > 1){
            geneA <- source_prot[1]
            geneA_name <- strsplit(geneA,':')[[1]][1]
            geneA_split_exon_list <- unlist(strsplit(unlist(strsplit(geneA,':x'))[2],split = '\\.'))
            geneB <- source_prot[2]
            geneB_name <- strsplit(geneB,':')[[1]][1]
            geneB_split_exon_list <- unlist(strsplit(unlist(strsplit(geneB,':x'))[2],split = '\\.'))
            neopep_names[i] <- paste(geneA_name,':x',
                                     geneA_split_exon_list[length(geneA_split_exon_list)],
                                     '-',geneB_name,':x',
                                     geneB_split_exon_list[length(geneB_split_exon_list)],
                                     sep = ''    )


         } else if (length(source_prot) == 1) {
            source_prot_gene_name <- strsplit(source_prot,':')[[1]][1]

            split_exon_list <- unlist(strsplit(strsplit(source_prot,':x')[[1]][2],split = '\\.'))
            neopep_names[i] <- paste(filtered_orf.df$gene[i],':x',split_exon_list[length(split_exon_list)],sep='')
         } else {
            message('error in source_prot formatting (filtered_orf.df)')
         }

      }
      for (p in unique(neopep_names)){
         matching_names <- neopep_names[neopep_names == p]
         if (length(matching_names) > 1){
            neopep_names[neopep_names == p] <- paste(matching_names,1:length(matching_names),sep='_')
         }

      }

      np.df <- data.frame('Name' = neopep_names,
                          'Seq' = filtered_orf.df$neopep,
                          'ORF' = filtered_orf.df$seq,
                          'event' = filtered_orf.df$event,
                          'Gene' = filtered_orf.df$gene,
                          'ProtName' = filtered_orf.df$source_prot,
                          'tx_name' = filtered_orf.df$tx_name,
                          'nReads' = filtered_orf.df$n_reads,
                          'nObvs' = filtered_orf.df$nObvs,
                          'sumConf' = filtered_orf.df$sumConf,
                          stringsAsFactors = FALSE
      )
   } else if (nrow(orf.df) > 0) {
      message('No neopeptides met filter criteria.')
      np.df <- data.frame()
   } else {
      'No matching ORFs found.'
      np.df <- data.frame()
   }
   reads_analysis <- FALSE
   if (reads_analysis){
      for (i_np in 1:nrow(np.df)){

         matching_np.df <- orf.df[orf.df$neopep == np.df[i_np,'Seq'] ,]
         if (nrow(matching_np.df) > 1){
            blat_map_mask <- !logical(length = nrow(matching_np.df))
            for (j in 1:nrow(matching_np.df)){
               gap_o <- mcols(blat.gr[blat.gr$query == matching_np.df$tx_name[j]])[['gap_o']]
               if (any(gap_o > 0)){
                  blat_map_mask[j] <- FALSE
               }

            }
            matching_np.df <- matching_np.df[blat_map_mask,]
            # np.df[i,'reads'] <- paste(sort(unique(unlist(strsplit(matching_np.df$seq_list,' ')))),collapse = ' ')
            np.df[i_np,'n_reads'] <- matching_np.df$n_reads
         } else if (nrow(matching_np.df) == 1){
            # np.df[i,'reads'] <- matching_np.df$seq_list
            np.df[i_np,'n_reads'] <- matching_np.df$n_reads
         }
      }
      dev_dir <- paste(unlist(strsplit(np.dir,'/'))[1:(grep('prot',unlist(strsplit(np.dir,'/')))-1)],collapse='/')
      reads_tx_hash.rdg <- file.path(dev_dir,'tx',paste('event',np.df[i_np,'DNA_eventId'],sep='_'),'vo','read_tx_hash.Rds')
      reads_tx.hash <- readRDS(reads_tx_hash.rdg)
      reads_list <- rownames(reads_tx.hash)[reads_tx.hash[,np.df[i_np,'tx_name']] == 1]
      # reads_list <- unique(unlist(strsplit(paste(as.character(tx.df$seq_list),collapse = ' '),' ')))
      # reads_list <- reads_list[order(reads_list)]
      n_reads <- length(reads_list)
      reads.df <- data.frame('ReadNum' = reads_list,stringsAsFactors = FALSE)
      reads.mat <- matrix(0,nrow = n_reads,ncol = nrow(np.df))
      for (i_np in 1:n_reads){
         reads.df[i_np,'in_np'] <- paste(np.df[grepl(reads.df$ReadNum[i_np],np.df[,'reads']),'Name'],collapse = ' ')
         for (j in 1:nrow(np.df)){
            reads.mat[i_np,j] <- grepl(np.df$Name[j],reads.df[i_np,'in_np'])
         }


      }
   }

   message(paste(nrow(orf.df),'TX found to be peptigenic;',nrow(np.df),'unique neopeptides found.'))



   if (length(np.df)>0){
      np.ss <- AAStringSet(x = np.df$Seq)
      names(np.ss) <- np.df$Name
      writeXStringSet(x = np.ss,filepath = neopep_fa_file)


      np.df$sampleId <- rep(as.character(file.prefix),nrow(np.df))

      for (i in 1:nrow(np.df)){
         tx_idx <-which(tx.df$tx_num == np.df[i,'tx_name'])
         np.df[i,'chrA'] <- unique(tx.df[tx_idx,'chrA'])
         np.df[i,'chrB'] <- unique(tx.df[tx_idx,'chrB'])
         np.df[i,'blat'] <- unique(tx.df[tx_idx,'blat'])
         if (unique(tx.df$chrA[tx_idx]) == unique(tx.df$chrB[tx_idx])){
            np.df[i,'intra.chromosomal.flag'] <- TRUE
         } else {
            np.df[i,'intra.chromosomal.flag'] <- FALSE
         }
      }

      write.table(x = np.df,file = neopep_csv_file,quote = FALSE,col.names=TRUE,sep=',',row.names = FALSE)
   }

   output <- list('orf.df' = orf.df,
                  'np.df' = np.df,
                  'tx_log' = tx_log)
   return(output)

}

findORF <- function(tx_seq,tx_name,diamond.df,blat.gr,np.dir,
                    np.log.file = NA,tx.blacklist = '',
                    out.file = NA,jct.gr = NA,kmer_cutoff = 5,kmer_size = 5,
                    min_length_novel = 1,overwrite_np = FALSE){
   if (is.na(np.log.file) ){
      np.log.file <- file.path(np.dir,'np.out')
   }
   if (file.exists(out.file) & !overwrite_np){
      return(NULL)
   }

   # edb <- AnnotationDbi::loadDb(file = dbfile(dbconn(EnsDb.Hsapiens.v105)))
   edb <- EnsDb(x = system.file('extdata/EnsDb.Hsapiens.v105.sqlite',package='EnsDb.Hsapiens.v105'))
   # edb <- useMySQL(x = EnsDb.Hsapiens.v105, user = '',pass = '',host = 'localhost')
   on.exit(RSQLite::dbDisconnect(dbconn(edb)))

   seqlevelsStyle(edb) <- "UCSC"


   # jct.gr = NA;kmer_cutoff = 5;kmer_size = 5;min_length_novel = 1
   method = 'blat' # or 'diamond'

   if (!exists('kmer_cutoff')){kmer_cutoff = 5; kmer_size = 5; min_length_novel = 1;}

   output <- list('log' = '','orf.df' = data.frame())
   # blat checkpoint ----
   if (!tx_name %in% blat.gr$query){
      output[['log']] <- 'no BLAT hit'
      cat(paste(tx_name,output[['log']]),file = np.log.file,sep='\n',append = TRUE)
      return(output)
   }

   # ORFfinder ----
   orf_prediction <- orfFinderJS(checkSeq = as.character(tx_seq))
   if (is.null(orf_prediction)){
      output[['log']] <- 'no ORFs found'
      cat(paste(tx_name,output[['log']]),file = np.log.file,sep='\n',append = TRUE)
      return(output)
   }


   if (!is.null(diamond.df) && any(diamond.df$query == tx_name)){
      tx_inDiamond <- sapply(X = 1:nrow(orf_prediction), USE.NAMES = FALSE,
                             FUN = function(i) any( orf_prediction[i,'startPos'] == diamond.df[diamond.df$query == tx_name,'query_start'] &
                                                       orf_prediction[i,'endPos'] == diamond.df[diamond.df$query == tx_name,'query_end']) &&
                                diamond.df[diamond.df$query == tx_name & diamond.df$query_start == orf_prediction[i,'startPos'] & diamond.df$query_end == orf_prediction[i,'endPos'],'percent_id'] >= 99)

      if (any(tx_inDiamond)){
         output[['log']] <- 'Perfect prot match'
         cat(paste(tx_name,output[['log']]),file = np.log.file,sep='\n',append = TRUE)
         return(output)
      }
   }




   ## Remove nested ORFs  ----
   rm_nested = FALSE
   if (rm_nested){
      nested_list <- logical(length = nrow(orf_prediction))
      for (i in 1:nrow(orf_prediction)){
         nested_list[i] <- try(any(grepl(orf_prediction[i,'seq'],orf_prediction[-i,'seq'])))
      }
      orf_prediction <- orf_prediction[which(!nested_list),]
      row.names(orf_prediction) <- 1:nrow(orf_prediction)
   }

   # Identify matching protein(s) ----
   if (method == 'diamond'){
      genes_in_range <- genes(x = edb,filter = GRangesFilter(jct.gr + 5000))

      tx_diamond.df <- diamond.df[diamond.df$query == tx_name & diamond.df$gene_name %in% genes_in_range$gene_name,]
      if (nrow(tx_diamond.df) == 0){
         log <- 'no matching proteins'
         return(list('orf.df' = data.frame(), 'log' = log))
      }
      tx_diamond.df[tx_diamond.df$query_strand == '-',c('query_start','query_end')] <- tx_diamond.df[tx_diamond.df$query_strand == '-',c('query_end','query_start')]
      ## Select top hit isoform for each gene ----

      tx.gr <- GRanges(seqnames = tx_diamond.df$query,ranges = IRanges(start = tx_diamond.df$query_start,end = tx_diamond.df$query_end),strand = tx_diamond.df$query_strand)
      tx.gr$bit_score = tx_diamond.df$bit_score; tx.gr$query_length <- tx_diamond.df$query_length; tx.gr$gp_match<- tx_diamond.df$gp_match

      #message(paste('tx.gr has length',length(tx.gr)))

      tx.range <- range(tx.gr-5,ignore.strand = FALSE)
      #message(paste('tx.range has length',length(tx.range)))
      top_list <- sapply(X = 1:length(tx.range),FUN = function(i) {
         tx.gr.ss <- tx.gr[tx.gr %over% tx.range[i]]
         tx_hit <- tx.gr.ss$gp_match[tx.gr.ss$bit_score == max(tx.gr.ss$bit_score)]
         return(tx_hit)
      }
      ,USE.NAMES = FALSE)
      top.gr <- tx.gr[tx.gr$gp_match %in% top_list,]
      top_diamond.df <- tx_diamond.df[tx_diamond.df$gp_match %in% top.gr$gp_match,]

      # top_diamond.df <- tx_diamond.df[tx_diamond.df$gp_match %in%
      #                                    sapply(X = unique(tx_diamond.df$gene_name),
      #                                           FUN = function(g) tx_diamond.df[tx_diamond.df$gene_name == g & tx_diamond.df$bit_score == max(tx_diamond.df$bit_score[tx_diamond.df$gene_name == g]),'gp_match'],USE.NAMES = FALSE),]
      # top_diamond.df[top_diamond.df$query_strand == '-',c('query_start','query_end')] <- top_diamond.df[top_diamond.df$query_strand == '-',c('query_end','query_start')]
      top_diamond.df <- top_diamond.df[!duplicated(top_diamond.df$target_seq) & !duplicated(top_diamond.df$trans_query_seq),]
      intruderPep_input.df <- top_diamond.df
      gene_prot <- AAStringSet(x = intruderPep_input.df$target_seq); names(gene_prot) <- intruderPep_input.df$gp_match
   } else if (method == 'blat'){




      e_blat.gr <- blat.gr[blat.gr$query == tx_name,]
      e_prot.df <- proteins(edb,
                            filter = AnnotationFilterList(GRangesFilter(value = c(e_blat.gr,invertStrand(e_blat.gr))),
                                                          GeneBiotypeFilter(value = 'protein_coding'),
                                                          TxSupportLevelFilter(value = 1)),
                            columns = c('gene_name','tx_biotype','protein_id','protein_sequence','tx_id','seq_name','seq_strand','tx_seq_start','tx_seq_end','tx_support_level'))
      # e_cds.gr <- unlist(cdsBy(x = edb.ucsc,by = 'tx'))
      # e_ol_cds.gr <- subsetByOverlaps(x = e_cds.gr,ranges = e_blat.gr,ignore.strand = TRUE)
      # suppressWarnings(cdsByProt.pa <- unlist(genomeToProtein(x = e_ol_cds.gr,db = edb.ucsc))); cdsByProt.pa <- cdsByProt.pa[which(mcols(cdsByProt.pa)[,'cds_ok'])]
      # tx_list <- unique(mcols(cdsByProt.pa)[,'tx_id'])
      # prot_list <- unique(names(cdsByProt.pa))
      # #protFromCDS.grl <- GRangesList()
      # protFromCDS.grl <- sapply(X = tx_list,FUN = function(i) cdsByProt.pa[mcols(cdsByProt.pa)[,'tx_id'] == i])





      gene_prot <- AAStringSet(x = e_prot.df$protein_sequence); names(gene_prot) <- paste(e_prot.df$protein_id,e_prot.df$gene_name,sep='|')

      # blat.gr headers: prot_id, gene_name, gp_match

   }

   if (length(gene_prot) == 0){
      output[['log']] <- 'No prot found'
      cat(paste(tx_name,output[['log']]),file = np.log.file,sep='\n',append = TRUE)
      return(output)
   }


   # if (nrow(top_diamond.df) > 1){
   #
   #    top.gr <- GRanges(seqnames = top_diamond.df$query,ranges = IRanges(start = top_diamond.df$query_start,end = top_diamond.df$query_end),strand = top_diamond.df$query_strand)
   #    top.gr$gene_name <- top_diamond.df$gene_name; top.gr$bit_score <- top_diamond.df$bit_score
   #    ol.gr <- top.gr[sapply(X = 1:length(top.gr),FUN = function(i) top.gr[i] %over% top.gr[-i])]
   #    to_discard <- ol.gr$gene_name[ol.gr$bit_score < max(ol.gr$bit_score)]
   #    if (length(to_discard) > 0){
   #       top_diamond.df <- top_diamond.df[-which(top_diamond.df$gene_name %in% to_discard),]
   #    }
   #
   # }



   # find ORF matching translated query seq ----
   orf_kmers <- allKmers(x = orf_prediction$seq,k = kmer_size)

   if (length(orf_kmers) < 3){
      # message('No unique kmers; exiting.')
      output[['log']] <- 'no unique kmers'
      cat(paste(tx_name,output[['log']]),file = np.log.file,sep='\n',append = TRUE)
      return(output)
   }


   aa_kmers <- allKmers(x = gene_prot,k = kmer_size,rename.duplicates = TRUE)



   kmer_hit_idx <- matrix()
   kmer_hits_pairs_tmp <- expand.grid(1:nrow(orf_kmers),1:nrow(aa_kmers))
   kmer_hits_pairs <- paste(kmer_hits_pairs_tmp$Var1,kmer_hits_pairs_tmp$Var2,sep=',')
   kmer_hits_pairs_mat <- matrix(kmer_hits_pairs,nrow = nrow(orf_kmers))


   # register(MulticoreParam(workers = 3))
   kmer_hits_out <- foreach(x = kmer_hits_pairs,.combine=c) %do% {sum(names(orf_kmers)[orf_kmers[as.integer(strsplit(x,',')[[1]][1]),] == 1] %in% names(aa_kmers)[aa_kmers[as.integer(strsplit(x,',')[[1]][2]),] == 1])}
   kmer_hits <- matrix(kmer_hits_out,nrow = nrow(orf_kmers))


   if (!any(kmer_hits > 0)){
      # message('No kmer hits.')
      output[['log']] <- 'No kmer hits'
      cat(paste(tx_name,output[['log']]),sep='\n',file = np.log.file,append = TRUE)
      return(output)
   }


   if (!any(kmer_hits >= kmer_cutoff)){
      # message('Insufficient kmer hits.')
      output[['log']] <- 'n_kmers < kmer_cutoff'
      cat(paste(tx_name,output[['log']]),sep='\n',file = np.log.file,append = TRUE)
      return(output)
   }

   gp_hits <- names(gene_prot)# [colSums(kmer_hits) > 0]
   gp_genes <- sapply(X = gp_hits,FUN = function(x) strsplit(x,'\\|')[[1]][2],USE.NAMES = FALSE)
   gene_hits <- unique(gp_genes)

   # associate w gene & number of kmers ----
   orf_prediction$n_kmers <- 0; orf_prediction$gp_match = ''
   orf_prediction$gene <- sapply(X = 1:nrow(orf_prediction), FUN = function(r) ifelse(test = any(kmer_hits[r,] >= kmer_cutoff),
                                                                                      yes = paste(unique(gp_genes[which(kmer_hits[r,] >= kmer_cutoff)]),collapse=','),
                                                                                      no = ''))

   orf_prediction$gp_match <- sapply(X = 1:nrow(orf_prediction),USE.NAMES = FALSE,
                                     FUN = function(r) {
                                        r_genes <- unique(gp_genes[kmer_hits[r,] >= kmer_cutoff])
                                        if (length(r_genes) == 1){
                                           return(gp_hits[kmer_hits[r,] == max(kmer_hits[r,])][1])
                                        } else if (length(r_genes) > 1){
                                           paste(sapply(X = r_genes,FUN = function(g) gp_hits[gp_genes == g][which(kmer_hits[r,gp_genes == g]== max(kmer_hits[r,gp_genes == g]))[1]],USE.NAMES = FALSE),collapse=',')
                                        } else if (length(r_genes) == 0){
                                           return('')
                                        }
                                     }
   )
   # TODO : adapt this so it's compatible with rm_nested = FALSE
   orf_prediction$n_kmers <- sapply(X = 1:nrow(orf_prediction),USE.NAMES = FALSE,
                                    FUN = function(r) ifelse(test = orf_prediction$gp_match[r] == '',yes = 0,
                                                             no = sum(kmer_hits[r,gp_hits %in% strsplit(orf_prediction$gp_match[r],split = ',')[[1]]])))


   #which(orf_prediction$n_kmers / orf_prediction$lengthAA == max(orf_prediction$n_kmers / orf_prediction$lengthAA))
   # orf_prediction$gene <- sapply(X = 1:nrow(orf_prediction),FUN = function(i) ifelse(test = orf_prediction[i,'gp_match'] != '',
   #                                                                                   yes = strsplit(orf_prediction$gp_match[i],'\\|')[[1]][2],
   #                                                                                   no = ''))

   # orf_prediction$prot <- sapply(X = 1:nrow(orf_prediction),FUN = function(i) ifelse(test = orf_prediction[i,'gp_match'] != '',
   #                                                                                   yes = strsplit(orf_prediction$gp_match[i],'\\|')[[1]][1],
   #                                                                                   no = ''))

   orf_match <- orf_prediction[orf_prediction$n_kmers == max(orf_prediction$n_kmers),]
   if (nrow(orf_match) > 1 & length(unique(orf_match$strand)) == 1 & length(unique(orf_match$frame)) == 1){ # pick shortest nested ORF
      orf_match <- orf_match[orf_match$lengthAA == min(orf_match$lengthAA),]
   } else if (nrow(orf_match) > 1){
      output$log <- 'Multiple plausible ORF hits'
      return(output)
      cat(paste(tx_name,output[['log']]),sep='\n',file = np.log.file,append = TRUE)
   }


   #kmer_match$cum_kmers <- sapply(X = kmer_match$orf,FUN = function(x) sum(kmer_match[kmer_match$orf == x,'n_kmers']))

   # pick top orf by number of kmers ----
   if (nrow(orf_match) == 1){

      kmer_match <- as.data.frame(which(kmer_hits >= kmer_cutoff, arr.ind=TRUE)) # The gene and ORFs with sufficient matching kmers
      names(kmer_match) <-  c('orf','prot')

      kmer_match <- kmer_match[kmer_match$orf == as.integer(row.names(orf_match)) & kmer_match$prot %in% which(gp_hits %in% strsplit(orf_match$gp_match,',')[[1]] ),]
      kmer_match$n_kmers <- sapply(X = paste(kmer_match$orf,kmer_match$prot,sep=','),FUN = function(pair_str) kmer_hits_out[kmer_hits_pairs == pair_str])
      kmer_match$prot_id <- sapply(X = kmer_match$prot, FUN = function(i) strsplit(names(gene_prot[i]),'\\|')[[1]][1])
      kmer_match$gene_name <- sapply(X = kmer_match$prot, FUN = function(i) strsplit(names(gene_prot[i]),'\\|')[[1]][2])
      # kmer_match <- kmer_match[sapply(X = unique(kmer_match$gene_name),FUN = function(g) which(kmer_match$gene_name == g & kmer_match$n_kmers == max(kmer_match$n_kmers[kmer_match$gene_name == g]))[1],USE.NAMES = FALSE),]
   } else if (nrow(orf_match) == 0){
      # message('No matching ORFs found.')
      output[['log']] <- 'No WT prot hits'
      # write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
      cat(paste(tx_name,output[['log']]),sep='\n',file = np.log.file,append = TRUE)
      return(output)
   } else {
      output[['log']] <- 'Multiple ORF hits'
      # write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
      cat(paste(tx_name,output[['log']]),sep='\n',file = np.log.file,append = TRUE)
      return(output)
   }


   out.df <- data.frame()
   orf_match$name <- rep(tx_name,nrow(orf_match))

   # change e_blat.gr strand according to orf_match ----

   if (orf_match$strand == '+'){
      e_blat.gr <- e_blat.gr[sort(x = e_blat.gr$query_start,index.return = TRUE)$ix]
   } else {
      # reset: e_blat.gr <- blat.gr[blat.gr$query == tx_name]
      e_blat.gr <- e_blat.gr[sort(x = e_blat.gr$query_start,index.return = TRUE,decreasing = TRUE)$ix]
      e_blat.gr <- invertStrand(e_blat.gr)

      e_blat.gr$rev_query_start <- nchar(tx_seq) - e_blat.gr$query_end + 1
      e_blat.gr$rev_query_end <- nchar(tx_seq) - e_blat.gr$query_start + 1
      e_blat.gr$query_start <- e_blat.gr$rev_query_start;e_blat.gr$query_end <- e_blat.gr$rev_query_end; e_blat.gr$rev_query_end <- NULL; e_blat.gr$rev_query_start <- NULL
   }
   e_blat.gr <- unique(e_blat.gr)


   # annotate e_blat.gr with corresponding gene structure of kmer_match hit(s) ----
   if (nrow(orf_match) > 0 & nrow(kmer_match) > 0){
      exons.gr <- unlist(exonsBy(x = edb,by = 'tx',
                                 filter = AnnotationFilterList(GRangesFilter(value = e_blat.gr),ProteinIdFilter(value = kmer_match$prot_id)),
                                 columns = c('gene_name','exon_id','protein_id','tx_id')))

      genes.gr <- genes(x = edb,
                        filter = AnnotationFilterList(GeneNameFilter(value = unique(exons.gr$gene_name[exons.gr$gene_name != ''])),
                                                      GeneBiotypeFilter(value = 'LRG_gene',not = TRUE)))


      e_blat.gr$gene_match <- sapply(X = 1:length(e_blat.gr),FUN = function(i) genes.gr$gene_name[genes.gr %over% e_blat.gr[i] ])

      e_blat.gr$cds_match <- sapply(X = 1:length(e_blat.gr),FUN = function(i) ifelse(test = overlapsAny(subject = exons.gr,query = e_blat.gr[i],ignore.strand = TRUE),
                                                                                     yes = subsetByOverlaps(x = exons.gr,ranges = e_blat.gr[i],ignore.strand = TRUE)$gene_name,
                                                                                     no = '*'))
      e_blat.gr$exon_rank <- sapply(X = 1:length(e_blat.gr),FUN = function(i) ifelse(test = overlapsAny(subject = exons.gr,query = e_blat.gr[i],ignore.strand = TRUE),
                                                                                     yes = subsetByOverlaps(x = exons.gr,ranges = e_blat.gr[i],ignore.strand = TRUE)$exon_rank,
                                                                                     no = '*'))
   }


   if (nrow(orf_match) == 1 & nrow(kmer_match) == 1){
      # truncation case ----
      if (method == 'diamond'){
         diamond_dat.df <- top_diamond.df[, c('percent_id','alignment_length','mismatches','gap_o','query_start','query_end','target_length','bit_score','prot_id','gene_name','gp_match')]
         orf_dat.df <- orf_match[,c('name','seq','strand','frame','startPos','endPos','lengthNT','lengthAA','n_kmers')]
         out.df <- cbind(orf_dat.df,diamond_dat.df)
      } else if (method == 'blat'){
         #e_blat.df <- as.data.frame(x = e_blat.gr)
         # prot_cds <- cdsBy(x = edb.ucsc,by = 'tx',filter = ProteinIdFilter(value = c(kmer_match$prot_id,'ENSP00000375067')),columns = c(listColumns(x = edb.ucsc,'cds'),'gene_name','exon_rank'))
         # prot_exons <- unlist(exonsBy(x = edb.ucsc,filter = ProteinIdFilter(value = c(kmer_match$prot_id,'ENSP00000375067')),columns = c(listColumns(x = edb.ucsc,'exons'),'gene_name'),by = 'tx'))
         #prot_exons$exon_rank <- sapply(X = prot_exons$exon_id,FUN = function(e) )
         # exons.gr <- unlist(exonsBy(x = edb.ucsc,by = 'tx',
         #                            filter = AnnotationFilterList(GRangesFilter(value = c(e_blat.gr,invertStrand(e_blat.gr))),ProteinIdFilter(value = kmer_match$prot_id)),
         #                            columns = c('gene_name','exon_id','protein_id')))

         out.df <- orf_match
         # out.df$tx_id <- unique(exons.gr$tx_id)
         out.df$exons <- paste(e_blat.gr$exon_rank,collapse=',')
      }
      gene_prot_hit <- gene_prot[names(gene_prot) == paste(kmer_match$prot_id,kmer_match$gene_name,sep='|')]
      out.df$prot <- strsplit(names(gene_prot_hit),'\\|')[[1]][1]

      orf_kmers_list <- names(orf_kmers)[which(orf_kmers[kmer_match[,'orf'],] > 0)]
      aa_kmers_list <- names(aa_kmers)
      overlapping_kmers <- orf_kmers_list[orf_kmers_list %in% aa_kmers_list]
      orf_kmer_loc <- sort(sapply(X = 1:length(overlapping_kmers),FUN = function(i) regexpr(overlapping_kmers[i],out.df$seq)[1],USE.NAMES = FALSE))
      orf_kmer_loc_delta <- diff(orf_kmer_loc)

      wt_kmer_loc <- sort(sapply(X = 1:length(overlapping_kmers),FUN = function(i) regexpr(overlapping_kmers[i],gene_prot_hit)[1],USE.NAMES = FALSE))
      wt_kmer_loc_delta <- diff(wt_kmer_loc)


      ol_start <- min(orf_kmer_loc)
      ol_end <- max(orf_kmer_loc) + kmer_size - 1
      out.df$break_pos <- ol_end
      out.df$leading_aa <- ol_start
      out.df$novel_length <- out.df$lengthAA - ol_end
      if (out.df$novel_length >= min_length_novel){
         out.df$novel_seq <- subseq(out.df$seq,ol_end - 9 + min_length_novel,out.df$lengthAA)
         out.log <- 'truncation'
      } else {
         # message('Insufficient novel sequence.')
         out.df$novel_seq <- ''
         out.log <- 'No novel seq'
         output[['log']] <- out.log
         cat(paste(tx_name,output[['log']]),file = np.log.file,append = TRUE)
         return(output)
         #output[['log']] <- 'No novel seq'
         #return(output)
      }

      # out.log = '1 neopeptide'
   } else if (nrow(orf_match) == 1 & nrow(kmer_match) == 2){



      # fusion case ----
      message(paste('FUSION!',tx_name))
      ## determine orientation ----

      tx_genes <- unique(e_blat.gr$cds_match[!e_blat.gr$cds_match %in% c('','*')])
      gene5 <- tx_genes[1]
      gene3 <- tx_genes[2]
      out.df <- orf_match
      #out.df$n_kmers <- 0

      out.df$gene <- paste(gene5,gene3,sep='-')
      out.df$prot <- paste(strsplit(strsplit(out.df$gp_match,'-')[[1]][1],'\\|')[[1]][1],
                           strsplit(strsplit(out.df$gp_match,'-')[[1]][2],'\\|')[[1]][1],sep='-')
      out.df$exons <- paste(paste(e_blat.gr$exon_rank[e_blat.gr$gene_match == gene5],collapse=','),
                            paste(e_blat.gr$exon_rank[e_blat.gr$gene_match == gene3],collapse=','),sep='-')

      # out.df$tx_id <- paste(unique(exons.gr$tx_id[exons.gr$gene_name == gene5]),unique(exons.gr$tx_id[exons.gr$gene_name == gene3]),sep='-')

      if (!grepl(gene5,strsplit(out.df$gp_match,',')[[1]][1])){
         out.df$gp_match <- paste(strsplit(out.df$gp_match,',')[[1]][c(2,1)],collapse='-')
      } else {
         out.df$gp_match <- sub(pattern = ',',replacement = '-',x = out.df$gp_match)
      }



      out.log = 'fusion'

      #gene_prot_hit <- gene_prot[names(gene_prot) %in% paste(kmer_match$prot_id,kmer_match$gene_name,sep='|')]
      gp_5p <- gene_prot[names(gene_prot) %in% strsplit(out.df$gp_match,'-')[[1]][1]]
      gp_3p <- gene_prot[names(gene_prot) %in% strsplit(out.df$gp_match,'-')[[1]][2]]

      orf_kmers_list <- names(orf_kmers)[which(orf_kmers[unique(kmer_match[,'orf']),] > 0)]
      aa5_kmers_list <- names(aa_kmers)[which(aa_kmers[names(gp_5p),] > 0)]
      aa3_kmers_list <- names(aa_kmers)[which(aa_kmers[names(gp_3p),] > 0)]

      overlapping_kmers_5 <- orf_kmers_list[orf_kmers_list %in% aa5_kmers_list]
      overlapping_kmers_3 <- orf_kmers_list[orf_kmers_list %in% aa3_kmers_list]

      orf_kmer_loc_5 <- sort(sapply(X = 1:length(overlapping_kmers_5),FUN = function(i) regexpr(overlapping_kmers_5[i],out.df$seq)[1],USE.NAMES = FALSE))
      orf_kmer_loc_3 <- sort(sapply(X = 1:length(overlapping_kmers_3),FUN = function(i) regexpr(overlapping_kmers_3[i],out.df$seq)[1],USE.NAMES = FALSE))

      #orf_kmer_loc_delta <- diff(orf_kmer_loc)

      wt5_kmer_loc <- sort(sapply(X = 1:length(overlapping_kmers_5),FUN = function(i) regexpr(overlapping_kmers_5[i],gp_5p)[1],USE.NAMES = FALSE))
      wt3_kmer_loc <- sort(sapply(X = 1:length(overlapping_kmers_3),FUN = function(i) regexpr(overlapping_kmers_3[i],gp_3p)[1],USE.NAMES = FALSE))
      # wt_kmer_loc_delta <- diff(wt_kmer_loc)
      out.df$break_pos <- as.numeric(paste(unique(max(orf_kmer_loc_5) + 4,min(orf_kmer_loc_3) ),sep='.'))
      out.df$leading_aa <- min(orf_kmer_loc_5)

      # need cols: "exons","break_pos","leading_aa","novel_length","novel_seq"

      out.df$novel_length = 2*9 - min_length_novel
      out.df$novel_seq = substr(x = out.df$seq,start = max(orf_kmer_loc_5) + 4 - (9-min_length_novel),stop = min(orf_kmer_loc_3) + (9 - min_length_novel))
      #out.df$tx_id <- unique(row.names(exons.gr))
      # out.df$exons <- paste(exons.gr$exon_rank,collapse=',')





   } else if (nrow(orf_match) > 1 & nrow(top_diamond.df) > 1){
      out.log = 'wut'
   } else if (nrow(orf_match) > 0 & nrow(top_diamond.df) == 0) {
      out.log = paste(nrow(orf_match),'ORF,', nrow(top_diamond.df),'prot')
   } else if (nrow(orf_match) == 0 & nrow(top_diamond.df) > 0) {
      out.log = paste(nrow(orf_match),'ORF,', nrow(top_diamond.df),'prot')
   } else {
      out.log = paste(nrow(orf_match),'ORF,', nrow(top_diamond.df),'prot')
   }

   #out.df$gp_match <- NULL

   #orf_match <- orf_match[orf_match$n_kmers == max(orf_match$n_kmers),]

   #output <- list('gene_prot' = gene_prot,'kmer_match' = kmer_match,'orf_prediction' = orf_match,'aa_kmers' = aa_kmers,'orf_kmers' = orf_kmers[kmer_match[,'orf'],])

   output <- list('log' = out.log,'orf.df' = out.df)
   if (!dir.exists(dirname(out.file))){dir.create(dirname(out.file),recursive = TRUE)}
   write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
   #writeLines(x = log,)



   return(output)
}




matchORF <- function(tx_seq,gene_prot,tx_name,diamond.df,sid,
                     #fivePrimeUTR.gr,hash.tab,breakpoint.gr,
                     kmer_size = 5, rm_nested = TRUE, kmer_cutoff = 5, fusion_margin = 10, min_length_novel = 2 ) {
   output <- list('orf.df' = data.frame(stringsAsFactors = FALSE),
                  'log' = '')
   orf.df <- data.frame()



   if (FALSE){


      #gene_prot <- gene_prot[names(gene_prot) %in% tx_diamond_match]

      if (length(gene_prot) == 0){
         message('No matching proteins.')
         output[['log']] <- 'No matching proteins'
         write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
         return(output)
      }

      if (any(width(gene_prot) <= 9)){
         gene_prot[width(gene_prot) <= 9] <- NULL
         if (length(gene_prot) == 0){
            message('No matching proteins of sufficient length.')
            output[['log']] <- 'gene_prot too short'
            write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
            return(output)
         }
      }

      if (as.character(tx_seq) == as.character(reverseComplement(DNAString(tx_seq)))){
         message('Sequence is a palindrome!')
         output[['log']] <- 'Palindromic TX'
         write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
         return(output)
      }


      if (is.factor(tx_seq)){tx_seq <- as.character(tx_seq)}

      # ORF prediction ----
      orf_prediction <- orfFinderJS(checkSeq = as.character(tx_seq))



      if (is.null(orf_prediction)){
         output[['log']] <- 'no ORFs found'
         message('No ORFs found.')
         write.table(x = out.df,file = out.file,col.names = FALSE,row.names = FALSE,quote = FALSE)
         return(output)
      }


      ## Remove nested ORFs and prot seqs ----

      if (rm_nested){
         nested_list <- logical(length = nrow(orf_prediction))
         for (i in 1:nrow(orf_prediction)){
            nested_list[i] <- try(any(grepl(orf_prediction[i,'seq'],orf_prediction[-i,'seq'])))
         }
         orf_prediction <- orf_prediction[which(!nested_list),]

      }


      ## Identify and tally kmer matches ----
      orf_kmers <- allKmers(x = orf_prediction$seq,k = kmer_size)

      if (length(orf_kmers) < 3){
         message('No unique kmers; exiting.')
         output[['log']] <- 'no unique kmers'
         return(output)
      }

      # there used to be a flag here that controlled for duplicate gene_prot entries, which breaks allKmers - or at least, it did.
      use_hash.tab <- FALSE
      if (use_hash.tab && !missing(hash.tab)){
         aa_kmers <- hash.tab
      } else {
         aa_kmers <- allKmers(x = gene_prot,k = kmer_size,rename.duplicates = TRUE)
         #aa_kmers <- aa_kmers[,colSums(aa_kmers) < length(gene_prot)]
      }

      kmer_hit_idx <- matrix()
      kmer_hits_pairs_tmp <- expand.grid(1:nrow(orf_kmers),1:nrow(aa_kmers))
      kmer_hits_pairs <- paste(kmer_hits_pairs_tmp$Var1,kmer_hits_pairs_tmp$Var2,sep=',')
      kmer_hits_pairs_mat <- matrix(kmer_hits_pairs,nrow = nrow(orf_kmers))



      kmer_hits_out <- foreach(x = kmer_hits_pairs,.combine=c) %do% {sum(names(orf_kmers)[orf_kmers[as.integer(strsplit(x,',')[[1]][1]),] == 1] %in% names(aa_kmers)[aa_kmers[as.integer(strsplit(x,',')[[1]][2]),] == 1])}
      kmer_hits <- matrix(kmer_hits_out,nrow = nrow(orf_kmers))
      # kmer_hits_in <- matrix(paste(kmer_hits_pairs$Var1,kmer_hits_pairs$Var2,sep=','),nrow = nrow(orf_kmers),byrow = FALSE)
      #kmer_hits_in <- c(paste(1:nrow(orf_kmers)))



      # kmer_hits <- matrix(nrow = nrow(orf_kmers),ncol = nrow(aa_kmers))
      # for (i in 1:nrow(orf_kmers)) { # orf = ith row
      #    for (j in 1:nrow(aa_kmers)) { # protein = jth column
      #       kmer_hits[i,j] <- sum(names(orf_kmers)[orf_kmers[i,] == 1] %in% names(aa_kmers)[aa_kmers[j,] == 1])
      #    }
      # }
      # Row number = which ORF
      # Col number = which protein

      ## Checkpoint: sufficient kmer hits ----

      if (!any(kmer_hits > 0)){
         # message('No kmer hits.')
         output[['log']] <- 'No kmer hits'
         return(output)
      }

      # kmer_cutoff = 5
      if (!any(kmer_hits >= kmer_cutoff)){
         message('Insufficient kmer hits.')
         output[['log']] <- 'n_kmers < kmer_cutoff'
         return(output)
      }

      kmer_match <- as.data.frame(which(kmer_hits >= kmer_cutoff, arr.ind=TRUE)) # The gene and ORFs with sufficient matching kmers
      names(kmer_match) <-  c('orf','prot')
   }

   findORF.out <- findORF(tx_seq = tx_seq,tx_name = tx_name,diamond.df = diamond.df,kmer_size = kmer_size,kmer_cutoff = kmer_cutoff)
   if (is.null(findORF.out)){return(output)}
   kmer_match <- findORF.out$kmer_match; orf_hit <- findORF.out$orf_prediction; gene_prot <- findORF.out$gene_prot; orf_kmers = findORF.out$orf_kmers; aa_kmers = findORF.out$aa_kmers
   ## Process each ORF in kmer_match ----
   #for (orf in unique(kmer_match[,'orf'])){ # for each ORF with kmers matching one of the WT protein sequences...
   #orf_hit <- orf_prediction[orf,]

   #matching_prot_idx <- kmer_match[kmer_match[,'orf'] == orf,'prot']

   matching_gene_prot <- as.data.frame(gene_prot)
   matching_gene_prot$gene <- sapply(X = strsplit(row.names(matching_gene_prot),split = '\\|'), '[',2)
   matching_gene_prot$prot_ref_name <- rownames(matching_gene_prot)
   matching_gene_prot$prot_id <- sapply(strsplit(rownames(matching_gene_prot),'\\|'),FUN = '[',select = 1)
   matching_gene_prot$n_kmers <- kmer_match$n_kmers
   matching_gene_prot$prot_length <- nchar(matching_gene_prot$x)
   matching_gene_prot$kmer_coverage <- (matching_gene_prot$n_kmers + 4)/matching_gene_prot$prot_length
   matching_gene_prot <- matching_gene_prot[order(matching_gene_prot$n_kmers,matching_gene_prot$kmer_coverage,matching_gene_prot$prot_length,decreasing = TRUE),]

   gene_hits <- matching_gene_prot$gene
   #gene_hits <- unique(matching_gene_prot$gene)


   if (length(gene_hits) == 0){
      msg <- 'No matching genes.'
      output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      return(output)

   } else if (length(gene_hits) == 1){
      ### Case: single gene matching ----


      i_orf_kmers <- names(orf_kmers)[which(orf_kmers > 0)]
      i_aa_kmers <- unique(names(aa_kmers)[aa_kmers[row.names(aa_kmers) == row.names(matching_gene_prot)[1],] > 0])
      overlapping_kmers <- i_orf_kmers[i_orf_kmers %in% i_aa_kmers]
      kmer_loc_in_orf <- rep(0,length(overlapping_kmers))
      for (j in 1:length(overlapping_kmers)){
         kmer_loc_in_orf[j] <- regexpr(overlapping_kmers[j],orf_hit[,'seq'])[1]
      }
      kmer_loc_in_wt_prot <- rep(0,length(overlapping_kmers))
      for (j in 1:length(overlapping_kmers)){
         kmer_loc_in_wt_prot[j] <- regexpr(overlapping_kmers[j],matching_gene_prot[1,'x'])[1]
      }

      ####  Check consistency with gene model hypothesis before adding to neopep.list ----
      gene_model_flag <- FALSE
      ol_start <- min(kmer_loc_in_orf)
      ol_end <- max(kmer_loc_in_orf) + kmer_size - 1

      novel_seq = ''

      #if (length(breakpoint.gr) > 0){seqlevelsStyle(breakpoint.gr) <- 'UCSC'}
      #exon_chunk <- strsplit(strsplit(matching_gene_prot$prot_ref_name[1],'_')[[1]][1],'x')[[1]][2]
      #exon_list <- strsplit(exon_chunk,'\\.')[[1]]
      #first_exon <- exon_list[1]
      #last_exon <- exon_list[length(exon_list)]



      matching_prot_id <- matching_gene_prot$prot_id

      prot_cds.gr <- unlist(cdsBy(by = 'tx',
                                  x = edb.ucsc,
                                  columns = 'gene_name',
                                  filter = AnnotationFilterList(ProteinIdFilter(value =  matching_prot_id))))
      #matching_fivePrimeUTR.gr <- fivePrimeUTR.gr[fivePrimeUTR.gr$protein_id %in% matching_prot_id]



      # first_exon.gr <- prot_cds.gr[prot_cds.gr$exon_rank == first_exon]
      # last_exon.gr <- prot_cds.gr[prot_cds.gr$exon_rank == last_exon]

      # suppressWarnings(breakpoint_check <- length(subsetByOverlaps(x = last_exon.gr,ranges = breakpoint.gr)) > 0)
      # if (breakpoint_check){
      #    #message('Breakpoint found within exon ',last_exon.gr$exon_rank,' of ',last_exon.gr$gene_name,'!')
      # }

      # suppressWarnings(fivepUTR_check <- length(subsetByOverlaps(x = first_exon.gr,ranges = matching_fivePrimeUTR.gr)) > 0)
      # if (fivepUTR_check){
      #    #message('TX might hit initiation codon of ',first_exon.gr$gene_name,'!')
      # }

      neo_length_check <- orf_hit$lengthAA - ol_end >= min_length_novel
      ol_start_check <- ol_start <= min_length_novel

      # complete_exon_check <- abs(max(kmer_loc_in_wt_prot) + 4 -
      #                               matching_gene_prot$prot_length[1]) <= 3 # FALSE:
      orf_hit$gene = gene_hits
      orf_hit$gene_5p = gene_hits

      orf_hit$source_prot = matching_gene_prot$prot_ref_name[1]
      orf_hit$kmer_cov = length(overlapping_kmers)/nchar(orf_prediction[orf,'seq'])
      orf_hit$n_kmers = length(overlapping_kmers)
      #orf_hit$fivepUTR <- fivepUTR_check
      #orf_hit$breakpoint <- breakpoint_check
      orf_hit$late_start_seq_flag <- !ol_start_check
      orf_hit$mid_exon_break_flag <- !complete_exon_check

      if (neo_length_check){


         try(novel_seq <- subseq(x = orf_hit$seq,start = ol_end - 9 + min_length_novel,end = nchar(orf_hit$seq)))
         break_pos <- ol_end
         orf_hit$break_pos <- break_pos
         orf_hit$neopep <- novel_seq
         orf_hit$np_aa_length <- nchar(novel_seq)
         orf_hit$name <- paste(tx_name,'_GN_',orf_hit$gene_5p,'_nKs_',orf_hit$n_kmers,sep='')
         message(paste('Neopeptide found matching',gene_hits,'with',length(overlapping_kmers),'kmers on the ',orf_hit$strand,' strand.'))
         rm(break_pos,novel_seq)
      } else {
         orf_hit$break_pos <- 0
         orf_hit$neopep <- ''
         orf_hit$np_aa_length <- 0
         orf_hit$name <- ''
         message('Insufficient novel sequence.')
         msg <- 'No novel seq'
         output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      }


      if(nrow(orf.df) == 0){
         orf.df <- orf_hit

      } else {
         orf.df <- rbind(orf.df,orf_hit)
      } # APPEND orf.df
      # if (neo_length_check) {
      #
      #
      #
      # } else if (!is.numeric(ol_end)){
      #    message('Kmer coverage error')
      #    msg <- "kmer coverage error"
      #    output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      # } else  if (!ol_start_check){
      #    message('Kmer overlap begins too far downstream in ORF, no promoter overlap')
      #    msg <- 'novel lead sequence'
      #    output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      #    # } else if (!complete_exon_check){
      #    # message('Evidence of frameshift-inducing seq error?')
      #    # msg <- 'incomplete exon'
      #    # output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      # } else {
      #    message('Kmer coverage inconsistent with gene model hypothesis')
      #    msg <- 'other truncation gene model violation'
      #    output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      #    #return(output)
      # }

   } else if (length(gene_hits) == 2){




      #message('WARNING: 5pUTR and breakpoint flags not functional in fusion case!')
      ### Case: two genes matching ----
      # if (juliaDebug) {
      #    # run
      #    tx_seq
      #    gene_prot
      #    tx_name
      #    rm_nested
      #    kmer_cutoff
      #    min_length_novel
      # }


      pos_A_gene <- matching_gene_prot[grepl(gene_hits[1],matching_gene_prot$gene),]
      pos_A_gene <- pos_A_gene[order(pos_A_gene$n_kmers,pos_A_gene$kmer_coverage,pos_A_gene$prot_length,decreasing = TRUE),]

      pos_B_gene <- matching_gene_prot[grepl(gene_hits[2],matching_gene_prot$gene),]
      pos_B_gene <- pos_B_gene[order(pos_B_gene$n_kmers,pos_B_gene$kmer_coverage,pos_B_gene$prot_length,decreasing = TRUE),]

      pos_A_gene <- pos_A_gene[1,]
      pos_B_gene <- pos_B_gene[1,]

      fusion_h1_total_kmers <- nchar(orf_prediction[orf,'seq']) - (kmer_size - 1)*2
      kmer_margin <- .05

      #### code to identify position of break -----
      i_orf_kmers <- unique(names(orf_kmers)[orf_kmers[orf,] > 0])
      i_A_aa_kmers <- unique(names(aa_kmers)[aa_kmers[row.names(pos_A_gene),] > 0])
      i_B_aa_kmers <- unique(names(aa_kmers)[aa_kmers[row.names(pos_B_gene)[1],] > 0])
      A_ol_kmers <- i_orf_kmers[i_orf_kmers %in% i_A_aa_kmers]
      B_ol_kmers <- i_orf_kmers[i_orf_kmers %in% i_B_aa_kmers]

      A_kmer_loc_in_prot <- rep(0,length(A_ol_kmers))
      for (j in 1:length(A_ol_kmers)){
         A_kmer_loc_in_prot[j] <- regexpr(A_ol_kmers[j],pos_A_gene$x)[1]
      }

      B_kmer_loc_in_prot <- rep(0,length(B_ol_kmers))
      for (j in 1:length(B_ol_kmers)){
         B_kmer_loc_in_prot[j] <- regexpr(B_ol_kmers[j],pos_B_gene$x)[1]
      }

      A_kmer_loc_in_orf <- rep(0,length(A_ol_kmers))
      for (j in 1:length(A_ol_kmers)){
         A_kmer_loc_in_orf[j] <- regexpr(A_ol_kmers[j],orf_hit[,'seq'])[1]
      }
      B_kmer_loc_in_orf <- rep(0,length(B_ol_kmers))
      for (j in 1:length(B_ol_kmers)){
         B_kmer_loc_in_orf[j] <- regexpr(B_ol_kmers[j],orf_hit[,'seq'])[1]
      }
      A_ol_start <- min(A_kmer_loc_in_orf)
      A_ol_end <- max(A_kmer_loc_in_orf) + kmer_size -1

      #A_break_pos

      B_ol_start <- min(B_kmer_loc_in_orf)
      B_ol_end <- max(B_kmer_loc_in_orf) + kmer_size - 1

      #B_break_pos <- min(B_kmer_loc_in_orf)
      # orf_hit$gene = ''
      # orf_hit$break_pos <- 0
      # orf_hit$neopep <- ''
      # orf_hit$np_aa_length <- 0
      #

      novel_seq <- ''
      fusion_name <- ''
      # default 6.8.21:  fusion_margin <- 10
      ab_orientation <- A_ol_end < B_ol_start
      ba_orientation <- B_ol_end < A_ol_start

      #max(pos_A_gene$prot_length

      complete_exon_check <- ab_orientation | ba_orientation





      if (ab_orientation){
         # complete exon check: gene A completes its exon
         complete_exon_check <- pos_A_gene$prot_length - max(A_kmer_loc_in_prot) <= kmer_size &
            min(B_kmer_loc_in_prot) <= kmer_size

         A_prot_cds.gr <- unlist(cdsBy(x = edb.ucsc,by = 'tx',filter = ProteinIdFilter(pos_A_gene$prot_id)))
         exon_list <- strsplit(strsplit(pos_A_gene$prot_ref_name,'_')[[1]][1],'\\.')[[1]]
         first_exon.gr <- GRanges()
         last_exon.gr <- GRanges()
         try(first_exon.gr <- A_prot_cds.gr[A_prot_cds.gr$exon_rank == strsplit(exon_list[1],':x')[[1]][2] ])
         try(last_exon.gr <- A_prot_cds.gr[A_prot_cds.gr$exon_rank == exon_list[length(exon_list)] ])
         # initiation site check: no leading ORF sequence


         fivepUTR_check <- any(fivePrimeUTR.gr$protein_id == first_exon.gr$protein_id &
                                  fivePrimeUTR.gr$exon_rank == first_exon.gr$exon_rank)

         breakpoint_check <- overlapsAny(last_exon.gr,breakpoint.gr)
         if (length(breakpoint_check) == 0){breakpoint_check = FALSE}
         ol_start_check = min(A_ol_start) <= kmer_size


         break_pos <- c(A_ol_end,B_ol_start)
         fusion_name <- paste(gene_hits,collapse='-')
         try(novel_seq <- subseq(orf_hit$seq,start = break_pos[1] - 9 + min_length_novel,end = break_pos[2] + 9 - min_length_novel))
         orf_hit$gene = fusion_name
         orf_hit$gene_5p = gene_hits[1]
         orf_hit$gene_3p = gene_hits[2]
         orf_hit$source_prot = paste(pos_A_gene$prot_ref_name,pos_B_gene$prot_ref_name,sep='-')
         orf_hit$break_pos <- as.numeric(paste(break_pos,collapse = '.'))
         orf_hit$neopep <- novel_seq
         orf_hit$np_aa_length <- nchar(novel_seq)

         orf_hit$fivepUTR <- fivepUTR_check
         orf_hit$breakpoint <- breakpoint_check
         orf_hit$late_start_seq_flag <- !ol_start_check
         orf_hit$mid_exon_break_flag <- !complete_exon_check


         orf_hit$kmer_cov = (pos_A_gene[1,'n_kmers'] + pos_B_gene[1,'n_kmers'] + 8)/nchar(orf_prediction[orf,'seq'])
         orf_hit$n_kmers <- pos_A_gene[1,'n_kmers'] + pos_B_gene[1,'n_kmers']

         orf_hit$name <- paste(tx_name,'_GN_',orf_hit$gene,'_nKs_',orf_hit$n_kmers,sep='')
         if(nrow(orf.df) == 0){orf.df <- orf_hit} else {orf.df <- rbind(orf.df,orf_hit)} # APPEND orf.df

      } else if (ba_orientation) {


         B_prot_cds.gr <- unlist(cdsBy(x = edb.ucsc,by = 'tx',filter = ProteinIdFilter(pos_B_gene$prot_id)))
         B_exon_list <- strsplit(strsplit(pos_B_gene$prot_ref_name,'_')[[1]][1],'\\.')[[1]]
         first_exon.gr <- B_prot_cds.gr[B_prot_cds.gr$exon_rank == strsplit(B_exon_list[1],':x')[[1]][2] ]
         last_exon.gr <- B_prot_cds.gr[B_prot_cds.gr$exon_rank == B_exon_list[length(B_exon_list)] ]
         # initiation site check: no leading ORF sequence
         if (length(first_exon.gr) == 0 | length(last_exon.gr) == 0){
            message('Warning: some exons from ensembl don"t match gene_prot hit')
         }

         fivepUTR_check <- any(fivePrimeUTR.gr$protein_id == first_exon.gr$protein_id &
                                  fivePrimeUTR.gr$exon_rank == first_exon.gr$exon_rank)

         breakpoint_check <- overlapsAny(last_exon.gr,breakpoint.gr)
         if (length(breakpoint_check) == 0){breakpoint_check = FALSE}
         ol_start_check = B_ol_start <= kmer_size
         complete_exon_check <- pos_B_gene$prot_length - max(B_kmer_loc_in_prot) <= kmer_size &
            min(A_kmer_loc_in_prot) <= kmer_size


         break_pos <- c(B_ol_end,A_ol_start)
         fusion_name <- paste(gene_hits[c(2,1)],collapse='-')
         try(novel_seq <- subseq(orf_hit$seq,start = break_pos[1] - 9 + min_length_novel,end = break_pos[2] + 9 - min_length_novel))
         orf_hit$gene = fusion_name
         orf_hit$gene_5p = gene_hits[2]
         orf_hit$gene_3p = gene_hits[1]
         orf_hit$source_prot = paste(pos_B_gene$prot_ref_name,pos_A_gene$prot_ref_name,sep='-')
         orf_hit$break_pos <- as.numeric(paste(break_pos,collapse = '.'))
         orf_hit$neopep <- novel_seq
         orf_hit$np_aa_length <- nchar(novel_seq)
         orf_hit$fivepUTR <- fivepUTR_check
         orf_hit$breakpoint <- breakpoint_check
         orf_hit$late_start_seq_flag <- ol_start_check
         orf_hit$mid_exon_break_flag <- complete_exon_check


         orf_hit$kmer_cov = (pos_A_gene[1,'n_kmers'] + pos_B_gene[1,'n_kmers'] + 8)/nchar(orf_prediction[orf,'seq'])
         orf_hit$n_kmers <- pos_A_gene[1,'n_kmers'] + pos_B_gene[1,'n_kmers']

         orf_hit$name <- paste(tx_name,'_GN_',orf_hit$gene,'_nKs_',orf_hit$n_kmers,sep='')
         if(length(orf.df) == 0){orf.df <- orf_hit} else {orf.df <- rbind(orf.df,orf_hit)}   # APPEND orf.df
      } else {

         msg <- 'fusion gene model violation'
         message(tx_name,':\t',msg,'\n',paste(orf_hit,collapse='\t'),
                 '\nA_ol:',paste(A_kmer_loc_in_orf,collapse=','),'\tB_ol:',paste(B_kmer_loc_in_orf,collapse=','))
         output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
         eid = strsplit(tx_name,'_')[[1]][2]
         captainsLog(sid = sid,eid = eid,tid = tx_name,msg = msg)

      }








      # new_np.df <- data.frame('orf_no' = orf_prediction[orf,'name'],
      #                         'seq' = orf_prediction[orf,'seq'],
      #                         'length' = nchar(orf_prediction[orf,'seq']),
      #                         'start_pos' = orf_prediction[orf,'startPos'],
      #                         'gene' = fusion_name,
      #                         'kmer_cov' = (pos_A_gene[1,'n_kmers'] + pos_B_gene[1,'n_kmers'] + 8)/nchar(orf_prediction[orf,'seq']),
      #                         'frame' = orf_prediction[orf,'frame'],
      #                         'break' = break_pos)

      # if(length(orf.df) == 0){orf.df <- orf_hit} else {orf.df <- rbind(orf.df,orf_hit)}





   } else {message('Too many gene hits.')

      msg <- 'more than 2 genes involved'
      output[['log']] <- ifelse(output[['log']] == '',msg,paste(output[['log']],msg,sep='...'))
      #why()
   } # close cases: gene_list sizes



   # if (nrow(orf.df) > 1 & length(unique(orf.df$gene)) == 1){
   #    orf.df <- orf.df[orf.df$n_kmers == max(orf.df$n_kmers),]
   # }

   #if (nrow(orf.df)>0){orf.df$weight <- rep(1/nrow(orf.df),nrow(orf.df))}
   # Light formatting ----
   # if (nrow(orf.df) > 0){
   #    orf.df$name <- paste(orf.df$tx_name,'_GN_',orf.df$gene,'_nKs_',orf.df$n_kmers,sep='')
   # }



   # return final output -----
   output[['orf.df']] <- orf.df
   if (output[['log']] == '' & nrow(output[['orf.df']]) > 0) {
      output[['log']] <- ifelse(nrow(orf.df) == 1,'1 neopeptide',paste(nrow(orf.df),'neopeptides'))
      message(output[['log']])

   } else if (output[['log']] == ''){
      output[['log']] <- 'log empty and no ORF output'
   }



   return(output)

}

orfFinderJS <- function(checkSeq,minLength=25) {

   transList <- rep("",6)
   checkList <- unlist(strsplit(checkSeq,""))
   checkList <- checkList[which(checkList%in%c("T","C","G","A"))]
   checkListRev <- checkList[length(checkList):1]
   checkListComp <- as.character(sapply(checkListRev,function(x) c("T","C","G","A")[which(c("A","G","C","T")==x)]))
   transList[1] <- paste(GENETIC_CODE[sapply(seq(1,3*floor(length(checkList)/3),3),function(x) paste(checkList[(x):(x+2)],collapse=""))],collapse="")
   transList[2] <- paste(GENETIC_CODE[sapply(seq(2,3*floor((length(checkList)-1)/3),3),function(x) paste(checkList[(x):(x+2)],collapse=""))],collapse="")
   transList[3] <- paste(GENETIC_CODE[sapply(seq(3,3*floor((length(checkList)-2)/3),3),function(x) paste(checkList[(x):(x+2)],collapse=""))],collapse="")
   transList[4] <- paste(GENETIC_CODE[sapply(seq(1,3*floor(length(checkList)/3),3),function(x) paste(checkListComp[(x):(x+2)],collapse=""))],collapse="")
   transList[5] <- paste(GENETIC_CODE[sapply(seq(2,3*floor((length(checkList)-1)/3),3),function(x) paste(checkListComp[(x):(x+2)],collapse=""))],collapse="")
   transList[6] <- paste(GENETIC_CODE[sapply(seq(3,3*floor((length(checkList)-2)/3),3),function(x) paste(checkListComp[(x):(x+2)],collapse=""))],collapse="")

   startValsTot <- NULL
   endValsTot <- NULL
   revListTot <- NULL
   frameListTot <- NULL
   orfTot <- NULL
   for(i in 1:6) {
      endList <- unique(c(which(unlist(strsplit(transList[i],""))=="*")-1,nchar(transList[1])))
      startList <- unique(c(min(1,endList),which(unlist(strsplit(transList[i],""))=="M")))
      startEnds <- sapply(startList,function(x) endList[min(which(endList>=x))])
      lengthList <- startEnds-startList+1
      whichLength <- which(lengthList>=minLength)
      if(length(whichLength)>0) {
         orfTemp <- sapply(1:length(whichLength),function(x) paste(unlist(strsplit(transList[i],""))[startList[whichLength[x]]:startEnds[whichLength[x]]],collapse=""))
         orfTot <- c(orfTot,orfTemp)
         startVals <-  as.integer(i>3)*(length(checkList)+1)+((-1)^as.integer(i>3))*(3*(startList[whichLength]-1)+(i-1)%%3+1)
         endVals <- as.integer(i>3)*(length(checkList)+1)+((-1)^as.integer(i>3))*(3*(startEnds[whichLength])+(i-1)%%3)
         startValsTot <- c(startValsTot,startVals)
         endValsTot <- c(endValsTot,endVals)
         revListTot <- c(revListTot,rep(c("+","-")[as.integer(i>3)+1],length(orfTemp)))
         frameListTot <- c(frameListTot,rep((i-1)%%3+1,length(orfTemp)))
      }
   }
   if (is.null(startValsTot)) {outPrint <- NULL} else {
      startValsSign <- startValsTot
      startValsSign[which(startValsTot<=3)] <- paste0("<",startValsSign[which(startValsTot<=3)])
      startValsSign[which(startValsTot>=(length(checkList)-2))] <- paste0(">",startValsSign[which(startValsTot>=(length(checkList)-2))])
      endValsSign <- endValsTot
      endValsSign[which(endValsTot<=3)] <- paste0("<",endValsSign[which(endValsTot<=3)])
      endValsSign[which(endValsTot>=(length(checkList)-2))] <- paste0(">",endValsSign[which(endValsTot>=(length(checkList)-2))])

      outPrint <- as.data.frame(cbind(paste0("ORF",1:length(orfTot)),orfTot,revListTot,frameListTot,startValsTot,endValsTot,startValsSign,endValsSign,abs(startValsTot-endValsTot)+1,nchar(orfTot)),stringsAsFactors = FALSE)
      colnames(outPrint) <- c("name","seq","strand","frame","startPos","endPos","startPosSign","endPosSign","lengthNT","lengthAA")
      for (colName in c('frame','startPos','endPos','lengthNT','lengthAA')){
         outPrint[,colName] <- as.numeric(outPrint[,colName])
      }
   }
   return(outPrint)
}


allKmers <- function(x,k,
                     return.adjacency = FALSE,
                     rename.duplicates = FALSE) {

   ux <- list()
   id.list <- character()
   kmer.list <- character()
   n <- length(x)
   x.str <- as.character(gsub('\\.','',x))
   for (i in 1:n) {
      if (is.null(names(x))) {
         id.list <- c(id.list,paste('TX',i,sep = '_'))
      } else {id.list <- c(id.list,names(x[i]))}
      kmers <- substring(x.str[i],1:(nchar(x.str[i]) - k + 1), k:nchar(x.str[i]))
      kmers <- kmers[!grepl('\\*',kmers)]
      kmer.list <- c(kmer.list,kmers)
      ux[[i]] <- kmers
   }

   kmer.list <- unique(kmer.list)
   kmer.freq <- data.frame(matrix(0,nrow = length(id.list),ncol = length(kmer.list)))


   if (rename.duplicates) {
      if (length(unique(id.list)) < length(id.list)) {message('There are duplicated IDs! Patching...')
         id.dups <- which(duplicated(id.list))
         id.list[duplicated(id.list)] <- paste(id.list[duplicated(id.list)],rep('_dup',sum(duplicated(id.list))),sep = '')
      }
   }
   row.names(kmer.freq) <- names(x)
   colnames(kmer.freq) <- kmer.list
   for (i in 1:length(id.list)) {
      kmer.table <- table(ux[[i]])
      for (km in names(kmer.table)) {

         kmer.freq[i,km] <- as.numeric(kmer.table[km])
      }
   }


   kmer.freq <- kmer.freq[,sort(colSums(kmer.freq),decreasing = TRUE,index.return = TRUE)$ix]
   if (return.adjacency) {

      adj.mat <- matrix(data = 0,nrow = length(kmer.list),ncol = length(kmer.list))
      for (j in 2:(nchar(x.str[i]) - k)) {
         prv <- ux[[i]][[j - 1]]
         nxt <- ux[[i]][[j + 1]]
         adj.mat[i,which(kmer.list == prv)] <- 1
         adj.mat[i,which(kmer.list == nxt)] <- 1
         row.names(adj.mat) <- kmer.list
         colnames(adj.mat) <- kmer.list
      }
      return(list("kmer.freq" = kmer.freq,"kmer.adj" = adj.mat))
   } else {return(kmer.freq)
   }

}

count9mers <- function(x){
   n = length(x)
   pep_list <- character()
   for (i in 1:n){
      for (j in 1:(nchar(x[i])-8)){
         pep_list <- c(pep_list,substr(x[i],start = j,stop = j + 8))
      }

   }
   sort(table(pep_list))
   pep_list <- unique(pep_list)
   return(pep_list)
}

findThreePrimeUTR <- function(iORF.df,tx.gr){
   gene_3p <- character(length = nrow(iORF.df))
   forward_gene_match <- genes(x = edb.ucsc,filter = GRangesFilter(value = tx.gr,feature = 'gene'))
   forward_gene_match_extended_window <- genes(x = edb.ucsc,filter = GRangesFilter(value = tx.gr + 1240000,feature = 'gene'))
   reverse_gene_match <- genes(x = edb.ucsc,filter = GRangesFilter(value = invertStrand(tx.gr),feature = 'gene'))
   reverse_gene_match_extended_window <- genes(x = edb.ucsc,filter = GRangesFilter(value = invertStrand(tx.gr + 1240000),feature = 'gene'))
   for (i_orf in 1:nrow(iORF.df)){
      gene_5p <- iORF.df[i_orf,'gene']
      gene_match.gr <- GRanges()
      if (gene_5p %in% forward_gene_match$gene_name){
         gene_match.gr <- forward_gene_match
      } else if (gene_5p %in% reverse_gene_match$gene_name){
         gene_match.gr <- reverse_gene_match
      } else if (gene_5p %in% forward_gene_match_extended_window$gene_name){
         gene_match.gr <- forward_gene_match_extended_window
      } else if (gene_5p %in% reverse_gene_match_extended_window$gene_name){
         gene_match.gr <- reverse_gene_match_extended_window
      }

      if (length(gene_match.gr) == 2){
         # easy case?
         gene_3p[i_orf] <- gene_match.gr$gene_name[gene_match.gr$gene_name != gene_5p]
      } else if (length(gene_match.gr == 1)){
         gene_3p[i_orf] <- 'no_gene'
      } else {
         # difficult case: which gene?
         gene_3p[i_orf] <- 'ambiguous'
      }
   }
   return(gene_3p)
}

jctPosToGRange <- function(event.gr = NULL,extendSearch = 10000){
   extendSearch <- 10000
   genes.gr <- genes(x = edb.ucsc,filter = GRangesFilter(event.gr + extendSearch))
   gene_window_padding <- 1000
   gene_windows.gr <- genes.gr; strand(gene_windows.gr) <- '*'; gene_windows.gr<- range(gene_windows.gr) + gene_window_padding
   nGenes <- length(unique(genes.gr$gene_name))
   min_window_dist <- 50000
   if(start(event.gr[1]) - start(gene_windows.gr[1]) < min_window_dist){start(gene_windows.gr[1]) <- start(event.gr[1]) - min_window_dist}
   if(start(event.gr[2]) - start(gene_windows.gr[2]) < min_window_dist){start(gene_windows.gr[2]) <- start(event.gr[2]) - min_window_dist}

   if(end(gene_windows.gr[1]) - end(event.gr[1]) < min_window_dist){end(gene_windows.gr[1]) <- end(event.gr[1]) + min_window_dist}
   if(end(gene_windows.gr[2]) - end(event.gr[2]) < min_window_dist){end(gene_windows.gr[2]) <- end(event.gr[2]) + min_window_dist}

   return(gene_windows.gr)

}
