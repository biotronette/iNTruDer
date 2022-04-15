# libraries & paths ----
#
# library(IRanges)
# library(Biostrings)
# library(ensembldb)
# library(seqinr)
# library(EnsDb.Hsapiens.v86)
# library(readr)
# library(tictoc)
#
# edb <- EnsDb.Hsapiens.v86
# edb.ucsc <- edb
# ensembldb::seqlevelsStyle(edb.ucsc) <- 'UCSC'
#
# blat_path <- '/research/bsi/tools/biotools/blat/36.5/bin/blat'
# gen_ref_path <- '/research/labs/experpath/vasm/shared/Genome/Human/GRCh38/GCRh38_chromosomesOnly_noMask.fna'
#
# # source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/utilities.R", echo=FALSE)
# source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/intruderPep.R", echo=FALSE)
# source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/intruderAssembly.R", echo=FALSE)
#
#


intruderPipe <- function(reads.filename = NA,
                         reads1 = NA,
                         reads2 = NA,
                         file.prefix = NA,
                         caseId = NA,
                         sampleIdDNA = NA,
                         sampleIdRNA = NA,
                         out.dir = NA,
                         alts.filename = NA,
                         arribaOut.filename = NA,
                         results.outfile = NA,
                         geneExpression.filename = NA,
                         chimeric.bam = NA,
                         aligned.bam = NA,
                         insert_length = NA,
                         reads.threshold = 0,
                         k = seq(11,31,by=4),
                         overwrite_tx = FALSE,
                         overwrite_np = FALSE,
                         overwrite_plots = FALSE,
                         polyA_threshold = 15,
                         min_length_novel = 2,
                         cleanup_tx_dir = FALSE,
                         write.output = TRUE,
                         ref_source = 'none',
                         overwrite_blat = FALSE,
                         cov_cutoff = 'dynamic',
                         flash.threshold = 50,
                         overwrite_ref = FALSE,
                         optimize_velvet = TRUE,
                         tpm_threshold = 30,
                         prot_ref_by = 'prot',
                         numCores = 12){
   # for MCF7 debugging:
   # aligned.bam = "/research/labs/experpath/vasm/shared/NextGen/Projects/Neoantigens/RNAseq/NE48501/GRCh38/maprseq3/output/bams/MCF7_RNA.Aligned.sortedByCoord.out.bam"
   # chimerid.bam = "/research/labs/experpath/vasm/shared/NextGen/Projects/Neoantigens/RNAseq/NE48501/GRCh38/intruder/MCF7_RNA.Chimeric.out.bam"


   # Things that used to be input variables that should maybe just be hard-coded ---
   # debug as of 5/20/21:
   # overwrite_tx = FALSE; overwrite_blat = FALSE; overwrite_plots = FALSE; overwrite_pep = FALSE;
   # write.bed = FALSE; polyA_threshold = 15; k = seq(11,31,by=4); reads.threshold = 0;
   # flash.threshold = 50; cov_cutoff = 'dynamic';
   # Initializing variables ----
   output <- list('seq.df' = data.frame(),
                  'tx.df' = data.frame(),
                  'orf.df' = data.frame(),
                  'np.df' = data.frame(),
                  'np.ss' = AAStringSet(),
                  'blat.gr' = GRanges())
   if (is.na(caseId)){caseId = sampleIdRNA}



   # wrangling the rnaSupportforDNASeqJcts file -----
   if (!dir.exists(dirname(out.dir))){
      dir.create(dirname(out.dir)); dir.create(out.dir)
   } else if (!dir.exists(out.dir)){dir.create(out.dir)}

   seq.df <- loadRNAReads(readsFileName = reads.filename)
   # TODO: verify that loadRNAreads handles case of not merging events correctly
   message('\n*********\nBeginning analysis of ',nrow(seq.df), ' spanning RNA RPs from ',sampleIdRNA,'...')


   write_jct_csv = FALSE
   if(write_jct_csv){
      event_list <- unique(seq.df$DNA_eventId)
      jct_split_dir <- paste(out.dir,'/jct_csv',sep='')
      if (!dir.exists(out.dir)){dir.create(out.dir)}
      if (!dir.exists(jct_split_dir)){dir.create(jct_split_dir)}
      for (j in 1:length(event_list)){
         write.csv(seq.df[seq.df$DNA_eventId == event_list[j],],
                   file = paste(jct_split_dir,'/',file.prefix,'_event_',event_list[j],'.csv',sep=''),
                   quote = FALSE,row.names = FALSE)
      }
   }

   if (nrow(seq.df) == 0){
      message('No RNA reads found for any junctions.')
      return(NULL)
   }

   output[['seq.df']] <- seq.df


   # Transcript wrangling ----

   tx.dir <- file.path(out.dir,'tx')
   if (!dir.exists(tx.dir)){dir.create(tx.dir)}

   tx.df.filename <- file.path(out.dir,paste(sampleIdRNA,'TX_final.csv',sep='_'))
   tx.df_tmp.filename <- file.path(out.dir,'tx',paste(sampleIdRNA,'TX.csv',sep='_'))
   tx.df <- data.frame()
   tx.gr <- GRanges()
   tx.gr_file <- file.path(tx.dir,'tx.alignment.csv')
   # if (is.na(insert_length) & !is.na(folderNum)){insert_length <- getInsertLength(folderNum = folderNum,read_length_source = 'BMD')}

   ## call intruderAssembly ----

   if (overwrite_tx || ! file.exists(tx.df.filename) & ! file.exists(tx.df_tmp.filename) ){



      # debug:
      # reads.threshold = 0;  cov_cutoff = 'dynamic';
      # overwrite_tx <- TRUE ;  cleanup_tx_dir = FALSE;
      # optimize_velvet = TRUE; out.dir = tx.dir; file.prefix = sampleIdRNA

      # for transrate: need aligned, chimera bam files
      #
      cov_cutoff = 'dynamic'; method = 'spades';optimize_velvet = FALSE;
      #out.dir = tx.dir;
      tic()
      assembly.out <- intruderAssembly(reads.filename = reads.filename,
                                       reads1 = reads1,reads2 = reads2,
                                       aligned.bam = aligned.bam,
                                       chimeric.bam = chimeric.bam,
                                       tx.dir = tx.dir,
                                       file.prefix = file.prefix,
                                       reads.threshold = reads.threshold,
                                       insert_length = insert_length,
                                       cov_cutoff = cov_cutoff,
                                       overwrite = overwrite_tx,
                                       overwrite_plots = overwrite_plots,
                                       numCores = numCores,
                                       tpm_threshold = tpm_threshold,
                                       # optimize_velvet = optimize_velvet,
                                       method = method)


      tx.df <- assembly.out$tx.df
      tx.gr <- assembly.out$tx.gr
      tx.toc <- toc()
      message(paste(nrow(tx.df),'TX assembled in',
                    paste(format((tx.toc$toc - tx.toc$tic)/60,digits = 3),'minutes')))
      #diamond.df <- assembly.out$diamond.df
      #quant.df <- assembly.out$quant.df


   } else if (file.exists(tx.df.filename) ) {

      try(tx.df <- read.table(file = tx.df.filename,header=TRUE,as.is=TRUE,sep=','))
      message('Loaded ',nrow(tx.df), ' assembled transcripts from ',tx.df.filename,'.')
   } else if (file.exists(tx.df_tmp.filename) && length(readLines(tx.df_tmp.filename))> 1){
      try(tx.df <- read.table(file = tx.df_tmp.filename,header=TRUE,as.is=TRUE,sep=','))
      message('Loaded ',nrow(tx.df), ' assembled transcripts from ',tx.df_tmp.filename,'.')
      try (tx.gr <- makeGRangesFromDataFrame(read.table(file = tx.gr_file, header = TRUE,sep=','),keep.extra.columns = TRUE))
   }
   if (nrow(tx.df) > 0){
      tx.df <- tx.df[order(tx.df$DNA_eventId,tx.df$tx_num),]
   }

   # Compare to reference proteins w diamond ----
   #tx.df <-
   tx.out.fa <- file.path(tx.dir,'transcripts.fa')
   diamond.df <- runDiamond(tx.dir = tx.dir,tx.out.fa)

   tx.df$inDiamond <- sapply(X = tx.df$tx_num,FUN = function(x) x %in% diamond.df$query)
   if (!'goodTranscript' %in% names(tx.df)){tx.df$goodTranscript <- tx.df$tx_num %in% getFastaHeader(fasta_file = tx.out.fa)}

   # Quantifying read abundance ----

   quant.df <- salmon(tx.dir = tx.dir, salmon_input.fa = tx.out.fa,reads1 = reads1,reads2 = reads2,tpm_threshold = tpm_threshold,overwrite = overwrite_tx)


   # Adding quant, mapping information to tx.df ----

   tx.df$salmon_tpm <- sapply(X = tx.df$tx_num, FUN = function(i) ifelse(i %in% quant.df$Name, quant.df[quant.df$Name == i,'TPM'],0))
   tx.df$numReads <- sapply(X = tx.df$tx_num, FUN = function(i) ifelse(i %in% quant.df$Name, quant.df[quant.df$Name == i,'NumReads'],0))

   quant.df <- read.table(file = file.path(tx.dir,'salmon','quant.sf'),sep='\t',header = TRUE,stringsAsFactors = FALSE)
   quant_pass.df <- quant.df[quant.df$TPM >= tpm_threshold,]


   # Filter TX before proceeding ----
   #tx.df$tpm

   output[['tx.df']] <- tx.df
   output[['tx.gr']] <- tx.gr
   tx_gr.csv <- file.path(tx.dir,'TX_mapping.tsv')
   #blat.df <- read.table(file = tx_gr.csv,sep=' ',header=TRUE,stringsAsFactors = FALSE)
   #blat.df[,c('start','end','width','alignment_length','mismatches','gap_o')]
   blat.gr <- makeGRangesFromDataFrame(df = read.table(file = tx_gr.csv,sep=' ',header=TRUE,stringsAsFactors = FALSE),keep.extra.columns = TRUE)
   #
   #
   #    if (method != 'spades'){
   #       ## call BLAT, blatToGRange ------
   #       blat.gr = GRanges()
   #
   #
   #       if (nrow(tx.df) > 0) {
   #          blat.dir <- file.path(out.dir,'blat')
   #          if (!dir.exists(blat.dir)){dir.create(blat.dir)}
   #          #tx.df.hits <- tx.df[tx.df$geneA_match > 1 | tx.df$geneB_match > 1,]
   #
   #
   #          blat_filename <- file.path(blat.dir,paste(sampleIdRNA,'_junction_tx.rds',sep=''))
   #
   #          if (overwrite_blat | !file.exists(blat_filename)){
   #             tx.ss <- DNAStringSet(tx.df$TX,use.names = TRUE)
   #             names(tx.ss) <- tx.df$tx_num
   #             message('Running blat...')
   #             tic()
   #             blat.df <- blat(stringset = tx.ss,filename = sampleIdRNA,directory = blat.dir)
   #             dt <- toc()
   #             message('BLAT ran in ',format((dt$toc - dt$tic)/60,digits = 3),' minutes; ',format((dt$toc - dt$tic)/nrow(tx.df),digits = 3),' seconds per TX.')
   #             message('Sorting blat results...')
   #             tic()
   #             blat.gr <- suppressWarnings(blatToGRange(blat.df = blat.df,tx.df = tx.df))
   #             blat.gr$event <- as.numeric(unlist(regmatches(blat.gr$query,gregexpr('(?<=J_)\\d+(?=_Loc)',blat.gr$query,perl=TRUE))))
   #             blat.gr <- blat.gr[order(blat.gr$event,blat.gr$query,blat.gr$query_start)]
   #
   #             dt <- toc()
   #             message('BLAT results filtered in ',format((dt$toc - dt$tic)/60,digits = 3),' minutes; ',format(x = (dt$toc - dt$tic)/nrow(tx.df),digits = 3),' seconds per TX.\n')
   #             blat_file <- file(blat_filename)
   #             saveRDS(object = blat.gr,file = blat_file)
   #             close(blat_file)
   #          } else {
   #             blat.gr <- readRDS(blat_filename)
   #
   #          }
   #
   #
   #          output[['blat.gr']] <- blat.gr
   #          tx.df$blat <- tx.df$tx_num %in% blat.gr$query
   #       }
   #    } else if (file.exists(tx_gr.csv)) {
   #
   #       blat.gr <- makeGRangesFromDataFrame(df = read.table(tx_gr.csv,header=TRUE,sep=','),keep.extra.columns = TRUE)
   #       blat.gr$event <- sapply(X = blat.gr$query,FUN = function(i) strsplit(i,'_')[[1]][2])
   #       tx.df$blat <- tx.df$tx_num %in% blat.gr$query
   #    } else {
   #       blat.gr <- GRanges()
   #    }


   ## Writing junction + TX plots ----


   blat_event_list <- unique(blat.gr$event)

   if (overwrite_plots){cl <- makeForkCluster(spec = getOption('cl.cores',numCores))
   parSapply.out <- parSapplyLB(cl = cl, X = blat_event_list,FUN = function(event) runTxViz(event.dir = file.path(tx.dir,paste0('event_',event)),overwrite_plots = overwrite_plots) )
   stopCluster(cl)}


   ## Write protein ref by exon ----


   # debug: overwrite_ref = FALSE; overwrite_plots = FALSE; overwrite_np = TRUE
   # out.dir <- file.path(ppdir,'dev')
   np.dir <- file.path(out.dir,'prot')
   #if (!dir.exists(np.dir)){dir.create(np.dir)}
   if (FALSE){
      np_ref.dir <- file.path(out.dir,'prot','ref')
      if (!dir.exists(np_ref.dir)){dir.create(np_ref.dir,recursive = TRUE)}

      # event_list <- as.integer(unique(blat.gr$event))
      if (prot_ref_by == 'exon'){
         if (length(event_list) > 0){
            message('Checking reference proteome files...')
            make_hash_table <- FALSE
            for (event in event_list){
               e_tx.gr <- blat.gr[grepl(paste('J',event,'',sep = '_'),blat.gr$query)] # removed merged tag
               e_tx.df <- tx.df[tx.df$DNA_eventId == event ,]#& tx.df$kmer_size == 'merged',]
               e_prot_ref_fa_filename <- file.path(np_ref.dir,paste('event',event,'prot_ref.fa',sep='_'))
               e_prot_hash_rds_filename <- file.path(np_ref.dir,paste('event',event,'prot_ref_hash.rds',sep='_'))
               if (length(e_tx.gr) > 0 && (!file.exists(e_prot_ref_fa_filename) | overwrite_ref)){

                  suppressMessages(suppressWarnings(protOut <- getProtExonSeq(e_tx.gr = e_tx.gr)))

                  e_ref_prot <- unique(protOut$gene_prot)
                  if (length(e_ref_prot) ==0  ){
                     # message('No CDS overlaps for event ',event)
                     unlink(e_prot_ref_fa_filename)
                     e_ref_prot <- AAStringSet()
                  } else {
                     message('Found ',length(protOut$gene_prot),' reference protein sequences for event ',event)
                     if (make_hash_table){
                        kmer_size = 5
                        hash.tab <- allKmers(x = e_ref_prot,k = kmer_size,rename.duplicates = TRUE)
                        hash_file <- file(e_prot_hash_rds_filename)
                        saveRDS(object = hash.tab,file = hash_file)
                        close(hash_file)
                     }


                     writeXStringSet(x = e_ref_prot,filepath = e_prot_ref_fa_filename)
                  }

               } else if (length(e_tx.gr) == 0 & file.exists(e_prot_ref_fa_filename)){
                  # message('No CDS overlaps for event ',event)
                  unlink(e_prot_ref_fa_filename)
                  e_ref_prot <- AAStringSet()
               }

            }
         } else {

            message('No events; skipping generation of reference proteome.')
         }
      } else if (prot_ref_by == 'gene'){
         foo <- sapply(X = blat_event_list,FUN = function(event) getProtSeq(diamond.df = diamond.df,blat.gr = blat.gr,event = event,
                                                                            out_fa = file.path(np_ref.dir,paste0('event_',event,'_prot_ref.fa')),overwrite_ref = overwrite_ref ))
         # parSapply.out <- mclapply(cl = cl, X = blat_event_list,FUN = function(event) getProtSeq(diamond.df = diamond.df,event = event,
         #                                                                                            out_fa = file.path(np_ref.dir,paste0('event_',event,'_prot_ref.fa')),overwrite_ref = overwrite_ref))
      }
   }


   # debug polyA_threshold = 15; min_length_novel = 2; merged_tx_only = TRUE;
   # Neopeptide Prediction ------
   neopep_fa_file <- file.path(np.dir,paste(sampleIdRNA,'_neopeptides.fa',sep=''))
   orf_csv_file <- file.path(np.dir,paste(sampleIdRNA,'_orfs.csv',sep=''))
   neopep_csv_file <- file.path(np.dir,paste(sampleIdRNA,'_neopeptides.csv',sep=''))

   tx_log_filename <- file.path(np.dir,paste(sampleIdRNA,'TX_to_NP_log.csv',sep='_'))

   np_input_tx.df <- tx.df[tx.df$tpm >= tpm_threshold &
                              # tx.df$inDiamond &
                              tx.df$goodTranscript,]

   if (nrow(np_input_tx.df) > 0 & (!file.exists(orf_csv_file) | overwrite_np)){
      # loop over tx_seq via findORF ----

      register(MulticoreParam(workers = numCores),default = TRUE)
      tic()
      # cl <- makeForkCluster(nnodes = numCores)

      log.list <- unlist(bplapply(X = 1:nrow(np_input_tx.df),
                                  FUN = function(t) findORF(tx_seq = np_input_tx.df[t,'tx'],
                                                            tx_name = np_input_tx.df[t,'tx_num'],
                                                            diamond.df = diamond.df,
                                                            blat.gr = blat.gr,
                                                            np.dir = np.dir,
                                                            out.file = file.path(np.dir,'tx_out',paste0(np_input_tx.df[t,'tx_num'],'_ORF.csv')))$log))
      # stopCluster(cl)
      toc()
      orf_files <- paste(file.path(np.dir,'tx_out'),list.files(path = file.path(np.dir,'tx_out')),sep='/')
      cat(c("name" ,"seq","strand","frame",
            "startPos" ,"endPos","startPosSign",
            "endPosSign","lengthNT","lengthAA",
            "n_kmers","gp_match","gene","exons","prot",
            "break_pos","leading_aa","novel_length","novel_seq" ),
          sapply(X = orf_files,FUN = function(i) paste0('\n',readLines(i)),USE.NAMES = FALSE),
          file = orf_csv_file,append = FALSE)
      orf.df <- read.table(file = orf_csv_file,fill = TRUE,header = TRUE)

      np_input_tx.df$np_out <- log.list

      orf.df$tpm <- sapply(X = orf.df$name,USE.NAMES = FALSE,FUN = function(tx_name) tx.df[tx.df$tx_num == tx_name,'tpm'])

      np.df <- orf.df[orf.df$novel_length > 0,]
      output[['np.df']] <- np.df


      ## Modify tx.df with np_output -----
      tx.df$np_output <- sapply(X = tx.df$tx_num, USE.NAMES = FALSE,
                                FUN = function(i) ifelse(i %in% np_input_tx.df$tx_num,np_input_tx.df$np_out[np_input_tx.df$tx_num == i],'N/A'))
      output$tx.df <- tx.df; output$orf.df = orf.df

      # tx.df$np_output <- intruderPep.out[['tx_log']]
      # write.table(x =data.frame('np_output'=log.list), file =  tx_log_filename,quote = FALSE,row.names = FALSE,col.names = TRUE)

      # Modify orf.df with seq info ----
      orf.df$breakpointA <- character(nrow(orf.df))
      orf.df$breakpointB <- character(nrow(orf.df))
      for (i_orf in 1:nrow(orf.df)){
         bpA <- unique(seq.df[seq.df$DNA_eventId == getEventFromTxNum(orf.df[i_orf,'name']),'DNA_posA'])
         bpB <- unique(seq.df[seq.df$DNA_eventId == getEventFromTxNum(orf.df[i_orf,'name']),'DNA_posB'])
         chrA <- unique(seq.df[seq.df$DNA_eventId == getEventFromTxNum(orf.df[i_orf,'name']),'chrA'])
         chrB <- unique(seq.df[seq.df$DNA_eventId == getEventFromTxNum(orf.df[i_orf,'name']),'chrB'])
         if (length(bpA) == 1){
            orf.df[i_orf,'breakpointA'] <- paste('chr',chrA,':',bpA,sep='')
         } else {
            orf.df[i_orf,'breakpointA'] <- paste('chr',chrA,':',bpA,sep='',collapse=';')
         }

         if (length(bpA) == 1){
            orf.df[i_orf,'breakpointB'] <- paste('chr',chrB,':',bpB,sep='')
         } else {
            orf.df[i_orf,'breakpointB'] <- paste('chr',chrB,':',bpB,sep='',collapse=';')
         }

      }
      # if (nrow(tx.df) > 0){
      #    tx.df.filename <- ''
      #       write.table(x = tx.df,file = tx.df.filename,quote = FALSE,col.names=TRUE,sep=',')
      # }
      output[['orf.df']] <- orf.df
      np.ss <- AAStringSet(x = np.df$novel_seq); names(np.ss) <- np.df$name
      writeXStringSet(x = np.ss,filepath = neopep_fa_file)
   } else if (nrow(tx.df) > 0 & file.exists(neopep_fa_file ) & file.exists(neopep_csv_file) & !overwrite_np){ # if simply loading output files (overwrite = false)
      #tx_log.df <- read.table(file = tx_log_filename,header=TRUE,as.is=TRUE,sep=',')
      np.ss <- readAAStringSet(filepath = neopep_fa_file)
      orf.df <- read.table(file = orf_csv_file,header=TRUE,as.is=TRUE,sep=',',fill = TRUE)
      np.df <- read.table(file = neopep_csv_file,header = TRUE, as.is = TRUE,sep = ',')
      # tx_log <- read.table(file = tx_log_filename,header = TRUE,row.names = NULL,sep=',')
      # tx.df$np_output <- tx_log[,]
      output[['orf.df']] <- orf.df
      output[['np.df']] <- np.df
      output[['np.ss']] <- np.ss

      if (!'np_output' %in% names(tx.df)){}
   } else if (nrow(tx.df) == 0){

   }
   # Create and populate event.df ----
   event.df <- unique(seq.df[,c('RNA_folderId','DNA_eventId',
                                'RNA_RPs',
                                'DNA_RPs',
                                'chrA','chrB',
                                'DNA_posA','DNA_posB',
                                'geneA','geneB','eventUniqueGenes')])
   event.df <- event.df[order(event.df$DNA_eventId),]
   row.names(event.df) <- NULL
   event.df$breakpointA <- character(nrow(event.df))
   event.df$breakpointB <- character(nrow(event.df))
   event.df$nTX <- numeric(length = nrow(event.df))
   event.df$nTX_filtered <- numeric(length = nrow(event.df))
   event.df$event_tpm <- numeric(length = nrow(event.df))
   # event.df$nRefProt <- numeric(length = nrow(event.df))
   event.df$nDiamond_hits <- numeric(length = nrow(event.df))
   event.df$nNeopep <- numeric(length = nrow(event.df))
   # event.df$ORF_n9mers <- numeric(length = nrow(event.df))
   event.df$NP_n9mers <- numeric(length = nrow(event.df))
   if (anyDuplicated(event.df$DNA_eventId) != 0){event.df <- event.df[!duplicated(event.df$DNA_eventId),]}


   for (i in 1:nrow(event.df)){
      event.df[i,'breakpointA'] <- convertGenomicVectorToString(chr = unique(seq.df$chrA[seq.df$DNA_eventId == event.df[i,'DNA_eventId']]),
                                                                pos = floor(mean(seq.df$DNA_posA[seq.df$DNA_eventId == event.df[i,'DNA_eventId']])))
      event.df[i,'breakpointB'] <- convertGenomicVectorToString(chr = unique(seq.df$chrB[seq.df$DNA_eventId == event.df[i,'DNA_eventId']]),
                                                                pos = floor(mean(seq.df$DNA_posB[seq.df$DNA_eventId == event.df[i,'DNA_eventId']])))
      i_tx.df <- tx.df[tx.df$DNA_eventId == event.df$DNA_eventId[i] ,]#& tx.df$kmer_size == 'merged',]

      i_ref_filename <- file.path(out.dir,'prot/ref',paste('event',event.df$DNA_eventId[i],'prot_ref.fa',sep='_'))


      event.df[i,'nTX'] <- nrow(i_tx.df)
      event.df[i,'nTX_filtered'] <- sum(i_tx.df$tx_num %in% np_input_tx.df$tx_num)
      event.df[i,'event_tpm'] <- sum(quant.df$TPM[sapply(quant.df$Name,getEventFromTxNum,USE.NAMES = FALSE) == event.df[i,'DNA_eventId']])
      # if (file.exists(i_ref_filename)){
      #    event.df[i,'nRefProt'] <- sum(grepl('>',readLines(i_ref_filename))) #ie 'event_122_prot_ref.fa'
      # } else {
      #    event.df[i,'nRefProt'] <- 0
      # }

      if (exists('orf.df')){
         event.df[i,'nDiamond_hits'] <- sum(sapply(X = diamond.df$query,getEventFromTxNum,USE.NAMES = FALSE) == event.df$DNA_eventId[i] )
         # event.df[i,'ORF_n9mers'] <- ifelse(event.df[i,'nORF_hits']>0,length(count9mers(x = orf.df$neopep[orf.df$event == event.df$DNA_eventId[i]])),0)
      }
      if (exists('np.df')){
         event.df[i,'nNeopep'] <- sum(sapply(X = np.df$name,getEventFromTxNum,USE.NAMES = FALSE) == event.df$DNA_eventId[i])
         event.df[i,'NP_n9mers'] <- ifelse(event.df[i,'nNeopep'] > 0,length(count9mers(x = np.df[sapply(np.df$name,getEventFromTxNum,USE.NAMES = FALSE) == event.df$DNA_eventId[i],'novel_seq'])),0)
      }
   }
   event_csv_file <- file.path(out.dir,paste(sampleIdRNA,'_event_summary.csv',sep=''))
   write.table(x = event.df,file = event_csv_file,quote = FALSE,col.names=TRUE,sep=',')

   if (exists('blat.df')){output[['blat.df']] <- blat.df}
   output[['event.df']] <- event.df
   output[['tx.df']] <- tx.df
   if (exists('orf.df')){output[['orf.df']] <- orf.df}
   if (exists('np.df')){output[['np.df']] <- np.df}

   # Create and populate case.df ----

   # rna_sampleId <- getSampleId(folderList = folderNum)
   dna_folderId <- unique(seq.df$DNA_folderId)
   dna_sampleId <- unique(seq.df$DNA_sampleId)
   dna_ppdir <- getPostProcessingDir(dna_folderId)
   jct_file_name <- alts.filename
   jct_file_lines <- readLines(jct_file_name)
   jct.df <- read.delim(file = jct_file_name,skip = grep('eventId',jct_file_lines) - 1,sep=',')

   # if (is.na(maprseq.out.dir) & !is.na(folderNum)){
   #    maprseq.out.dir <- file.path(dirname(getPostProcessingDir(folder = folderNum)),'maprseq3')
   # }
   geneExpression.df <- read.delim(file = geneExpression.filename,sep=',')

   arriba.df <- loadArriba(out.tsv = arribaOut.filename)
   output[['arriba.df']] <- arriba.df


   if (!dir.exists(file.path(out.dir,'reports'))){dir.create(file.path(out.dir,'reports'))}

   # Run caseSummary() ----
   case.df <- caseSummary(caseId = caseId,
                          sampleIdRNA = sampleIdRNA,
                          sampleIdDNA = sampleIdDNA,
                          tx.df = tx.df,
                          jct.df = jct.df,
                          orf.df = orf.df,
                          np.df = np.df,
                          blat.gr = blat.gr,
                          geneExpression.df = geneExpression.df,
                          arriba.df = arriba.df,
                          event.df = event.df)
   write.table(case.df,file = results.outfile,quote = FALSE,sep='\t',row.names = FALSE,col.names = TRUE)
   # case_tmp.df <- read.table(case_summary_file,header = TRUE,sep='\t',quote = '')
   output[['case.df']] <- case.df
   output[['arriba.df']] <- arriba.df
   output[['jct.df']] <- jct.df
   output[['geneExpression.df']] <- geneExpression.df

   if (cleanup_files){

   }

   saveRDS(output,file = file.path(out.dir,paste(sampleIdRNA,'.rds',sep='')))
   return(output)
}

blat <- function(stringset,filename,directory,
                 ref_stringset_file = gen_ref_path,
                 allow_small_segments = FALSE){
   # debug
   # blat.dir <- file.path(out.dir,'blat');    filename = sampleIdRNA;
   # tx.ss <- DNAStringSet(tx.df$TX,use.names = TRUE) ;   tx.df.names <- tx.df$tx_num;   names(tx.ss) <- tx.df.names
   # stringset = tx.ss; filename = sampleIdRNA;directory = blat.dir


   # if (allow_small_segments) {
   #    ref_input <- paste(directory,'/',filename,'.fa',sep = '')
   #    writeXStringSet(ref_stringset,ref_input,append=FALSE,format='fasta')
   #    #faTo2bit.path <- paste('/',c(strsplit(blat_path,'/')[[1]][2:5],'faToTwoBit'),collapse='',sep='')
   #    reference_path <- paste(directory,'/',filename,'.fa',sep = '')
   #    # faTo2bit.cmd <- paste(faTo2bit.path,ref_input,reference_path)
   #    # system(faTo2bit.cmd)
   # }
   reference_path <- ref_stringset_file

   infile <- paste(directory,'/',filename,'.in',sep = '')
   outfile <- paste(directory,'/',filename,'.psl',sep = '')
   writeXStringSet(stringset, infile, append = FALSE, format = "fasta")
   if (allow_small_segments){segment_options <- '-tileSize=6 -minMatch=1 -minScore=1'} else {segment_options <- ''}
   blat.cmd <- paste(blat_path,reference_path,infile,'-out=blast8 -minIdentity=90 -maxGap=0',segment_options,outfile) # TODO: fix blat_path variable
   system(blat.cmd)
   psl_names <- c('match','mis_match','rep_match','Ns','Q_gap_count','Q_gap_bases','T_gap_count','T_gap_bases','strand',
                  'Q_name','Q_size','Q_start','Q_end','T_name','T_size','T_start','T_end','block_count','blockSizes',
                  'qStarts','Tstarts')
   blast8_names <- c('query','chr','percent_id','alignment_length','mismatches','gap_o',
                     'query_start','query_end','subject_start','subject_end','E_value','bit_score')
   blat_results <- read.table(outfile,header = FALSE,sep = '\t',col.names = blast8_names,stringsAsFactors = FALSE)
   if (any(grepl('NC',blat_results$chr))){
      blat_results <- blat_results[grepl('NC',blat_results$chr),]
      blat_results$chr <- as.numeric(regmatches(blat_results$chr,regexpr('(?<=0{4})\\d{2}',blat_results$chr,perl = TRUE)))
   } else {
      names(blat_results) <- c('query','chr','percent_id','alignment_length','mismatches','gap_o',
                               'query_start','query_end','subject_start','subject_end','E_value','bit_score')
   }

   blat_results$strand <- ifelse(blat_results$subject_end - blat_results$subject_start > 0,'+','-')
   # for (j in unique(ss_names)) {
   #   j.df <- blat_results[blat_results$query == j,]
   #   j.df <- j.df[j.df$percent_id > percent_id & j.df$bit_score > bit_score,]
   #   if (nrow(j.df) > 0){
   #     if (sum(j.df$alignment_length) - unique(j.df$))
   #   }
   # }
   #blat_sort <- order(blat_results$Q_name,blat_results$Q_start,blat_results$Q_end)
   return(blat_results)
}



blatToGRange <- function(blat.df,tx.df = NA,
                         percent_id = 95,
                         bit_score = 100,
                         processAsReads = FALSE,
                         verbose = FALSE,
                         map_cov_cutoff = .99){
   # debug: blat.df = blat.df;tx.df = tx.df
   # percent_id = 95;   bit_score = 100;   processAsReads = FALSE;   verbose = FALSE;   map_cov_cutoff = .7
   narrowBy = 6
   blat.gr <- GRanges()

   event_list <- unique(tx.df$DNA_eventId)

   for (i_event in event_list) {
      #message(j)


      tx_list <- unique(as.character(tx.df$tx_num[tx.df$DNA_eventId == i_event]))
      for (i_tx in tx_list) {
         #message(t)
         tx_idx <- which(tx.df$DNA_eventId == i_event & tx.df$tx_num == i_tx)
         tx_length <- unique(tx.df[tx_idx,'tx_length'])
         tx_conf <- unique(tx.df[tx_idx,'tx_conf'])

         if (processAsReads) {
            query_name <- paste(i_event,i_tx,sep = '_')
         } else if (!grepl('event',i_tx)){
            query_name <- i_tx #paste('event',formatC(j,width = 3,flag = '0'),t,sep = '_')
         } else {query_name <- i_tx}

         if (any(grepl('chr',blat.df$chr)) & !grepl('chr',unique(tx.df$chrA[tx_idx]))){
            chrA <- paste('chr',unique(tx.df$chrA[tx_idx]),sep='')
            chrB <- paste('chr',unique(tx.df$chrB[tx_idx]),sep='')
         } else {
            chrA = unique(tx.df$chrA[tx_idx])
            chrB = unique(tx.df$chrB[tx_idx])
         }

         unfiltered_blat_hits <- unique(blat.df[blat.df$query == query_name &
                                                   blat.df$chr %in% c(chrA,chrB),])

         blat_hits <- unique(blat.df[blat.df$query == query_name &
                                        blat.df$chr %in% c(chrA,chrB) &
                                        blat.df$percent_id >= percent_id &
                                        blat.df$alignment_length > 2*narrowBy &
                                        blat.df$mismatches <= blat.df$alignment_length*0.1,])
         # blat.gr <- GRanges()
         if (nrow(blat_hits) > 0){
            blat_hits$map_cov <- blat_hits$alignment_length/tx_length

            # call blatFilter ----
            #debug: blat_hits = blat_hits;map_cov_cutoff = map_cov_cutoff;tx_length = tx_length
            blatFilter_out.gr <- blatFilter(blat_hits = blat_hits,
                                            map_cov_cutoff = map_cov_cutoff,
                                            tx_length = tx_length,
                                            narrowBy = narrowBy)
         } else {
            blatFilter_out.gr <- GRangesList()
         }

         blat.gr <- c(blat.gr,blatFilter_out.gr)




      }
   }
   blat.gr$chr = as.character(seqnames(blat.gr))
   miss_list <- tx.df$tx_num[!tx.df$tx_num %in% unique(mcols(blat.gr)[,'query'])]
   if (length(miss_list) > 0){message(paste('Couldn\'t find hits covering ',map_cov_cutoff*100,'% of tx for ',
                                            length(miss_list),' queries/TX, ',format(length(miss_list)/nrow(tx.df)*100,digits = 3),'% of total.',sep = ''))}
   return(blat.gr)
}


blatFilter <- function(blat_hits,tx_length,narrowBy = 6,map_cov_cutoff = .99){
   if (missing(tx_length)) {tx_length <- max(blat_hits$query_end)}
   filtered_output <- data.frame()
   blat_hits_input <- blat_hits
   blat_hits$old_chr <- blat_hits$chr
   blat_hits$chr <- 42
   # narrowBy <- 6


   nHits = 0
   blat.gr <- GRanges()
   if (nrow(blat_hits) == 1) {
      # if (blat_hits[1,'map_cov'] < map_cov_cutoff) {
      #    blat_hits <- data.frame()
      # } else {
      filtered_output <- blat_hits
      nHits = 1
      # }


   } else if (length(blat_hits) > 1) {

      blat_sort <- order(blat_hits$bit_score,decreasing = TRUE)
      blat_hits <- blat_hits[blat_sort,]


      if (!any(countOverlaps(makeGRangesFromDataFrame(blat_hits,start.field = 'query_start',end.field = 'query_end')) > 1) &
          abs(sum(blat_hits$alignment_length) - tx_length) / tx_length < (1 - map_cov_cutoff)) {
         filtered_output <- blat_hits
         nHits <- nrow(blat_hits)
      } else  { # if blat_hits doesn't meet map_cov_cutoff threshold

         nHits <- 0

         top_hits <- makeGRangesFromDataFrame(blat_hits[1,],
                                              start.field='query_start',
                                              end.field='query_end',
                                              keep.extra.columns =TRUE)
         other_hits <- makeGRangesFromDataFrame(blat_hits[2:nrow(blat_hits),],
                                                start.field='query_start',
                                                end.field='query_end',
                                                keep.extra.columns = TRUE)
         #narrowBy <-min( floor(width(top_hits)*0.05))
         remainders <- subsetByOverlaps(x = other_hits,ranges = top_hits-narrowBy,invert = TRUE,ignore.strand=TRUE)
         #remainders <- subsetByOverlaps(x = other_hits,ranges = narrow(top_hits,start = narrowBy + 1,end = width(top_hits) - narrowBy),invert = TRUE,ignore.strand=TRUE)
         #remainders <- remainders[width(remainders) <= tx_length - sum(width(top_hits)) + 4]
         #remainders <- remainders
         cum_map_cov <- abs(sum(width(top_hits))) / tx_length
         if (length(remainders) > 0){
            while (length(remainders) > 0 &
                   cum_map_cov < 1 &
                   sum(width(top_hits)) < tx_length
            ){


               top_hits <- c(top_hits,remainders[1])
               new_cum_map_cov <- abs(sum(width(top_hits))) / tx_length
               remainders <- subsetByOverlaps(x = remainders,ranges = top_hits-narrowBy,invert=TRUE,ignore.strand=TRUE)

               if (new_cum_map_cov >= map_cov_cutoff | length(remainders) == 0) {
                  nHits <- length(top_hits)
                  filtered_output <- as.data.frame(top_hits)
                  filtered_output$query_start <- blat_hits$start
                  filtered_output$query_end <- blat_hits$end
                  # for (i in 1:length(filtered_output)){
                  #   filtered_output[i,'old_chr'] <- blat_hits_input[blat_hits$query == filtered_output[i,'query'],'chr']
                  # }
                  # filtered_output$old_chr <- blat_hits$old_chr



                  # break
               }

               cum_map_cov <- new_cum_map_cov
               #remainders <- remainders[width(remainders) <= tx_length - sum(width(top_hits)) + 3 + 1]

            }
         } else if (length(remainders) == 0 ){
            nHits <- length(top_hits)
            filtered_output <- as.data.frame(top_hits)
            filtered_output$query_start <- blat_hits$start
            filtered_output$query_end <- blat_hits$end
         }


         if (nHits == 0) {
            #miss_list <- c(miss_list,query_name)

            filtered_output <- data.frame()
         } #else {blat_hits <- blat_hits[1:nHits,]}
      } # if first element or all blat_hits meets map_cov_cutoff threshold after simple QC, no further modification needed
      if (nHits == 0){
         filtered_output <- data.frame()
      }


   }

   if(nHits  == 0) {
      return(NULL)
   } else {
      if (any(filtered_output$strand == '-')) {
         filtered_output[filtered_output$strand == '-',c('subject_start','subject_end')] <- filtered_output[filtered_output$strand == '-',c('subject_end','subject_start')]
      }
      if ('chr' %in% names(filtered_output)){
         filtered_output$seqnames <- filtered_output$old_chr
         #filtered_output$chr <- NULL
      }
      #filtered_output$seqnames <- filtered_output$old_chr
      #filtered_output$chr <- NULL
      # TODO: check definition of query_start and query_end in output
      if (!'query_start' %in% names(filtered_output) || any(is.na(filtered_output[,'query_start']))){
         filtered_output$query_start <- filtered_output$start
         filtered_output$query_end <- filtered_output$end
      }


      filtered_output$start <- filtered_output$subject_start
      filtered_output$end <- filtered_output$subject_end
      filtered_output <- filtered_output[order(filtered_output$query_start),]
      #filtered_output$query_start <-
      if (any(filtered_output$start > filtered_output$end)){
         filtered_output[filtered_output$start > filtered_output$end,c('start','end')] <- filtered_output[filtered_output$start > filtered_output$end,c('end','start')]
      }
      filtered_output$seqnames <- NULL
      row.names(filtered_output) <- NULL
      blat.gr <- makeGRangesFromDataFrame(filtered_output,
                                          seqnames.field = 'old_chr',
                                          start.field = 'start',
                                          end.field = 'end',
                                          strand.field = 'strand',
                                          keep.extra.columns = TRUE)
      blat.gr$query_map_cov <- sum(blat.gr$map_cov)
      blat.gr <- blat.gr[order(blat.gr$query,blat.gr$query_start),]

   }


   #else {blat.gr <- GRanges()}
   return(blat.gr)

}

getProtExonSeq <- function(e_tx.gr,exon_margin = 2){
   # debug: tx.gr <- e_tx.gr
   if (length(e_tx.gr) == 0){
      return(list('gene_prot' = AAStringSet()))
   }

   cds.full.gr <- unlist(cdsBy(x = edb.ucsc,
                               by = 'tx',
                               filter = GRangesFilter(value = c(e_tx.gr,invertStrand(e_tx.gr)),type = 'any'),
                               columns = c('gene_name','protein_id','exon_idx','tx_id','tx_biotype')))

   cds.gr <- GRanges()
   if (length(cds.full.gr) > 0){
      cds.gr = subsetByOverlaps(x = cds.full.gr,ranges = e_tx.gr + 2000,ignore.strand = TRUE)
      cds.gr <- cds.gr[!grepl('LRG_',names(cds.gr))]
   }

   plot_cds <- FALSE
   if (plot_cds & length(cds.gr) > 0){
      tx_id_nums <- Rle(names(cds.gr)); runValue(tx_id_nums) <- 1:length(runValue(tx_id_nums))
      cds.gr$tx_num <- as.vector(tx_id_nums)
      cds.gr$tx_id <- names(cds.gr)

      for (i_chr in unique(seqnames(cds.gr))){
         x.gr <- cds.gr[seqnames(cds.gr) == i_chr]
         ggplot(data = x.gr) + geom_arrowrect(aes(color=tx_biotype,fill=gene_name,y = tx_num),
                                              stat = 'identity') + scale_y_discrete(name = 'Ensembl tx_id',breaks = names(x.gr))
      }
   }

   proteins <- data.frame()

   gene_prot <- AAStringSet()
   cds.pa <- GRanges()
   if (length(cds.gr) == 0) { return(list('gene_prot' = gene_prot,'cds.pa' = cds.pa,'proteins' = proteins)) }



   unique_cds.gr <- GRangesList()

   use_flank <- FALSE

   if (use_flank){
      cds_plus_flank.gr <- GRanges()
      for (i_cds in 1:length(cds.gr)){
         ex_net_radius <- 2
         cds_plus_flank.gr <- c(flanking_cds.gr,cds.full.gr[cds.full.gr$protein_id ==
                                                               mcols(cds.gr)[i_cds,'protein_id'] &
                                                               abs(cds.full.gr$exon_rank - mcols(cds.gr)[i_cds,'exon_rank']) <
                                                               ex_net_radius,])
      }
      try(unique_cds.gr <- unique(cds_plus_flank.gr))
   } else {
      try(unique_cds.gr <- unique(cds.gr))
   }
   #full_protAnnotation <-  suppressWarnings(unlist(genomeToProtein(x = exons.gr,edb.ucsc)))

   gene_name <- unique(cds.gr$gene_name)
   #exon_tx_ol.gr <- subsetByOverlaps(ranges = tx.gr,x = exons.gr,ignore.strand  = TRUE)
   #unique_prot <- unique(exon_tx_ol.gr$protein_id)
   #unique_cds.gr <- unique(cds.gr)
   if (length(unique_cds.gr) > 0){
      cds_unfiltered.pa <- suppressWarnings(unlist(genomeToProtein(x = unique_cds.gr,edb.ucsc))   )
   } else {
      cds_unfiltered.pa <- suppressWarnings(unlist(genomeToProtein(x = cds.gr,edb.ucsc)))
   }

   cds_filtered.pa <- cds_unfiltered.pa[mcols(cds_unfiltered.pa)[,'tx_id'] %in% names(cds.gr)]
   if (!is.null(names(cds_filtered.pa))){
      mcols(cds_filtered.pa)[,'protein_id'] <- names(cds_filtered.pa)
   }

   for (i_cds in which(start(cds_filtered.pa) == -1)){ #
      cds.pid <- unique(mcols(cds.gr[names(cds.gr) == mcols(cds_filtered.pa)[i_cds,'tx_id']])[,'protein_id'])
      if (length(cds.pid) == 1){
         cds_prot_matches.prot <- proteins(edb.ucsc,filter = ProteinIdFilter(cds.pid))
         cds_prot_matches.prot$length <- nchar(cds_prot_matches.prot$protein_sequence)

         upstream_exon.gr <- cds.full.gr[cds.full.gr$protein_id == cds.pid &
                                            cds.full.gr$exon_rank == (mcols(cds_filtered.pa)[i_cds,'exon_rank'] - 1)]

         prot_cds.gr <- cds.full.gr[cds.full.gr$protein_id == cds.pid]

         if (length(upstream_exon.gr) > 0){

            suppressWarnings(upstream_exon.pa <- unlist(genomeToProtein(x = upstream_exon.gr,db = edb.ucsc)))
            if (!is.null(names(upstream_exon.pa))){
               mcols(upstream_exon.pa)[,'protein_id'] <- names(upstream_exon.pa)
               upstream_exon.pa <- upstream_exon.pa[mcols(upstream_exon.pa)[,'protein_id'] == cds.pid]


               if (length(upstream_exon.pa) > 0 && end(upstream_exon.pa) < cds_prot_matches.prot$length){
                  pid.ir <- IRanges(start = end(upstream_exon.pa) + 1,end= cds_prot_matches.prot$length)
                  names(pid.ir) <- cds_prot_matches.prot$protein_id
                  #rev_cds_from_prot.pa <- GRanges()
                  try(rev_cds_from_prot.pa <- unlist(proteinToGenome(x = pid.ir,db = edb.ucsc))[[cds.pid]])
                  rev_cds_check <- length(rev_cds_from_prot.pa) > 0 &&
                     !any(mcols(cds_filtered.pa)[,'exon_id'] == mcols(rev_cds_from_prot.pa)[,'exon_id'] & start(cds_filtered.pa) != -1)

                  if (rev_cds_check){
                     if (mcols(rev_cds_from_prot.pa)[,'tx_id'] == mcols(cds_filtered.pa)[i_cds,'tx_id'] &
                         mcols(rev_cds_from_prot.pa)[,'exon_rank'] == mcols(cds_filtered.pa)[i_cds,'exon_rank']){
                        end(cds_filtered.pa[i_cds]) <- mcols(rev_cds_from_prot.pa)[,'protein_end']
                        start(cds_filtered.pa[i_cds]) <- mcols(rev_cds_from_prot.pa)[,'protein_start']
                        #mcols(cds_filtered.pa)[i,'ref_patch_flag'] <- TRUE
                        mcols(cds_filtered.pa)[i_cds,'protein_id'] <- cds.pid
                     }
                  }


               }
            }
         }
      }



   }



   if (length(cds_filtered.pa) > 0 & 'protein_id' %in% names(mcols(cds_filtered.pa))){
      if (any(names(cds_filtered.pa) == '' & mcols(cds_filtered.pa)[,'protein_id'] != '')){
         names(cds_filtered.pa) <-  mcols(cds_filtered.pa)[,'protein_id']
      }
      cds.pa <- cds_filtered.pa[start(cds_filtered.pa) != -1]
   }

   if (length(cds.pa) > 0){

      for (i_cds in 1:length(cds.pa)){
         i_exon_id <- mcols(cds.pa)[['exon_id']][i_cds]
         name_match <- unique(cds.gr$gene_name[ cds.gr$exon_id == i_exon_id])
         # name_match <- unique(exons.gr$gene_name[which(exons.gr$exon_id == mcols(cds.pa)[['exon_id']][i])])
         #message(name_match)
         if (length(name_match) == 0){
            #mcols(cds.pa)[['gene_name']][i] <- ''
         } else {
            mcols(cds.pa)[['gene_name']][i_cds] <- name_match
         }

      }

      prot_ids <- mcols(cds.pa)[,'protein_id']
      names(cds.pa) <- mcols(cds.pa)[,'protein_id']
      proteins <- as.data.frame(proteins(edb.ucsc,
                                         columns = c("gene_name","protein_id","protein_sequence"),
                                         filter = ProteinIdFilter(prot_ids)),stringsAsFactors = FALSE)
      # if (length(unique(names(cds.pa)))> 1){
      #    message(print(tx_name,'has multiple genes matching its protein annotation! wwwhhhoooooaaaaaa'))
      # } else {
      #    mcols(cds.pa)[['prot_seq']] <- subseq(x = rep(proteins$protein_sequence[proteins$protein_id == unique(names(cds.pa))],length(cds.pa)),
      #                                                  start = start(cds.pa),end=end(cds.pa))
      # }
      #prot_ids <- unique(names(cds.pa))

      # gene_prot <- AAStringSet()
      if (length(prot_ids) == 1){
         gene_prot <- protToExonSeqs(cds.pa,proteins,prot_ids)

      } else if (length(prot_ids) > 1){
         #gene_prot <- AAStringSet()
         for (prot in prot_ids){
            # call protToExonSeqs() ----
            # cds.pa <- cds.pa; prot = prot_ids[1]

            prot.ss <- protToExonSeqs(cds.pa,proteins,prot)
            #prot.ss <- unique(prot.ss)
            if (length(gene_prot) == 0){
               gene_prot <- prot.ss
            } else if (length(prot.ss) > 0){gene_prot <- c(gene_prot,prot.ss)}
         }
         gene_prot <- unique(gene_prot)

         # TODO: handle fusion cases of protein annotation
      } else if (length(prot_ids) == 0){
         message('No matching proteins found.')
         #gene_prot = getUniprotSeq(gene_list,pe_filtered = TRUE)
      }
   } else if (length(cds.gr) == 0) {
      message('No CDS overlaps found.')
      #gene_prot <- getUniprotSeq(gene_list,pe_filtered = TRUE)
   }

   if (length(gene_prot) > 0){
      gene_prot <- unique(gene_prot)
      gene_prot[width(gene_prot) <= 9] <- NULL
   }

   return(list('gene_prot' = gene_prot,'cds.pa' = cds.pa,'proteins' = proteins))#,'exons.gr' = exons.gr))
   #return(gene_prot)
}


get5pUTRbyOverlap <- function(tx.gr){ #output:  fivePrimeUTR.gr
   fputr.gr <- GRanges()
   fputr.grl <- fiveUTRsByTranscript(x = edb.ucsc,
                                     columns = c('gene_name','protein_id'),
                                     filter = GRangesFilter(c(tx.gr,invertStrand(tx.gr)),type = 'any'))
   fivePrimeUTR.gr <- GRanges()
   if (length(fputr.grl) > 0){
      for (i_np in 1:length(fputr.grl)){

         if (unique(strand(fputr.grl[[i_np]]))=='+'){
            last.gr <- fputr.grl[[i_np]][end(fputr.grl[[i_np]]) == max(end(fputr.grl[[i_np]]))]
            utr_end.gr <- GRanges(seqnames = seqnames(last.gr),ranges = IRanges(start = end(last.gr),end = end(last.gr) + 1),strand = strand(last.gr))
            mcols(utr_end.gr) <- mcols(last.gr)
         } else {
            last.gr <- fputr.grl[[i_np]][start(fputr.grl[[i_np]]) == min(start(fputr.grl[[i_np]]))]
            utr_end.gr <- GRanges(seqnames = seqnames(last.gr),ranges = IRanges(start = start(last.gr)-1,end = start(last.gr) ),strand = strand(last.gr))
            mcols(utr_end.gr) <- mcols(last.gr)
         }
         if (length(fivePrimeUTR.gr) == 0){fivePrimeUTR.gr <- utr_end.gr} else {fivePrimeUTR.gr <- c(fivePrimeUTR.gr,utr_end.gr)}

      }
   }

   return(fivePrimeUTR.gr)
}

getPromotersByOverlap <- function(tx.gr){
   promoters.gr <- GRanges()
   try(promoters.gr <- promoters(x = edb.ucsc,
                                 upstream = 1,
                                 downstream = 1,
                                 use.names = TRUE,
                                 columns = c('gene_name','protein_id'),
                                 filter = GRangesFilter(c(tx.gr,invertStrand(tx.gr)),type = 'any')))
   if (length(promoters.gr) > 0){promoters.gr <- promoters.gr[!is.na(promoters.gr$protein_id)]}
   return(promoters.gr)
}

getBreakpointByOverlap <- function(tx.gr){
   breakpoint.gr <- GRanges()
   if (length(tx.gr) == 1){
      message('TX ',unique(tx.gr$query),' only maps to one region; skipping breakpoint prediction.')
      return(breakpoint.gr)
   }
   # foo <- disjointExons(x = edb.ucsc,
   #                      filter = GRangesFilter(value = c(tx.gr,invertStrand(tx.gr)),type = 'any'),
   #                      includeTranscripts = TRUE,
   #                      aggregateGenes = FALSE)
   cds.tmp.list <- cdsBy(x = edb.ucsc,
                         by = 'tx',
                         filter = GRangesFilter(value = c(tx.gr,invertStrand(tx.gr)),type = 'any'),
                         columns = c('gene_name','gene_id','tx_id','protein_id'))
   if (length(cds.tmp.list) > 0){
      cds_mask_list <- logical(length = length(cds.tmp.list))
      for (i_np in 1:length(cds.tmp.list)){
         cds_mask_list[i_np] <- cds.tmp.list[i_np] %within% cds.tmp.list[-i_np]
      }
      cds.full.gr <- unlist(cds.tmp.list); cds.full.gr <- cds.full.gr[!grepl('LRG_',cds.full.gr$protein_id) ]
      cds.ss.gr <- unlist(cds.tmp.list[!cds_mask_list]); cds.ss.gr <- cds.ss.gr[!grepl('LRG_',cds.ss.gr$protein_id) ]


      cds.tmp <- cds.ss.gr
   } else {cds.tmp <- GRanges()}




   # dj_exons.gr <- disjointExons(x = edb.ucsc,
   #                              filter = AnnotationFilterList(GRangesFilter(value = cds.tmp,type = 'any'),
   #                                                            ProteinIdFilter(unique(cds.tmp$protein_id))),
   #                              includeTranscripts = FALSE,
   #                              aggregateGenes = FALSE)
   # dj_exons.gr$gene_name <- character(length(dj_exons.gr))
   # for (i in 1:length(dj_exons.gr)){
   #    dj_exons.gr$gene_name[i] <- unique(cds.tmp[cds.tmp$gene_id == unlist(dj_exons.gr$gene_id[i])]$gene_name)
   #    # dj_exons.gr$
   # }

   if (length(cds.tmp) > 0){





      cds_ol_tx.gr <- subsetByOverlaps(x = cds.tmp,ranges = tx.gr,ignore.strand = TRUE,type = 'any')
      tx_ol_cds.gr <- subsetByOverlaps(ranges = cds.tmp,x = tx.gr,ignore.strand = TRUE,type = 'any')

      # tx_within_cds.gr <- subsetByOverlaps(ranges = cds.tmp,x = tx.gr,ignore.strand = TRUE,type = 'any')
      tx_ol_cds.gr$mid_exon_break <- logical(length = length(tx_ol_cds.gr))
      tx_ol_cds.gr$break_gene <- character(length = length(tx_ol_cds.gr))
      #tx_ol_cds.gr$break_exon <- character(length = length(tx_ol_cds.gr))

      if (length(tx_ol_cds.gr) > 0){
         for (i_np in 1:length(tx_ol_cds.gr)){



            matching_cds.gr <- subsetByOverlaps(x = cds.tmp,
                                                ranges = tx_ol_cds.gr[i_np],
                                                ignore.strand = TRUE)
            if (length(unique(matching_cds.gr)) > 1 & length(unique(matching_cds.gr$gene_name)) == 1){



               matching_cds_superset_range <- range(matching_cds.gr)
               matching_cds.gr <- matching_cds.gr[ranges(matching_cds.gr) == ranges(matching_cds_superset_range)]
            }

            fw_strand_ds_exon_break_flag <- as.character(strand(matching_cds.gr)) == '+' &
               abs(end(matching_cds.gr) - end(tx_ol_cds.gr[i_np])) > 2

            rev_strand_ds_exon_break_flag <- as.character(strand(matching_cds.gr)) == '-' &
               abs(start(tx_ol_cds.gr[i_np]) - start(matching_cds.gr)) > 2
            # tx_ol_cds.gr$exon_hit[i] <- paste()

            ol_flag <- any(fw_strand_ds_exon_break_flag,rev_strand_ds_exon_break_flag)
            if (!ol_flag %in% c(TRUE,FALSE)){
               ol_flag <- FALSE
               why()
            }
            tx_ol_cds.gr$mid_exon_break[i_np] <- ol_flag
            if (ol_flag){
               tx_ol_cds.gr$break_gene[i_np] <- paste(unique(matching_cds.gr$gene_name[fw_strand_ds_exon_break_flag | rev_strand_ds_exon_break_flag]),collapse='; ')
            }

         }
         tx_ol_cds.gr <- tx_ol_cds.gr[tx_ol_cds.gr$mid_exon_break]
      }

      breakpoint.gr <- GRanges()
      if (length(tx_ol_cds.gr) > 0){
         for (i_np in 1:length(tx_ol_cds.gr)){
            matching_cds.gr <- subsetByOverlaps(x = cds.tmp,
                                                ranges = tx_ol_cds.gr[i_np],
                                                ignore.strand = TRUE)
            cds_strand <- as.character(unique(strand(matching_cds.gr)))
            if (cds_strand == '+' ){

               message('TX ',unique(tx.gr$query),': Incomplete exon mapping found in ',
                       paste(unique(tx_ol_cds.gr[i_np]$break_gene),collapse=' and '),'.')
               i_breakpoint.gr <- GRanges(seqnames = seqnames(tx_ol_cds.gr[i_np]),
                                          ranges = IRanges(start = end(tx_ol_cds.gr[i_np]),
                                                           end = end(tx_ol_cds.gr[i_np]) + 1))
               mcols(i_breakpoint.gr) <- mcols(tx_ol_cds.gr[i_np])
               breakpoint.gr <- c(breakpoint.gr,i_breakpoint.gr)
            } else if ( cds_strand == '-'){
               message('TX ',unique(tx.gr$query),': Incomplete exon mapping found in ',
                       paste(unique(tx_ol_cds.gr[i_np]$break_gene),collapse=' and '),'.')

               i_breakpoint.gr <- GRanges(seqnames = seqnames(tx_ol_cds.gr[i_np]),
                                          ranges = IRanges(start = start(tx_ol_cds.gr[i_np]),end = start(tx_ol_cds.gr[i_np]) + 1))
               mcols(i_breakpoint.gr) <- mcols(tx_ol_cds.gr[i_np])
               breakpoint.gr <- c(breakpoint.gr,i_breakpoint.gr)
            }
         }
      }
   } else { # if no CDS matches found
      breakpoint.gr <- GRanges()
   }

   return(breakpoint.gr)
}

#breakpoint.gr <- c(unique(breakpoint.gr),)



protToExonSeqs <- function(cds.pa,proteins,prot_id){
   # prot_id = prot_ids[1]
   # cds.pa <- cds.pa
   gene_prot <- AAString()
   fullProtSeq <- proteins$protein_sequence[proteins$protein_id == prot_id]
   i_cds.pa <- unique(cds.pa[names(cds.pa) == prot_id,])
   i_cds.pa <- i_cds.pa[order(start(i_cds.pa),decreasing = FALSE)]
   prot_exons <- exons(x = edb.ucsc,filter = ProteinIdFilter(prot_id),columns = c('exon_idx'))

   if ( length(i_cds.pa) > 1 &&  length(i_cds.pa) > length(unique(mcols(i_cds.pa)[['exon_id']]))){
      rm_mask <- logical(length = length(i_cds.pa))
      for (i_np in 1:length(i_cds.pa)){
         i_i_cds.pa <- i_cds.pa[mcols(i_cds.pa)[['exon_id']] == mcols(i_cds.pa)[['exon_id']][i_np] ]
         rm_mask[i_np] <- width(i_cds.pa[i_np]) < max(width(i_i_cds.pa))
      }
      i_cds.pa <- i_cds.pa[!rm_mask, ]
   }

   #prot_seq_by_exon <- unlist(genomeToProtein(x = prot_exons,db = edb.ucsc))
   #prot_seq_by_exon <- prot_seq_by_exon[names(prot_seq_by_exon) == prot_id & start(prot_seq_by_exon) != -1,]
   n = length(i_cds.pa)
   np_ex_list <- ''
   if (n == 0){
      return(gene_prot)
   } else if (n > 1){
      gene_name <- unique(mcols(i_cds.pa)[['gene_name']])

      for (n_ex in 1:n){

         for (i_ex_start in 1:(n - n_ex + 1)){
            i_subset <- i_ex_start:(i_ex_start + n_ex - 1)
            new_prot <- AAStringSet(x = subseq(x = fullProtSeq,
                                               start = min(start(i_cds.pa[i_subset])),
                                               end = max(end(i_cds.pa[i_subset]))))
            names(new_prot) <- paste(gene_name,':x',paste(mcols(i_cds.pa)[['exon_rank']][i_subset],collapse='.'),'_',prot_id,sep='')
            if (length(gene_prot) > 0){gene_prot <- c(gene_prot,new_prot)} else {gene_prot <- new_prot}
         }


      }



   } else { # n = 1
      gene_name <- unique(mcols(i_cds.pa)[['gene_name']])
      np_start_pos <- start(i_cds.pa)
      np_end_pos <- end(i_cds.pa)
      gene_prot <- AAStringSet(x = subseq(x = fullProtSeq,start = np_start_pos,end = np_end_pos))
      names(gene_prot) <- paste(gene_name,':x',mcols(i_cds.pa)[['exon_rank']],'_',prot_id,sep = '')
   }

   return(gene_prot)
}


arriba <- function(chimeric.bam,aligned.bam,junction.csv,blacklist.file,annotation.file,ref.file,out.tsv,
                   overwrite.arriba = FALSE){
   if (!file.exists(out.tsv) | overwrite.arriba){
      out.dir <- paste(head(strsplit(out.tsv,'/')[[1]],-1),collapse='/')
      tmp_jct_file <- file.path(out.dir,'junction_tmp.tsv')
      jct_file_lines <- readLines(junction.csv)

      jcts_full.df <- read.delim(junction.csv,skip=grep('eventId',jct_file_lines) - 1,sep=',')
      junctions.df <- data.frame('pos1' = paste(convertChromToCharacter(ichr = jcts_full.df$chrA),':',jcts_full.df$posA,sep =''),
                                 'pos2' = paste(convertChromToCharacter(ichr = jcts_full.df$chrB),':',jcts_full.df$posB,sep =''),
                                 'dirA' = ifelse(sign(jcts_full.df$mposA)== 1,'-','+'),
                                 'dirB' = ifelse(sign(jcts_full.df$mposB)== 1,'-','+'))
      write.table(x = junctions.df,file = tmp_jct_file,row.names = FALSE,col.names = FALSE,quote=FALSE,sep='\t')

      arriba.cmd <- paste('/research/labs/experpath/vasm/shared/NextGen/tools/arriba/arriba -c',chimeric.bam,
                          '-x',aligned.bam,
                          '-g',annotation.file,
                          '-a',ref.file,
                          '-o',out.tsv,
                          '-b', blacklist.file,
                          '-d',tmp_jct_file,
                          '-f blacklist')
      system(arriba.cmd)
      unlink(tmp_jct_file)
   }

   arriba.df <- loadArriba(out.tsv)
   return(arriba.df)
}

loadArriba <- function(out.tsv){
   if (!is.na(out.tsv) & file.exists(out.tsv)){
      arriba.df <- read.delim(file = out.tsv,sep='\t')
      names(arriba.df) <- sub('X.','',names(arriba.df))
      arriba.df$genomic_breakpoint1 <- character(nrow(arriba.df))
      arriba.df$genomic_breakpoint2 <- character(nrow(arriba.df))
      for (i in 1:nrow(arriba.df)){
         arriba.df[i,'genomic_breakpoint1'] <- ifelse(arriba.df$closest_genomic_breakpoint1[i] == '.',
                                                      '',strsplit(arriba.df$closest_genomic_breakpoint1[i],'\\(')[[1]][1])
         arriba.df[i,'genomic_breakpoint2'] <- ifelse(arriba.df$closest_genomic_breakpoint2[i] == '.',
                                                      '',strsplit(arriba.df$closest_genomic_breakpoint2[i],'\\(')[[1]][1])
      }
      arriba.df$confidence <- factor(x = arriba.df$confidence, levels = c('high','medium','low'))
   } else {
      arriba.df <- data.frame()
   }
   return(arriba.df)
}

caseSummary <- function(caseId,sampleIdRNA,sampleIdDNA,tx.df,jct.df,orf.df,event.df,blat.gr,np.df,
                        geneExpression.df = NA,
                        arriba.df = NA){
   n <- nrow(tx.df) + sum(tx.df$np_output == '2 neopeptides') + sum(!unique(jct.df$eventId) %in% tx.df$DNA_eventId) + 20
   case.df <- data.frame('caseId' = rep(caseId,n),
                         'sampleIdRNA' = rep(sampleIdRNA,n),
                         'sampleIdDNA' = rep(sampleIdDNA,n),
                         'eventId' = integer(n),
                         'RNA_RPs' = integer(n),
                         'DNA_RPs' = integer(n),
                         'breakpointA' = character(n),
                         'breakpointB' = character(n),
                         'junctionType' = character(n),
                         'junctionSupport'= character(n),
                         'geneA' = character(n),
                         'readsGeneA' = integer(n),
                         'rpkmGeneA' = numeric(n),
                         'geneB' = character(n),
                         'readsGeneB' = integer(n),
                         'rpkmGeneB' = numeric(n),
                         'txName' = character(n),
                         'txSeq' = character(n),
                         'txLength' = integer(n),
                         'txScore' = numeric(n),
                         'txSpanContigs' = integer(n),
                         'txTPM' = character(n),
                         # 'txBlatCov' = numeric(n),
                         'orfOut' = character(n),
                         'orfGeneA' = character(n),
                         'orfGeneA_reads' = integer(n),
                         'orfExonHitA' = character(n),
                         'orfGeneB' = character(n),
                         'orfGeneB_reads' = integer(n),
                         'orfExonHitB' = character(n),
                         'orfKmers' = integer(n),
                         'orfPos' = character(n),
                         'orfSeq' = character(n),
                         'orfName' = character(n),
                         # 'orfFilter' = character(n),
                         'neopepSeq' = character(n),
                         'neopepName' = character(n),
                         'arribaPeptideMatch' = character(n),
                         'arribaGenes' = character(n),
                         'arribaSites' = character(n),
                         'arribaType' = character(n),
                         'arribaFrame' = character(n),
                         'arribaConf' = character(n),
                         'arribaSeq' = character(n))

   # tx_cols <- c('RNA_RPs','geneA','geneB','TX','tx_num')
   iRow <- 1
   event_list <- sort(unique(jct.df$eventId))
   # jct.df <- jct.df[!duplicated(jct.df$eventId),]
   for (iE in event_list){ # loop over iE ----
      jct_idx <- which(jct.df$eventId == iE & jct.df$Nassoc == max(jct.df$Nassoc[jct.df$eventId == iE]))[1]

      if (iE %in% event.df$DNA_eventId){
         e_spanning_RPs = event.df[event.df$DNA_eventId == iE,'RNA_RPs']
         e_DNA_RPs = event.df[event.df$DNA_eventId == iE,'DNA_RPs']
         e_bpA <- convertGenomicVectorToString(chr = event.df[event.df$DNA_eventId == iE,'chrA'],pos = event.df[event.df$DNA_eventId == iE,'DNA_posA'])
         e_bpB <- convertGenomicVectorToString(chr = event.df[event.df$DNA_eventId == iE,'chrB'],pos = event.df[event.df$DNA_eventId == iE,'DNA_posB'])
         e_geneA <- event.df[event.df$DNA_eventId == iE,'geneA']
         e_geneB <- event.df[event.df$DNA_eventId == iE,'geneB']
      } else {
         e_spanning_RPs = 0
         e_DNA_RPs <- max(jct.df[jct_idx,'Nassoc'])
         e_bpA <- convertGenomicVectorToString(chr = jct.df[jct_idx,'chrA'],pos = jct.df[jct_idx,'posA'])
         e_bpB <- convertGenomicVectorToString(chr = jct.df[jct_idx,'chrB'],pos = jct.df[jct_idx,'posB'])
         e_geneA <- paste(unique(jct.df[jct_idx,'geneA']),collapse=';')
         e_geneB <- paste(unique(jct.df[jct_idx,'geneA']),collapse=';')
      }
      junctionType = jct.df[jct_idx,'junctionType']
      junctionSupport = jct.df[jct_idx,'support']


      if (is.data.frame(geneExpression.df)){
         if (grepl(';',e_geneA)){
            e_geneA_reads <- sum(geneExpression.df[geneExpression.df$GeneName %in% strsplit(e_geneA,';')[[1]],'Count'])
            e_geneA_rpkm <- sum(geneExpression.df[geneExpression.df$GeneName %in% strsplit(e_geneA,';')[[1]],'RPKM_by_gene'])
         } else if (any(grepl(e_geneA,geneExpression.df$GeneName))) {
            e_geneA_reads <- sum(geneExpression.df[geneExpression.df$GeneName == e_geneA,'Count'])
            e_geneA_rpkm <- sum(geneExpression.df[geneExpression.df$GeneName == e_geneA,'RPKM_by_gene'])
         } else {
            e_geneA_reads <- 0
            e_geneA_rpkm = 0
         }

         if (grepl(';',e_geneB)){
            e_geneB_reads <- sum(geneExpression.df[geneExpression.df$GeneName %in% strsplit(e_geneB,';')[[1]],'Count'])
            e_geneB_rpkm <- sum(geneExpression.df[geneExpression.df$GeneName %in% strsplit(e_geneB,';')[[1]],'RPKM_by_gene'])
         } else if (any(grepl(e_geneB,geneExpression.df$GeneName))) {
            e_geneB_reads <- sum(geneExpression.df[geneExpression.df$GeneName == e_geneB,'Count'])
            e_geneB_rpkm <- sum(geneExpression.df[geneExpression.df$GeneName == e_geneB,'RPKM_by_gene'])
         } else {
            e_geneB_reads <- 0
            e_geneB_rpkm <- 0
         }
      } else {
         to_zero <- c('e_geneA_reads','e_geneA_rpkm','e_geneB_reads','e_geneB_rpkm')
         sapply(X = to_zero,FUN = 'assign',value = 0,pos = 1)
      }

      jct_arribaSeq = ''
      jct_arribaGenes = ''
      jct_arribaConf = ''
      jct_arribaSites = ''
      jct_arribaType = ''
      jct_arribaFrame = ''
      jct_arribaPeptideMatch = ''
      jct_arriba_idx <- integer()
      if (is.data.frame(arriba.df)){
         jct_arriba_idx <- which((arriba.df$genomic_breakpoint1 == e_bpA & arriba.df$genomic_breakpoint2 == e_bpB) |
                                    (arriba.df$genomic_breakpoint1 == e_bpB & arriba.df$genomic_breakpoint2 == e_bpA))
         if (length(jct_arriba_idx) == 1){
            jct_arribaGenes <- paste(arriba.df[jct_arriba_idx,c('gene1','gene2')],collapse='-->')
            jct_arribaConf <- as.character(arriba.df[jct_arriba_idx,'confidence'])
            jct_arribaSites <- paste(arriba.df[jct_arriba_idx,c('site1','site2')],collapse='-->')
            jct_arribaType <- arriba.df$type[jct_arriba_idx]
            jct_arribaFrame <- arriba.df$reading_frame[jct_arriba_idx]
            jct_arribaSeq <- arriba.df[jct_arriba_idx,'peptide_sequence']
         } else if (length(jct_arriba_idx) > 1){


            jct_arribaGenes <- character(length(jct_arriba_idx))
            jct_arribaConf <- character(length(jct_arriba_idx))
            jct_arribaSites <- character(length(jct_arriba_idx))
            jct_arribaType <- character(length(jct_arriba_idx))
            jct_arribaFrame <- character(length(jct_arriba_idx))
            jct_arribaSeq <- character(length(jct_arriba_idx))
            for (i_idx in 1:length(jct_arriba_idx)){
               jct_arribaGenes[i_idx] <- paste(arriba.df[jct_arriba_idx[i_idx],c('gene1','gene2')],collapse='-->')
               jct_arribaConf[i_idx] <- as.character(arriba.df[jct_arriba_idx[i_idx],'confidence'])
               jct_arribaSites[i_idx] <- paste(arriba.df[jct_arriba_idx[i_idx],c('site1','site2')],collapse='-->')
               jct_arribaFrame[i_idx] <- arriba.df$reading_frame[jct_arriba_idx[i_idx]]
               jct_arribaType[i_idx] <- arriba.df$type[jct_arriba_idx[i_idx]]
               jct_arribaSeq[i_idx] <- arriba.df[jct_arriba_idx[i_idx],'peptide_sequence']
            }
            jct_arribaGenes <- paste(jct_arribaGenes,collapse='|')
            jct_arribaConf <- paste(jct_arribaConf,collapse='|')
            jct_arribaSites <- paste(jct_arribaSites,collapse='|')
            jct_arribaType <- paste(jct_arribaType,collapse='|')
            jct_arribaFrame <- paste(jct_arribaFrame, collapse='|')
            jct_arribaSeq <- paste(jct_arribaSeq,collapse='|')
         }
      }


      e_tx.df <- tx.df[tx.df$DNA_eventId == iE,]

      if (nrow(e_tx.df) > 0){
         for (iTX in 1:nrow(e_tx.df)){ # loop over iTX ----
            # arriba_idx = integer(0)
            # arribaSeq <- ''
            # arribaPeptideMatch = ''
            # arribaGenes = ''
            # arribaConf = ''
            # arribaSites = ''
            # arribaType = ''
            # arribaFrame = ''



            txName <- e_tx.df[iTX,'tx_num']
            txSeq <- e_tx.df[iTX,'tx']
            txScore <- e_tx.df[iTX,'score']
            txLength <- nchar(txSeq)
            txSpanContigs <- e_tx.df[iTX,'nSpanningContigs']
            txTPM <- e_tx.df[iTX,'tpm']
            # t_blat.gr <- blat.gr[blat.gr$query == txName]
            # txBlat <- paste(seqnames(t_blat.gr),paste(ranges(t_blat.gr)),sep=':',collapse=';')
            # # foo <- subsetByOverlaps(x = unlist(cdsBy(x = edb.ucsc,filter = GRangesFilter(c(t_blat.gr,invertStrand(t_blat.gr))))),ranges = t_blat.gr,ignore.strand = TRUE)
            # if ('query_map_cov' %in% names(t_blat.gr)){
            #    txBlatCov <- unique(t_blat.gr$query_map_cov)
            # } else {
            #    txBlatCov <- sum(t_blat.gr$map_cov)
            # }

            t_orf.df <- orf.df[orf.df$name == e_tx.df[iTX,'tx_num'],]
            orfOut <- ifelse(!e_tx.df[iTX,'np_output'] %in% c('N/A',''),e_tx.df[iTX,'np_output'],'')
            if (nrow(t_orf.df) > 0){

               for (iORF in 1:nrow(t_orf.df)){
                  # arriba_idx = integer(0)
                  # arribaSeq <- ''
                  # arribaPeptideMatch = ''
                  # arribaGenes = ''
                  # arribaConf = ''
                  # arribaSites = ''
                  # arribaType = ''
                  # arribaFrame = ''


                  orfKmers <- t_orf.df[iORF,'n_kmers']
                  gene_str <- strsplit(t_orf.df[iORF,'gene'],'-')[[1]]
                  orfGeneA <- gene_str[1]
                  if (is.data.frame(geneExpression.df) && any(geneExpression.df$GeneName == orfGeneA)){
                     orfGeneA_reads <- sum(geneExpression.df[geneExpression.df$GeneName == orfGeneA[1],'Count'])
                  } else {orfGeneA_reads <- 0}
                  orfExonHitA <- strsplit(t_orf.df[iORF,'exons'],'-')[[1]][1]


                  if (length(gene_str) > 1){
                     orfGeneB <- gene_str[2]
                     if (is.data.frame(geneExpression.df) & any(geneExpression.df$GeneName == orfGeneB)){
                        orfGeneB_reads <- sum(geneExpression.df[geneExpression.df$GeneName == orfGeneB,'Count'])
                     } else {orfGeneB_reads <- 0}
                     orfExonHitB <- strsplit(t_orf.df[iORF,'exons'],'-')[[1]][2]
                  } else if (length(gene_str) == 1) {
                     orfGeneB = ''
                     orfExonHitB = ''
                     orfGeneB_reads = 0
                  }

                  orfSeq <- t_orf.df[iORF,'seq']
                  orfPos <- paste(t_orf.df[iORF,'strand'],':',t_orf.df[iORF,'startPosSign'],'-',t_orf.df[iORF,'endPosSign'],sep='')
                  # orfFilter <- t_orf.df[iORF,'ex_crit']
                  orfName <- t_orf.df[iORF,'name']
                  neopepSeq <- t_orf.df[iORF,'novel_seq']

                  t_np.df <- np.df[np.df$tx_name == e_tx.df[iTX,'tx_num'],]
                  if (nrow(t_np.df) == 1){
                     # neopepSeq <- t_np.df$Seq
                     neopepName <- t_np.df$Name
                  } else {
                     neopepName = ''
                  }

                  # arriba match? ----
                  arriba_idx <- integer(0)



                  if (is.data.frame(arriba.df)){
                     # arriba_idx <- which((arriba.df$genomic_breakpoint1 == e_bpA & arriba.df$genomic_breakpoint2 == e_bpB ) |
                     #                        (arriba.df$genomic_breakpoint1 == e_bpB &  arriba.df$genomic_breakpoint2 == e_bpA) )
                     arriba_idx <- which((arriba.df$genomic_breakpoint1 == e_bpA & arriba.df$gene1 == orfGeneA & arriba.df$genomic_breakpoint2 == e_bpB ) |
                                            (arriba.df$genomic_breakpoint1 == e_bpB & arriba.df$gene1 == orfGeneA & arriba.df$genomic_breakpoint2 == e_bpA) )
                  }

                  if (length(arriba_idx) == 1){
                     arribaGenes <- paste(arriba.df[arriba_idx,c('gene1','gene2')],collapse='-->')
                     arribaConf <- as.character(arriba.df[arriba_idx,'confidence'])
                     arribaSites <- paste(arriba.df[arriba_idx,c('site1','site2')],collapse='-->')
                     arribaType <- arriba.df$type[arriba_idx]
                     arribaFrame <- arriba.df$reading_frame[arriba_idx]
                     arribaSeq <- arriba.df[arriba_idx,'peptide_sequence']

                     peptide_sequence <- toupper(sub('\\*','',sub('\\|','',arriba.df[arriba_idx,'peptide_sequence'])))
                     if (peptide_sequence != '.'){
                        if (neopepSeq != ''){
                           pep_match <- countPattern(pattern = orfSeq,subject = peptide_sequence,max.mismatch = 2)
                           arribaPeptideMatch = 'YES'
                        } else if (nchar(orfSeq) >= nchar(peptide_sequence)){
                           pep_match <- countPattern(pattern = peptide_sequence,subject = subseq(orfSeq,start = nchar(orfSeq) - nchar(peptide_sequence)+1),max.mismatch = 2)
                           arribaPeptideMatch <- ifelse(pep_match > 0,'YES','NO')
                        } else {
                           pep_match <- countPattern(pattern = orfSeq,subject = subseq(peptide_sequence,start = nchar(peptide_sequence) - nchar(orfSeq)+1),max.mismatch = 2)
                           arribaPeptideMatch <- ifelse(pep_match > 0,'YES','NO')
                        }
                     } else {arribaPeptideMatch = 'NO'}


                     # if (arribaPeptideMatch == 'YES'){
                     #    arribaSeq <- arriba.df[arriba_idx,'peptide_sequence']
                     # }
                  } else if (length(arriba_idx) > 1){


                     peptide_sequence <- toupper(sub('\\*','',sub('\\|','',arriba.df[arriba_idx,'peptide_sequence'])))
                     peptide_match <- character(length(arriba_idx))
                     # arribaSeq <- arriba.df[arriba_idx,'peptide_sequence']
                     for (i_idx in 1:length(arriba_idx)){

                        if (peptide_sequence[i_idx] != '.'){
                           if (nchar(orfSeq) >= nchar(peptide_sequence[i_idx])){
                              pep_match <- countPattern(pattern = peptide_sequence[i_idx],subject = subseq(orfSeq,start = nchar(orfSeq) - nchar(peptide_sequence[i_idx])+1),max.mismatch = 2)
                              peptide_match[i_idx] <- ifelse(pep_match > 0,'YES','NO')
                           } else {
                              pep_match <- countPattern(pattern = orfSeq,subject = subseq(peptide_sequence[i_idx],start = nchar(peptide_sequence[i_idx]) - nchar(orfSeq)+1),max.mismatch = 2)
                              peptide_match[i_idx] <- ifelse(pep_match > 0,'YES','NO')
                           }
                        } else {
                           peptide_match[i_idx] <- 'NO'
                        }
                     }
                     if (sum(peptide_match == 'YES') == 1){
                        arriba_idx = arriba_idx[peptide_match == 'YES'][1]
                        arribaPeptideMatch <- 'YES'
                        arribaSeq <- arriba.df$peptide_sequence[arriba_idx][peptide_match == 'YES']
                        arribaGenes <- paste(arriba.df[arriba_idx,c('gene1','gene2')],collapse='-->')
                        arribaSeq <- arriba.df[arriba_idx,'peptide_sequence']
                        arribaConf <- as.character(arriba.df[arriba_idx,'confidence'])
                        arribaSites <- paste(arriba.df[arriba_idx,c('site1','site2')],collapse='-->')
                        arribaType <- arriba.df$type[arriba_idx]
                        arribaFrame <- arriba.df$reading_frame[arriba_idx]

                     } else if (sum(peptide_match == 'YES') > 1){
                        # message('Event ',iE,' TX ',txName,' ORF ',orfName,' has multiple Arriba hits')
                        arriba_idx = arriba_idx[peptide_match == 'YES'][1]

                        arriba_idx = arriba_idx[peptide_match == 'YES'][1]
                        arribaPeptideMatch <- 'YES'
                        arribaSeq <- arriba.df$peptide_sequence[arriba_idx]
                        arribaGenes <- paste(arriba.df[arriba_idx,c('gene1','gene2')],collapse='-->')
                        arribaSeq <- arriba.df[arriba_idx,'peptide_sequence']
                        arribaConf <- as.character(arriba.df[arriba_idx,'confidence'])
                        arribaSites <- paste(arriba.df[arriba_idx,c('site1','site2')],collapse='-->')
                        arribaType <- arriba.df$type[arriba_idx]
                        arribaFrame <- arriba.df$reading_frame[arriba_idx]


                     } else { # no peptide matches, multiple arriba_idx
                        arriba_idx <- arriba_idx[1]
                        arribaPeptideMatch <- 'NO'
                        arribaSeq <- arriba.df$peptide_sequence[arriba_idx]
                        arribaGenes <- paste(arriba.df[arriba_idx,c('gene1','gene2')],collapse='-->')
                        arribaConf <- as.character(arriba.df[arriba_idx,'confidence'])
                        arribaSites <- paste(arriba.df[arriba_idx,c('site1','site2')],collapse='-->')
                        arribaType <- arriba.df$type[arriba_idx]
                        arribaFrame <- arriba.df$reading_frame[arriba_idx]
                     }

                  } else {# cases: | arriba_idx | > 0
                     arriba_idx <- integer(0)
                     arribaPeptideMatch <- ''
                     arribaSeq <- jct_arribaSeq
                     arribaGenes <- jct_arribaGenes
                     arribaConf <- jct_arribaConf
                     arribaSites <- jct_arribaSites
                     arribaType <- jct_arribaType
                     arribaFrame <- jct_arribaFrame
                     arribaSeq <- jct_arribaSeq
                  }


                  case.df[iRow,'eventId'] <- iE
                  case.df[iRow,'RNA_RPs'] <- e_spanning_RPs
                  case.df[iRow,'DNA_RPs'] <- e_DNA_RPs
                  case.df[iRow,'junctionType'] <- junctionType
                  case.df[iRow,'junctionSupport'] <- junctionSupport
                  case.df[iRow,'breakpointA'] <- e_bpA
                  case.df[iRow,'breakpointB'] <- e_bpB
                  case.df[iRow,'geneA'] <- e_geneA
                  case.df[iRow,'readsGeneA'] <- e_geneA_reads
                  case.df[iRow,'rpkmGeneA'] <- e_geneA_rpkm
                  case.df[iRow,'geneB'] <- e_geneB
                  case.df[iRow,'readsGeneB'] <- e_geneB_reads
                  case.df[iRow,'rpkmGeneB'] <- e_geneB_rpkm
                  case.df[iRow,'txName'] <- txName
                  case.df[iRow,'txSeq'] <- txSeq
                  case.df[iRow,'txScore'] <- txScore
                  case.df[iRow,'txLength'] <- txLength
                  case.df[iRow,'txSpanContigs'] <- txSpanContigs
                  case.df[iRow,'txTPM'] <- txTPM
                  # case.df[iRow,'txBlatCov'] <- txBlatCov
                  case.df[iRow,'orfOut'] <- orfOut
                  case.df[iRow,'orfSeq'] <- orfSeq
                  case.df[iRow,'orfGeneA'] <- orfGeneA
                  case.df[iRow,'orfGeneA_reads'] <- orfGeneA_reads
                  case.df[iRow,'orfKmers'] <- orfKmers
                  case.df[iRow,'orfExonHitA'] <- orfExonHitA
                  case.df[iRow,'orfGeneB'] <- orfGeneB
                  case.df[iRow,'orfGeneB_reads'] <- orfGeneB_reads
                  case.df[iRow,'orfExonHitB'] <- orfExonHitB
                  case.df[iRow,'orfPos'] <- orfPos
                  # case.df[iRow,'orfFilter'] <- orfFilter
                  case.df[iRow,'orfName' ] <- orfName
                  case.df[iRow,'neopepSeq'] <- neopepSeq
                  case.df[iRow,'neopepName'] <- neopepName
                  case.df[iRow,'arribaGenes'] <- arribaGenes
                  case.df[iRow,'arribaConf'] <- arribaConf
                  case.df[iRow,'arribaSites'] <- arribaSites
                  case.df[iRow,'arribaFrame'] <- arribaFrame
                  case.df[iRow,'arribaType'] <- arribaType
                  case.df[iRow,'arribaSeq'] <- arribaSeq
                  case.df[iRow,'arribaPeptideMatch'] <- arribaPeptideMatch




                  iRow = iRow + 1
               } # loop over iORF in 1:nrow(t_orf.df)
            } else { #if nrow(t_orf.df) == 0
               case.df[iRow,'eventId'] <- iE
               case.df[iRow,'RNA_RPs'] <- e_spanning_RPs
               case.df[iRow,'DNA_RPs'] <- e_DNA_RPs
               case.df[iRow,'junctionType'] <- junctionType
               case.df[iRow,'junctionSupport'] <- junctionSupport
               case.df[iRow,'breakpointA'] <- e_bpA
               case.df[iRow,'breakpointB'] <- e_bpB
               case.df[iRow,'geneA'] <- e_geneA
               case.df[iRow,'readsGeneA'] <- e_geneA_reads
               case.df[iRow,'rpkmGeneA'] <- e_geneA_rpkm
               case.df[iRow,'geneB'] <- e_geneB
               case.df[iRow,'readsGeneB'] <- e_geneB_reads
               case.df[iRow,'rpkmGeneB'] <- e_geneB_rpkm

               case.df[iRow,'txName'] <- txName
               case.df[iRow,'txSeq'] <- txSeq
               case.df[iRow,'txLength'] <- txLength
               case.df[iRow,'txScore'] <- txScore
               case.df[iRow,'txSpanContigs'] <- txSpanContigs

               case.df[iRow,'txTPM'] <- txTPM
               # case.df[iRow,'txBlatCov'] <- txBlatCov
               case.df[iRow,'orfOut'] <- orfOut
               # case.df[iRow,'orfSeq'] <- ''
               # case.df[iRow,'orfGeneA'] <- ''
               # case.df[iRow,'orfGeneA_reads'] <- ''
               # case.df[iRow,'orfKmers'] <- 0
               # case.df[iRow,'orfExonHitA'] <- ''
               # case.df[iRow,'orfGeneB'] <- ''
               # case.df[iRow,'orfGeneB_reads'] <- 0
               # case.df[iRow,'orfExonHitB'] <- ''
               # case.df[iRow,'orfPos'] <- ''
               # case.df[iRow,'orfFilter'] <- ''
               # case.df[iRow,'orfName' ] <- ''
               # case.df[iRow,'neopepSeq'] <- ''
               # case.df[iRow,'neopepName'] <- ''

               case.df[iRow,'arribaGenes'] <- jct_arribaGenes
               case.df[iRow,'arribaConf'] <- jct_arribaConf
               case.df[iRow,'arribaSites'] <- jct_arribaSites
               case.df[iRow,'arribaFrame'] <- jct_arribaFrame
               case.df[iRow,'arribaType'] <- jct_arribaType
               case.df[iRow,'arribaSeq'] <- jct_arribaSeq

               iRow = iRow + 1
            } # if nrow(e_tx.df) > 0 and nrow(e_orf.df) == 0
         } # loop over iTX
      } else if (nrow(e_tx.df) == 0 & length(jct_arriba_idx) > 0) { # if nrow(e_tx.df) == 0
         case.df[iRow,'eventId'] <- iE
         case.df[iRow,'RNA_RPs'] <- e_spanning_RPs
         case.df[iRow,'DNA_RPs'] <- e_DNA_RPs
         case.df[iRow,'junctionType'] <- junctionType
         case.df[iRow,'junctionSupport'] <- junctionSupport
         case.df[iRow,'breakpointA'] <- e_bpA
         case.df[iRow,'breakpointB'] <- e_bpB
         case.df[iRow,'geneA'] <- e_geneA
         case.df[iRow,'readsGeneA'] <- e_geneA_reads
         case.df[iRow,'rpkmGeneA'] <- e_geneA_rpkm
         case.df[iRow,'geneB'] <- e_geneB
         case.df[iRow,'readsGeneB'] <- e_geneB_reads
         case.df[iRow,'rpkmGeneB'] <- e_geneB_rpkm


         case.df[iRow,'arribaGenes'] <- jct_arribaGenes
         case.df[iRow,'arribaConf'] <- jct_arribaConf
         case.df[iRow,'arribaSites'] <- jct_arribaSites
         case.df[iRow,'arribaFrame'] <- jct_arribaFrame
         case.df[iRow,'arribaType'] <- jct_arribaType
         case.df[iRow,'arribaSeq'] <- jct_arribaSeq

         iRow = iRow + 1
      } else { # if nrow(e_tx.df) == 0

         case.df[iRow,'eventId'] <- iE
         case.df[iRow,'RNA_RPs'] <- e_spanning_RPs
         case.df[iRow,'DNA_RPs'] <- e_DNA_RPs
         case.df[iRow,'junctionType'] <- junctionType
         case.df[iRow,'junctionSupport'] <- junctionSupport
         case.df[iRow,'breakpointA'] <- e_bpA
         case.df[iRow,'breakpointB'] <- e_bpB
         case.df[iRow,'geneA'] <- e_geneA
         case.df[iRow,'readsGeneA'] <- e_geneA_reads
         case.df[iRow,'rpkmGeneA'] <- e_geneA_rpkm
         case.df[iRow,'geneB'] <- e_geneB
         case.df[iRow,'readsGeneB'] <- e_geneB_reads
         case.df[iRow,'rpkmGeneB'] <- e_geneB_rpkm
         iRow = iRow + 1
      }
   } # loop over iE
   case.df <- case.df[1:iRow,]
   return(case.df)
}


plotJunctionTX <- function(tx.gr,e_reads.gr,pos.gr,out.file = NA,
                           breakpoints.gr = NA,
                           protein_coding_only = FALSE,highlight = FALSE,collapse_tx = 'longest',
                           plot_reads = FALSE){
   if (is.na(out.file)){write_to_pdf = FALSE} else {write_to_pdf = TRUE}
   # debug: tx.gr = tx.gr; pos.gr = A.gr; protein_coding_only = FALSE;  highlight = TRUE;  collapse_tx = 'longest'
   options(ucscChromosomeNames = TRUE)
   if (any(seqlevelsStyle(tx.gr) != 'UCSC')){seqlevelsStyle(tx.gr) <- 'UCSC'}
   if (any(seqlevelsStyle(pos.gr) != 'UCSC')){seqlevelsStyle(pos.gr) <- 'UCSC'}
   tx_ol.gr <- subsetByOverlaps(x = tx.gr,ranges = pos.gr,ignore.strand = TRUE)



   if (protein_coding_only){
      gr <- getGeneRegionTrackForGviz(edb.ucsc,
                                      filter = GeneBiotypeFilter('protein_coding'),
                                      chromosome = as.character(seqnames(pos.gr)),
                                      start = start(pos.gr),end = end(pos.gr),
                                      featureIs = 'gene_biotype')
   } else {

      gr <- getGeneRegionTrackForGviz(edb.ucsc,
                                      chromosome = seqnames(pos.gr),
                                      start = start(pos.gr),end = end(pos.gr),
                                      featureIs = 'gene_biotype')
   }

   #seqlevelsStyle(gr) <- 'UCSC'
   #collapse_transcripts_var <- collapse_tx

   grt <- GeneRegionTrack(range = gr,
                          exon = gr$exon_rank,
                          feature = gr$feature,
                          strand = gr$strand,
                          transcript = gr$transcript,
                          gene = gr$gene,
                          symbol = gr$symbol,
                          name = paste('Gene Region'),
                          transcriptAnnotation = 'symbol',shape = 'smallArrow',
                          collapseTranscripts = 'longest')

   gat <- GenomeAxisTrack()

   idt <- IdeogramTrack(chromosome = as.character(seqnames(pos.gr)),genome = 'hg38')

   if (length(tx_ol.gr) == 0){
      message('No TX map to this region.')

   } else if (length(tx_ol.gr) > 0) {


      new_ids <- sapply(X = strsplit(tx_ol.gr$query,split = paste(event,'_',sep=''))
                        ,FUN = '[',2)


      #merged_var <- ifelse(grepl('merged',tx_ol.gr$query),'merged','kmer')
      txTrack <- AnnotationTrack(range = tx_ol.gr,
                                 name = 'Assembled Transcripts',
                                 id = new_ids,
                                 strand = tx_ol.gr$strand,
                                 group = new_ids,
                                 groupAnnotation = 'id',
                                 genome =  'hg38',
                                 fontsize.group=12,
                                 featureAnnotation = 'feature',
                                 transcriptAnnotation = 'symbol',shape = 'smallArrow',
                                 showFeatureId = FALSE,
                                 match = 'red',no_match = 'gray')#,
      #feature = merged_var)
   }

   if (length(e_reads.gr) == 0){
      message('No read mappings found.')
   } else if (plot_reads) {

      new_ids <- sapply(strsplit(x = e_reads.gr$query,split = '_'),FUN='[',3)

      e_reads.gr$flag <- logical(length = length(e_reads.gr))

      readsTrack <- AlignmentsTrack(range = e_reads.gr,
                                    strand = rep('*',length(e_reads.gr)),
                                    genome = 'hg38',
                                    id = new_ids,
                                    type = 'coverage',
                                    stackHeight = 0.1,
                                    coverageHeight = 0.08,
                                    minCoverageHeight = 0,
                                    name = 'RNA spanRP Coverage')

   }

   cds_hits <- subsetByOverlaps(ranges = tx_ol.gr,x = gr[gr$feature == 'protein_coding'],ignore.strand = TRUE)
   #exons <- exonsByOverlaps(x = gr,ranges = )

   #filter = TxBiotypeFilter('protein_coding'))

   obj_start <- min(min(start(tx_ol.gr)),min(start(gr)))
   #obj_end <- max(max(end(pos.gr)),max(end(tx_ol.gr))+ width(range(tx_ol.gr,ignore.strand = TRUE))*.2)
   obj_end <- max(end(pos.gr),max(end(tx_ol.gr)))
   title_delta <- (obj_end - obj_start ) * .2
   plot_start <- obj_start - title_delta
   plot_end <- obj_end + title_delta

   # Assemble track list & plot ----
   trackList <- list(idt,gat)

   if ( plot_reads) { trackList <- c(trackList,readsTrack)}
   if (highlight && length(cds_hits) > 0){

      htrack <- HighlightTrack(trackList = list(txTrack,grt),start = start(cds_hits),end = end(cds_hits),chromosome =  seqnames(cds_hits),
                               col = "#98e9d6",fill="#98e9d6")
      trackList <- c(trackList,htrack)
   }  else if (length(tx_ol.gr) > 0){
      trackList <- c(trackList,txTrack,grt)
   } else {
      trackList <- c(trackList,grt)
   }

   idtSize <- 0.075; gatSize = 0.1; rtSize <- 0.3; grtSize <- length(unique(gr$symbol))*0.1
   txSize <- 1-idtSize - gatSize - rtSize - grtSize
   trackSizes <- c(idtSize,gatSize,rtSize,txSize,grtSize)
   plotTracks(trackList,sizes = trackSizes,
              from = plot_start,to = plot_end,chromosome = seqlevels(pos.gr))

}

runTxViz <- function(event.dir,event.tx.fa = NA,overwrite_plots = FALSE){
   event.blat.gr <- makeGRangesFromDataFrame(df = read.table(file = file.path(event.dir,'good.transcripts.alignments.tsv'),sep=' ',stringsAsFactors = FALSE,header = TRUE),
                                             keep.extra.columns = TRUE)
   if (!exists('event.tx.fa') || is.null(event.tx.fa)){event.tx.fa <- file.path(event.dir,'good.transcripts.fa')}
   jct.gr <- makeGRangesFromDataFrame(df = read.table(file.path(event.dir,'event.txt'),sep=' '))
   if (overlapsAny(jct.gr[1],jct.gr[2])){jct.gr <- range(jct.gr)}
   event.bam <- file.path(event.dir,'event.bam')
   span.bam <- file.path(event.dir,'span.bam')

   if (length(jct.gr) == 1){
      out.pdf <- file.path(event.dir,'region.pdf')
      if (!file.exists(out.pdf) | overwrite_plots){
         txViz(tx.gr = event.blat.gr,tx.fa = event.tx.fa,pos.gr = jct.gr,all.bam = event.bam,span.bam = span.bam,out.pdf = out.pdf)
      }


   } else if (length(jct.gr) == 2){
      outA.pdf <- file.path(event.dir,'regionA.pdf')
      if (!file.exists(outA.pdf) | overwrite_plots){
         txViz(tx.gr = event.blat.gr,tx.fa = event.tx.fa,pos.gr = jct.gr[1],all.bam = event.bam,span.bam = span.bam,out.pdf = outA.pdf)
      }

      outB.pdf <- file.path(event.dir,'regionB.pdf')
      if (!file.exists(outB.pdf) | overwrite_plots){
         txViz(tx.gr = event.blat.gr,tx.fa = event.tx.fa,pos.gr = jct.gr[2],all.bam = event.bam,span.bam = span.bam,out.pdf = outB.pdf)
      }


   }

}

txViz <- function(tx.gr,tx.fa,pos.gr,all.bam,span.bam,out.pdf){
   # aligned.bam = '/research/labs/experpath/vasm/shared/NextGen/Projects/Neoantigens/RNAseq/NE48501/GRCh38/intruder/tx/event_145/aligned.bam'
   # chimeric.bam = '/research/labs/experpath/vasm/shared/NextGen/Projects/Neoantigens/RNAseq/NE48501/GRCh38/intruder/tx/event_145/chimeric.bam'
   # tx.fa = fa_for_viz
   options(ucscChromosomeNames = TRUE)
   if (any(seqlevelsStyle(pos.gr) != 'UCSC')){seqlevelsStyle(pos.gr) <- 'UCSC'};
   if (any(seqlevelsStyle(tx.gr) != 'UCSC')){seqlevelsStyle(tx.gr) <- 'UCSC'};


   #reads.gr <- as(readGAlignmentPairs(file = reads.bam,use.names=TRUE),'GRanges')
   alignedReadsTrack <- AlignmentsTrack(range = all.bam,genome = 'GRCh38',name = 'All reads',type = 'coverage')
   spanningReadsTrack <- AlignmentsTrack(range = span.bam,genome = 'GRCh38',name = 'Spanning reads',type = 'coverage')
   #readsTrack <- AlignmentsTrack(range = reads.bam)




   tx_ol.gr <- subsetByOverlaps(x = tx.gr,ranges = pos.gr,ignore.strand = TRUE)
   suppressWarnings(cds.gr <- cdsBy(x = edb.ucsc,filter = GRangesFilter(value = pos.gr), by = 'gene'))

   if (!is.null(cds.gr)){
      suppressWarnings(gr <- getGeneRegionTrackForGviz(edb.ucsc,
                                                       chromosome = as.character(seqnames(pos.gr)),
                                                       start = start(pos.gr),end = end(pos.gr),
                                                       featureIs = 'gene_biotype'))
   } else {gr <- GRangesList()}



   grt <- GeneRegionTrack(range = gr,
                          exon = gr$exon_rank,
                          feature = gr$feature,
                          strand = gr$strand,
                          transcript = gr$transcript,
                          gene = gr$gene,
                          symbol = gr$symbol,
                          name = paste('Gene Region'),
                          transcriptAnnotation = 'symbol',shape = 'smallArrow',
                          collapseTranscripts = 'longest')

   gat <- GenomeAxisTrack()

   idt <- IdeogramTrack(chromosome = as.character(seqnames(pos.gr)),genome = 'hg38')

   txTrack <- AnnotationTrack(range = tx_ol.gr,
                              name = 'Assembled Transcripts',
                              stacking = 'full',
                              id = tx_ol.gr$query,
                              strand = tx_ol.gr$strand,
                              group = tx_ol.gr$query,
                              groupAnnotation = 'id',
                              genome =  'hg38',
                              fontsize.group=12,
                              featureAnnotation = 'feature',
                              transcriptAnnotation = 'symbol',shape = 'arrow',
                              showFeatureId = FALSE,
                              match = 'red',no_match = 'gray')#,
   trackSizes = c(.075,.075,
                  .2,.2,
                  .3,
                  .15)
   trackList <- list(idt,gat,alignedReadsTrack,spanningReadsTrack,txTrack,grt)
   pdf(out.pdf)
   plotTracks(trackList,sizes = trackSizes,extend.right = 3000,extend.left = 3000)
   dev.off()
   # plotTracks(trackList,sizes = trackSizes)#,
   #            #from = plot_start,to = plot_end,
   #            chromosome = seqlevels(pos.gr))
}

getProtSeq <- function(diamond.df,event,out_fa,overwrite_ref = TRUE){
   e_diamond.df <- diamond.df[grepl(pattern = paste0('J_',event,'_'),x = diamond.df$query),]
   if (nrow(e_diamond.df) > 0 & (!file.exists(out_fa) | overwrite_ref)){
      e_prot_ids <- sapply(X = e_diamond.df$target,FUN = function(i) strsplit(grep('ENSP',unlist(strsplit(i,'\\|')),value=TRUE),'\\.')[[1]][1],USE.NAMES = FALSE)
      e_gene_names <- sapply(X = e_diamond.df$target,FUN = function(i) unlist(strsplit(i,'\\|'))[7],USE.NAMES = FALSE)
      proteins <- proteins(edb.ucsc,columns = c("gene_name","protein_id","protein_sequence"),
                           filter = ProteinIdFilter(e_prot_ids))
      prot.ss <- AAStringSet(x = proteins$protein_sequence); names(prot.ss) <-paste(proteins$protein_id,proteins$gene_name,sep='|')
      writeXStringSet(x = prot.ss,filepath = out_fa,format='fasta')
   } else if (nrow(e_diamond.df) == 0 & file.exists(out_fa)){
      unlink(out_fa)
   }
}

getEventFromTxNum <- function(tx_num){
   return(as.integer(strsplit(tx_num,'_')[[1]][2]))

}

runGMAP <- function(input.fa,gmap.out.dir){
   gmap_genome_dir <- '/research/labs/experpath/vasm/shared/Genome/Human/referenceFiles/gmap/hg38'
   if (!dir.exists(gmap.out.dir)){dir.create(gmap.out.dir,recursive = TRUE)}
   gmap.out <- file.path(gmap.out.dir,'gmap.out')
   gmap.cmd <- paste('conda run -n gmap_sandbox gmap -d hg38 -P','-t',numCores,'-F -T -A -f 2 --nofails',input.fa,'>',gmap.out)
   system(gmap.cmd)
   gmap.gff3 <- read.table(file = gmap.out,sep='\t',col.names = c('seqid','source','type','start','end','score','strand','phase','attributes'))
   attribute_list <- c('ID','Name','Parent','Target','Dir','coverage','identity','matches','mismatches','indels','unknowns')
   for (att in attribute_list){
      gmap.gff3[,att] <- sapply(X = gmap.gff3$attributes,FUN = function(i) sub(paste0(att,'='),'',grep(paste0(att,'='),unlist(strsplit(i,';')),value=TRUE)))
   }
   # gmap.gff3$ID <- sapply(X = gmap.gff3$attributes,FUN = function(i) sub('ID=','',grep('ID=',unlist(strsplit(i,';')),value=TRUE)))
   # gmap.gff3$name <- sapply(X = gmap.gff3$attributes,FUN = function(i) sub('Name=','',grep('Name=',unlist(strsplit(i,';')),value=TRUE)))
}

runDiamond <- function(tx.dir,input.fa){
   # Align against reference TX ----


   tx_ref_file <- '/research/labs/experpath/vasm/shared/Genome/Human/referenceFiles/transcriptome/gencode.v39.transcripts.fa'
   diamond_gencode_ref_path <- '/research/labs/experpath/vasm/shared/NextGen/judell/data/Reference/gencode.v39.prot'
   tx_ref_2bit <- '/research/labs/experpath/vasm/shared/Genome/Human/referenceFiles/transcriptome/gencode.v39.transcripts.2bit'
   diamond.dir <- file.path(tx.dir,'diamond')
   if (!dir.exists(diamond.dir)){dir.create(diamond.dir)}
   diamond.outfile <- file.path(diamond.dir,'matches.tsv')
   diamond.cmd <- paste('diamond blastx -f 6 qseqid sseqid pident length mismatch gapopen qstrand qframe qlen qstart qend qcovhsp slen sstart send scovhsp evalue bitscore full_sseq qseq_translated',
                        '-iterate',#'-e 1e-8 --top 0.1',
                        '-d',diamond_gencode_ref_path,
                        '-q',input.fa,
                        '-o', diamond.outfile)
   system(diamond.cmd,ignore.stdout = TRUE, ignore.stderr =  TRUE)

   #  Additional Output Options:
   # qseq means Aligned part of query sequence
   # qseq_translated means Aligned part of query sequence (translated)
   # full_qseq means Query sequence

   blast6_header <- c('query','target','percent_id','alignment_length','mismatches','gap_o','query_strand','query_frame','query_length',
                      'query_start','query_end','query_cov_HSP','target_length','target_start','target_end','target_cov_HSP','E_value','bit_score','target_seq','trans_query_seq')
   diamond.df <- read.table(diamond.outfile,header = FALSE,sep = '\t',stringsAsFactors = FALSE,col.names = blast6_header)
   diamond.df$prot_id <- sapply(X = diamond.df$target,FUN = function(x) strsplit(strsplit(x,'\\|')[[1]][1],'\\.')[[1]][1])
   diamond.df$gene_name <- sapply(X = diamond.df$target,FUN = function(x) strsplit(x,'\\|')[[1]][7])
   diamond.df$gp_match <- paste(diamond.df$prot_id,diamond.df$gene_name,sep='|')
   filterDiamond = FALSE
   if (filterDiamond){
      diamond_wt_tx.fa <- file.path(diamond.dir,'diamond_wt.fa')
      diamond_wt_list_file <- file.path(diamond.dir,'diamond_wt.list')

      diamond_wt.df <- diamond.df[diamond.df$target_cov_HSP >= 90 & diamond.df$query_cov_HSP < 90,]

      writeLines(diamond_wt_list,diamond_wt_list_file)
      top_hits.cmd <- paste('seqtk subseq',unfiltered_tx.fa,diamond_wt_list_file,'>',diamond_wt_tx.fa)
      system(top_hits.cmd)

      diamond_out_tx.fa <- file.path(tx.dir,'diamond_tx.fa')
      fw_to_write <- diamond.df$query[diamond.df$query_strand == '+']
      fw_list_file <- file.path(diamond.dir,'diamond_fw.list')
      writeLines(fw_to_write,con=fw_list_file)
      rev_to_write <- diamond.df$query[diamond.df$query_strand == '-']
      rev_list_file <- file.path(diamond.dir,'diamond_rev.list')

      writeLines(rev_to_write,con=rev_list_file)
      rev.fa <- file.path(diamond.dir,'rc.fa')
      diamond_rc.cmd <- paste('seqtk subseq',unfiltered_tx.fa,rev_list_file,'| seqtk seq -r - >' ,rev.fa)
      system(diamond_rc.cmd)

      fwd.fa <- file.path(diamond.dir,'fw.fa')
      diamond_fw.cmd <- paste('seqtk subseq',unfiltered_tx.fa,fw_list_file,'>',fwd.fa)#,' || cat',fwd.fa,rev.fa,' >',diamond_out_tx.fa)
      system(diamond_fw.cmd)
      diamond_cat.cmd <- paste('cat',fwd.fa,rev.fa,'>',diamond_out_tx.fa)
      system(diamond_cat.cmd)

      # diamond_out_tx.fa <- diamond_out_tx.fa
   }
   return(diamond.df)
}
