# libraries -----

library(parallel)
library(plotrix)

#' @param readsfilename file created by rnaSupportForDNAseqJcts
#' @param file.prefix caseId
#' @param tx.dir path to the folder for the output
#' @param optimize_velet minimum number of reads to run pipeline
#' 

# readsFilename,out.dir,file.prefix,insert_length,
# flash.threshold = 50,
# k = seq(15,31,by=4),
# reads.threshold = 0,
# cov_cutoff = 'dynamic',
# take.rc = FALSE,
# cleanup = FALSE,
# overwrite_tx = TRUE,
# optimize_velvet = TRUE

spec = matrix(c(
   'reads.file', 'f',"character","location of file containing RNA spanning reads pairs",NA,
   'k','k','character','hash lengths for Velveth',paste(seq(15,31,by=4),collapse=','),
   'reads.threshold','r','integer','minimum number of spanning read pairs required for assembly',0,
   'insert.length','i','numeric','insert length of RNA library',NA,
   'cov_cutoff','c','character','coverage cutoff for velvet','dynamic',
   'flash.threshold','fl','integer','number of read pairs to trigger quick combining of read pairs',50),
   byrow = TRUE,ncol = 5)

debug = FALSE
if (debug){
   readsFilename = '/research/labs/experpath/vasm/shared/NextGen/Projects/ExVivoEndom/RNAseq/EM28502/GRCh38/svar-1/dev/EM102_RNA_supportForDNAseqJcts_reads_EM102_MP.csv'
   tx.dir = '/research/labs/experpath/vasm/shared/NextGen/Projects/ExVivoEndom/RNAseq/EM28502/GRCh38/svar-1/dev/tx'
   file.prefix = 'EM102_RNA'
   k = paste(seq(15,31,by=4),collapse=',')
   reads.threshold = 0
   cov_cutoff = 'dynamic'
   flash.threshold = 50
   insert_length = NULL
}

velvet_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/velvet/1.2.10'
oases_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/oases'
flash_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/FLASH-1.2.11/flash'
perl_path <- '/research/labs/experpath/vasm/shared/NextGen/judell/lib/perl5/perls/perl-ju/bin/perl'
vo_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/VelvetOptimiser/2.2.6/VelvetOptimiser.pl'
old_path <- Sys.getenv("PATH")
if (!grepl('velvet',old_path)){Sys.setenv(PATH = paste(old_path, velvet_path, sep = ":"))}
rw.pause.time <- 0.01

# Functions -----

intruderAssembly <- function(reads.filename,out.dir,file.prefix,
                             insert_length = NA,
                             flash.threshold = 50,
                             k = seq(15,31,by=4),
                             reads.threshold = 0,
                             cov_cutoff = 'dynamic',
                             # take.rc = FALSE,
                             cleanup = FALSE,
                             optimize_velvet = TRUE) {
   
   # Setup variables ----
   # k = seq(15,31,by=4); COV_CUTOFF = 'dynamic'; reads.threshold = 0;
   #   flash.threshold = 50;  overwrite_tx = TRUE; optimize_velvet = TRUE
   COV_CUTOFF <- cov_cutoff
   seq.df <- loadRNAReads(reads.filename)
   message(paste('Loaded ',nrow(seq.df), ' spanning RNA RPs','...',sep=''))
   
   if (nrow(seq.df) == 0){
      message('No RNA reads found for any junctions.')
      return(NULL)
   }
   
   # output[['seq.df']] <- seq.df
   
   
   if (!dir.exists(out.dir)){dir.create(out.dir)}
   
   tx.df.filename <- file.path(out.dir,paste(file.prefix,'TX.csv',sep='_'))
   
   event_list <- sort(unique(seq.df$DNA_eventId[seq.df$RNA_RPs >= reads.threshold]))
   
   if (length(event_list) == 0){message('No events meet reads threshold.')}
   
   check_prot <- FALSE
   
   
   output_tx <- data.frame()#data.frame(unique(inputSeqDF[,!names(inputSeqDF) %in% c('seqIdNum','ReadA','ReadB')]))
   df_header <- c('DNA_sampleId','DNA_eventId','RNA_RPs','chrA','chrB','DNA_posA','DNA_posB','geneA','geneB')
   
   
   
   
   # Loop over event list ------
   for (event in event_list) {
      
      e_seq.df <- seq.df[seq.df$DNA_eventId == event,]
      event.tx.dir <- file.path(out.dir,paste('event',event,sep='_'))
      if (!dir.exists(event.tx.dir)){dir.create(event.tx.dir)}
      
      
      n_spann_RPs <- nrow(e_seq.df)
      
      
      # Pre-process event reads ---- 
      if (!optimize_velvet){prehash_reads = TRUE} else {prehash_reads = FALSE}
      readsPP.out <- spanningRPsPreProcessing(seq.df = seq.df,event = event,
                                              event.tx.dir = event.tx.dir,
                                              file.prefix = file.prefix,
                                              flash.threshold = flash.threshold,
                                              remove.polyA = FALSE,
                                              prehash_reads = prehash_reads)
      #toc()
      flash.flag = readsPP.out$n_combined > 0
      # test assembly ----
      tx <- velvetOases(iJct = event,kmer = 15,
                        event.tx.dir = event.tx.dir,
                        file.prefix = file.prefix,
                        flash.flag = flash.flag)
      
      velvet_check <- length(tx) > 0 | n_spann_RPs > 15
      unlink(file.path(event.tx.dir,'k_15'),recursive = TRUE)
      
      if (nrow(e_seq.df)>1){
         e_seq.df$posAmin = min(e_seq.df$posAmin)
         e_seq.df$posAmax = min(e_seq.df$posAmax)
         e_seq.df$posBmin = min(e_seq.df$posBmin)
         e_seq.df$posBmax = min(e_seq.df$posBmax)
      }
      unique_e_seq.df <- unique(e_seq.df[,df_header])
      
      # event.readPairs.fa <- readsPP.out[['event.readPairs.fa']]
      # event.flashed.reads.fa <- readsPP.out[['event.flashed.reads.fa']]
      if (COV_CUTOFF == 'dynamic'){
         cov_cutoff = max(0,0.024*n_spann_RPs - 0.2*log10(n_spann_RPs))
      } else {cov_cutoff = COV_CUTOFF}
      
      #cov_cutoff <- min(2,floor(n_spann_RPs / 5))
      # toggle VO -----
      #optimize_velvet = TRUE
      
      if (velvet_check & optimize_velvet){
         #tic()
         # Velvet Optimiser ----
         if (any(grepl('k',list.dirs(event.tx.dir)))){
            unlink(grep('k',list.dirs(event.tx.dir),value=TRUE),recursive = TRUE)
         }
         
         
         tic()
         
         tx <- velvetOptimiser(event.tx.dir = event.tx.dir,flash.flag = flash.flag,file.prefix = file.prefix,
                               insert_length = insert_length,oases_cov_cutoff = cov_cutoff)
         vo.toc <- toc(quiet=TRUE)
         
         # tic()
         if (length(tx) > 0){
            vo_out.dir <- file.path(event.tx.dir,'vo')
            tx.info <- strsplit(names(tx),'_')
            tx.names <- strsplit(names(tx),'/')
            if (!dir.exists(vo_out.dir)){
               velvet_out_dir <- file.path(event.tx.dir,grep('k',list.files(event.tx.dir),value=TRUE))
               kmer_size = strsplit(grep('k',list.files(event.tx.dir),value=TRUE),'_')[[1]][2]
               data <- read.table(file.path(velvet_out_dir,'stats.txt'),header = TRUE)
               exp_cov <- weighted.mean(data$short1_cov, data$lgth)
               
               contig.out <- contigParser(tx_out_dir = velvet_out_dir,merged = FALSE)
            } else {
               vo_log.file <- file.path(vo_out.dir,grep('logfile',list.files(vo_out.dir),value = TRUE))
               vo_log <- readLines(vo_log.file)
               
               final_assembly_info = vo_log[(grep('Final optimised assembly details:',vo_log) + 2):length(vo_log)]
               kmer_size = sub('Velvet hash value: ','',grep('Velvet hash value',final_assembly_info,value = TRUE))
               velvetg_run_param <- grep('Velvetg parameter',final_assembly_info,value=TRUE)
               exp_cov <- unlist(regmatches(velvetg_run_param,gregexpr('(?<=-exp_cov )[\\.\\d]+',velvetg_run_param,perl=TRUE)))
               cov_cutoff <- unlist(regmatches(velvetg_run_param,gregexpr('(?<=-cov_cutoff )[\\.\\d]+',velvetg_run_param,perl=TRUE)))
               
               contig.out <- contigParser(tx_out_dir = vo_out.dir,merged = FALSE)
            }
            
            # message(length(tx),' transcript(s) assembled (hash length: ',kmer_size,') in ',format(vo.toc$toc-vo.toc$tic,nsmall = 1),' seconds.')
            
            
            
            reads_tx.hash <- contig.out$read_tx.hash
            for (i_tx in 1:length(tx)){
               
               n_reads <- sum(reads_tx.hash[,i_tx])
               total_used_reads <- sum(rowSums(reads_tx.hash) > 0)
               tx_conf <- as.numeric(tx.info[[i_tx]][which(grepl('Conf',tx.info[[i_tx]]))+1])
               
               df.data <- data.frame(unique_e_seq.df,
                                     'kmer_size' = kmer_size,
                                     'cov_cutoff' = cov_cutoff,
                                     'exp_cov' = exp_cov,
                                     'TX' = toString(tx[i_tx]),
                                     'tx_num' = tx.names[[i_tx]][1],
                                     'tx_conf' = tx_conf,
                                     'tx_length' = width(tx[i_tx]),
                                     'n_reads' = n_reads,
                                     # 'seq_list' = reads_list,
                                     # 'total_event_reads' = 2*n_spann_RPs,
                                     'total_used_reads' = total_used_reads)
               if (length(output_tx) == 0) {output_tx <- df.data} else {output_tx <- rbind(output_tx,df.data)}
               
            }
         } 
         message('Event: ',event,'\tRPs: ',format(n_spann_RPs,width=5),'\t',length(tx),' transcripts assembled.\tk = ',kmer_size,'\texp_cov: ',format(as.numeric(exp_cov),digits = 3),'\tcov_cutoff: ',format(as.numeric(min(3,cov_cutoff)),digits = 3))
         
      } else if (velvet_check) {
         tic()
         # Regular velvet ------
         message(paste('Assembling ',n_spann_RPs,' read pairs at event ',event,'...',sep=''))
         ## Loop over kmers ----
         for (kmer in k){ 
            kmer.tx.dir <- file.path(event.tx.dir,paste('k_',kmer,sep=''))
            if (!dir.exists(kmer.tx.dir)) {dir.create(kmer.tx.dir)}
            tx.file <- file.path(kmer.tx.dir,'transcripts.fa')
            
            
            # debug: 
            # cov_cutoff = 0; take.rc = FALSE; 
            # short.file <- event.readPairs.fa ; short.paired.file <- event.flashed.reads.fa
            
            
            # 
            #                   x = seq.df; iJct = event;  k = kmer; event.tx.dir = event.tx.dir;
            #                   file.prefix = file.prefix; cov_cutoff = cov_cutoff; take.rc = FALSE;
            #                   flash.flag = flash.flag; 
            # 
            velvetg_cov_cutoff = 0; oases_cov_cutoff = cov_cutoff
            
            tx <- velvetOases(iJct = event,
                              k = kmer,
                              event.tx.dir = event.tx.dir,
                              file.prefix = file.prefix,
                              velvetg_cov_cutoff = velvetg_cov_cutoff,
                              oases_cov_cutoff = oases_cov_cutoff,
                              flash.flag = flash.flag,
                              insert_length = insert_length)
            
            ## Remove transcripts that only use reference sequence ----
            # TODO: make all read_list objects read directly from read_tx_hash.tab
            
            
            # if(length(tx) > 0 && length(tx) < sum(grepl('ENST',contig.out$reads.df$Read))){
            #    print(paste('Event ',j,': More reference TX than RNA reads!',sep=''))
            # }
            if (length(tx) > 0) {
               # debug: tx_out_dir = kmer.tx.dir;merged = FALSE
               contig.out <- contigParser(tx_out_dir = kmer.tx.dir,merged = FALSE)
               read_tx.hash <- contig.out$read_tx.hash
               total_used_reads <- sum(rowSums(read_tx.hash) > 0)
               
               reads.df <- data.frame()
               transcripts.df <- data.frame()
               # reads.df <- contig.out$reads.df[contig.out$reads.df$nTX  > 0,]
               # transcripts.df <- contig.out$transcripts.df
               
               # loc_tx_str <- sub('Locus','Loc',sub('Transcript','TX',transcripts.df$Name))
               # if (any(!grepl(paste(loc_tx_str,collapse = '|'),names(tx)))){
               #    message('event',event,'kmer',k, 'TX - contig name mismatch?')
               #    tx <- tx[grepl(paste(loc_tx_str,collapse = '|'),names(tx))]
               # }
               
               
               tx.info <- strsplit(names(tx),'_')
               tx.names <- strsplit(names(tx),'/')
               
               
               for (i_tx in 1:length(tx)) {
                  tx_conf <- as.numeric(tx.info[[i_tx]][which(grepl('Conf',tx.info[[i_tx]]))+1])
                  # seq.df_data <- unique(seq.df[seq.df$DNA_eventId == event,df_header])
                  
                  
                  
                  n_reads = sum(read_tx.hash[,i_tx])
                  #reads_list = rownames(read_tx.hash)[read_tx.hash[,i_tx] == 1]
                  #total_used_reads <- sum(read_tx.hash)
                  
                  
                  
                  df.data <- data.frame(unique_e_seq.df,
                                        'kmer_size' = kmer,
                                        'cov_cutoff' = cov_cutoff,
                                        # 'exp_cov' = exp_cov,
                                        'TX' = toString(tx[i_tx]),
                                        'tx_num' = tx.names[[i_tx]][1],
                                        'tx_conf' = tx_conf,
                                        'tx_length' = width(tx[i_tx]),
                                        'n_reads' = n_reads,
                                        # 'seq_list' = reads_list,
                                        # 'total_event_reads' = 2*sum(seq.df$DNA_eventId == event),
                                        'total_used_reads' = total_used_reads)
                  
                  if (length(output_tx) == 0) {output_tx <- df.data} else {output_tx <- rbind(output_tx,df.data)}
               } # loop: building output_tx one tx at a time
               
            }  # check length of tx for this velvet-oases run
            
         } # looping over kmer in k
         # note: output_tx is defined BEFORE event (j)
         
         # Write assembled novel TX to file in event dir for assembly ----
         tx.ss <- DNAStringSet(x = output_tx$TX[output_tx$DNA_eventId == event ],use.names = TRUE)
         names(tx.ss) <- output_tx$tx_num[output_tx$DNA_eventId == event & output_tx$kmer_size != 'merged']
         tx.fa_file <- file.path(event.tx.dir,'raw_assembled_transcripts.fa')
         writeXStringSet(x = tx.ss,filepath = tx.fa_file) 
         
         # Merge assembled TX ----
         merged.tx.dir <- paste(event.tx.dir,'merged',sep='/')
         if (!dir.exists(merged.tx.dir)) {dir.create(merged.tx.dir)}  
         log.file <- file.path(merged.tx.dir,format(Sys.time(),'oases_%Y_%H%M%S_R.log'))
         
         velveth.cmd <- paste(velvet_path,'/velveth ',merged.tx.dir,' 27 -long ', tx.fa_file ,' >> ',log.file,sep='')
         velvetg.cmd <- paste(velvet_path,'/velvetg ',merged.tx.dir,' -read_trkg yes -conserveLong yes >> ',log.file,sep='')
         oases.cmd <- paste(oases_path, '/oases ',merged.tx.dir,' -merge yes >> ',log.file,sep='')
         
         system(velveth.cmd)
         system(velvetg.cmd)
         system(oases.cmd,ignore.stdout = TRUE)
         
         transcripts.fa_filename <- file.path(merged.tx.dir,'transcripts.fa')
         reformatFastaHeader(filename = transcripts.fa_filename,jct_no = event,k = 'merged')
         
         tx <- readDNAStringSet(filepath = paste(merged.tx.dir,'/transcripts.fa',sep = ''))
         merged_tx <- data.frame()
         # Parsing contigs of merged TX ----
         if (length(tx) > 0 ) {
            contig.out <- contigParser(tx_out_dir = merged.tx.dir,merged = TRUE)
            merged_read_tx.hash <- contig.out$merged_read_tx.hash
            rna_read_list <- rownames(merged_read_tx.hash)
            
            
            
            #message(paste('Event ',j,'; ',sum(seq.df$DNA_eventId == j),' associated reads, ',length(tx),' transcript(s)',sep = ''))
            tx.info <- strsplit(names(tx),'_')
            tx.names <- sapply(strsplit(names(tx),'/'),'[',1)
            total_used_reads <- sum(rowSums(merged_read_tx.hash) > 0)
            for (i_tx in 1:length(tx)) {
               
               n_reads <- sum(merged_read_tx.hash[,tx.names[i_tx]] == 1)
               
               
               ## format output data ----
               df.data <- data.frame(unique_e_seq.df,
                                     'kmer_size' = 'merged',
                                     'cov_cutoff' = cov_cutoff,
                                     'TX' = toString(tx[i_tx]),
                                     'tx_num' = tx.names[[i_tx]][1],#unlist(strsplit(tx.info[[i]][which(grepl('TX',tx.info[[i]]))+1],'/'))[1],
                                     'tx_conf' =  as.numeric(tx.info[[i_tx]][which(grepl('Conf',tx.info[[i_tx]]))+1]),
                                     'tx_length' = width(tx[i_tx]),
                                     'n_reads' = n_reads,
                                     # 'seq_list' = reads_list,
                                     # 'total_event_reads' = 2*sum(seq.df$DNA_eventId == event),
                                     'total_used_reads' = total_used_reads)
               
               if (length(merged_tx) == 0){merged_tx <- df.data} else {merged_tx <- rbind(merged_tx,df.data)}
            } # 
            
            
         } else {# merged tx > 0 check
            
         }
         dt <- toc(quiet=TRUE)
         # message('Assembled and merged ',nrow(merged_tx),' transcript(s) in ',format(dt$toc-dt$tic,nsmall=1),' seconds.')
         message('Event: ',event,',\tRPs: ',format(n_spann_RPs,width=5),' \t\t ',nrow(merged_tx),'merged transcripts.')
         
         # toc()
         
      } else {
         message('Event: ',event,'\tRPs: ',format(n_spann_RPs,width=5),' \t\t Velvet pre-check failed; no transcripts assembled.')
      }# use velvet optimiser?
      
   } # looping over event in event_list
   
   if (optimize_velvet){tx.df <- output_tx} else {tx.df <- merged_tx}
   
   ## Cleanup -----
   if (cleanup){
      merged.tx.dir <- paste(out.dir,'merged',sep='/')
      unlink(merged.tx.dir,recursive = TRUE)
      ##kmer_dirs <- paste(tx.dir,'/k_',k,sep='')
      for (kmer in k){
         kmer_dir <- paste(out.dir,'/k_',k,sep='')
         unlink(kmer_dir,recursive = TRUE)
      }
   }
   
   
   ## output ----
   tx_csv_file <- file.path(out.dir,paste(file.prefix,'TX_out.csv',sep='_'))
   if (file.exists(tx_csv_file)){unlink(tx_csv_file)}
   write.csv(x = tx.df,file = tx_csv_file,quote = FALSE,row.names = FALSE)
   return(tx.df)
   
}



# TODO: verify correct elimination of case where merge_overlapping_event = FALSE - some DNA_posA are not unique? MCF7 718?
loadRNAReads <- function(readsFileName,merge_overlapping_events = FALSE){
   #debug: readsFileName <- read_file
   case_dir <- strsplit(readsFileName,'/')
   seq.df <- data.frame()
   for (i_file in 1:length(readsFileName)) {
      
      tmp.df <- read.csv(file = readsFileName[i_file],sep = ',',stringsAsFactors = FALSE,fill = TRUE,header = FALSE)
      toSkip <- sum(grepl('#',tmp.df[,]))
      new.df <- read.csv(file = readsFileName[i_file],sep = ',',stringsAsFactors = FALSE,skip = toSkip)
      
      if (i_file == 1) {seq.df <- new.df} else {seq.df <- rbind(seq.df,new.df)}
   }
   colnames(seq.df) <- sub('MP','DNA',colnames(seq.df))
   #if (!exists(event_list)){event_list <- unique(seq.df$DNA_eventId)}
   
   
   
   if (nrow(seq.df) > 0) {
      
      event_list <- unique(seq.df$DNA_eventId)
      
      if (merge_overlapping_events){
         
         read_ol_mat <- matrix(0,nrow = length(event_list),ncol = length(event_list))
         event_rename_list <- seq.df$DNA_eventId
         
         #b_ol_mat <- matrix(0,nrow = length(event_list),ncol = length(event_list))
         for (i_event in 1:length(event_list)-1){
            for (j_event in i_event:length(event_list)){
               
               
               i_event_reads <- seq.df$seqID[seq.df$DNA_eventId == event_list[i_event]]
               j_event_reads <- seq.df$seqID[seq.df$DNA_eventId == event_list[j_event]]
               shared_reads <- sum(j_event_reads %in% i_event_reads)
               
               if (i_event != j_event && all(j_event_reads %in% i_event_reads)){
                  #event_subset_list <- 
                  event_rename_list[event_rename_list == event_list[j_event]] <- rep(event_list[i_event],
                                                                                     sum(event_rename_list == event_list[j_event]))
                  # TODO: message / log indicating which events were merged and why
               }
               # foo <- findOverlaps(seq.A.gr[i_event],seq.A.gr[j_event])
               #read_ol_mat[i_event,j_event] <- shared_reads/length(i_event_reads)
            }
            
         }
         #table(seq.df$DNA_eventId)
         #table(event_rename_list)
         seq.df$DNA_eventId = event_rename_list
         event_list <- unique(seq.df$DNA_eventId)
         
         
         
         #seq.df <- unique(seq.df)
         
         
         for (j in 1:length(event_list)) {
            jct <- seq.df[seq.df$DNA_eventId == event_list[j],]
            seq.df$RNA_RPs[which(seq.df$DNA_eventId == event_list[j])] <- rep(length(unique(jct$seqID)),nrow(jct))
            #seq.df$gene_prediction[which(seq.df$DNA_eventId == event_list[j])] <- rep(length(unique(jct$seqID)),nrow(jct))
            
            # Warning: commenting out the below may cause back compatibility issues if running on older (ie pre 2020) pipeline output
            #if (any(is.na(jct$gene_prediction))) {
            #   if (any(is.na(jct$gene_prediction) & sum(!is.na(jct$DNA_gene_prediction)) > 0)) {
            #      seq.df$gene_prediction <-  seq.df$DNA_gene_prediction
            #   }
            #   #jct$gene_prediction[is.na(jct$gene_prediction)] <- rep(unique(jct$gene_prediction[!is.na(jct$gene_prediction)]),sum(is.na(jct$gene_prediction)))
            #} 
            
            for (h in c('geneA','geneB','gene_prediction')) {
               old_gene_list <- unique(seq.df[seq.df$DNA_eventId == event_list[j],h])
               if (length(old_gene_list) > 1) {
                  gene_list <- sort(unique(unlist(strsplit(old_gene_list,';'))))
                  new_gene_list <- paste(gene_list[!gene_list %in% 'no_gene'],collapse = ';')
                  seq.df[seq.df$DNA_eventId == event_list[j],h] <- rep(new_gene_list,sum(seq.df$DNA_eventId == event_list[j]))
               }
               
            }
            if (nrow(unique(jct[,c('posAmin','posAmax')])) > 1) {
               new_range <- c(min(jct[,'posAmin']),max(jct[,'posAmax']))
               nRows_to_rep <- sum(seq.df$DNA_eventId == event_list[j])
               seq.df[seq.df$DNA_eventId == event_list[j],'posAmin'] <- rep(new_range[1],nRows_to_rep)
               seq.df[seq.df$DNA_eventId == event_list[j],'posAmax'] <- rep(new_range[2],nRows_to_rep)
               #seq.df[seq.df$DNA_eventId == event_list[j],c('posAmin','posAmax')] <- do.call('rbind',replicate(nRows_to_rep,new_range,simplify = FALSE))
            }
            if (nrow(unique(jct[,c('posBmin','posBmax')])) > 1) {
               new_range <- c(min(jct[,'posBmin']),max(jct[,'posBmax']))
               nRows_to_rep <- sum(seq.df$DNA_eventId == event_list[j])
               seq.df[seq.df$DNA_eventId == event_list[j],'posBmin'] <- rep(new_range[1],nRows_to_rep)
               seq.df[seq.df$DNA_eventId == event_list[j],'posBmax'] <- rep(new_range[2],nRows_to_rep)
               #seq.df[seq.df$DNA_eventId == event_list[j],c('posAmin','posAmax')] <- do.call('rbind',replicate(nRows_to_rep,new_range,simplify = FALSE))
            }
            gene_list <- unique(c(
               unlist(strsplit(seq.df$geneA[seq.df$DNA_eventId == event_list[j]],';')),
               unlist(strsplit(seq.df$geneB[seq.df$DNA_eventId == event_list[j]],';'))
            ))
            seq.df$eventUniqueGenes[seq.df$DNA_eventId == event_list[j]] <- 
               sum(gene_list != 'no_gene')
            
         }
      } else {
         for (j in 1:length(event_list)){
            jct <- seq.df[seq.df$DNA_eventId == event_list[j],]
            seq.df$RNA_RPs[which(seq.df$DNA_eventId == event_list[j])] <- rep(length(unique(jct$seqID)),nrow(jct))
            gene_A_list <- unique(
               unlist(strsplit(seq.df$geneA[seq.df$DNA_eventId == event_list[j]],';')))
            gene_B_list <- unique(
               unlist(strsplit(seq.df$geneB[seq.df$DNA_eventId == event_list[j]],';')))
            gene_list <- unique(c(gene_A_list,gene_B_list))
            
            seq.df$geneA[which(seq.df$DNA_eventId == event_list[j])] <- paste(gene_A_list,collapse=';')
            seq.df$geneB[which(seq.df$DNA_eventId == event_list[j])] <- paste(gene_B_list,collapse=';')
            
            if (length(unique(jct$DNA_posA))>1){
               seq.df$DNA_posA[which(seq.df$DNA_eventId == event_list[j])] <- unique(jct$DNA_posA[jct$DNA_RPs == max(jct$DNA_RPs )])
               seq.df$DNA_posB[which(seq.df$DNA_eventId == event_list[j])] <- unique(jct$DNA_posB[jct$DNA_RPs == max(jct$DNA_RPs )])
               seq.df$DNA_RPs[which(seq.df$DNA_eventId == event_list[j])] <- max(jct$DNA_RPs )
            }
            
            seq.df$eventUniqueGenes[seq.df$DNA_eventId == event_list[j]] <- 
               sum(gene_list != 'no_gene')
            
         }
      }
      
      seq.df$chrA <- as.character(seq.df$chrA)
      seq.df$chrA[seq.df$chrA == '23'] <- 'X'
      seq.df$chrA[seq.df$chrA == '24'] <- 'Y'
      seq.df$chrA[seq.df$chrA == '25'] <- 'M'
      
      seq.df$chrB[seq.df$chrB == '23'] <- 'X'
      seq.df$chrB[seq.df$chrB == '24'] <- 'Y'
      seq.df$chrB[seq.df$chrB == '25'] <- 'M'
      
      
   }
   return(seq.df)
}



trimTrailingPolyA <- function(x,nA){
   polyA_regexp = paste('A{',nA,',}',sep='')
   m <- gregexpr(polyA_regexp,x,perl = TRUE)[[1]]
   m.pos <- m[length(m)]
   m.len <- attr(m,'match.length')[length(m)]
   if (m.pos + m.len - 1 == nchar(x)){
      y <- substr(x, start = 1,stop = m.pos - 1)
   } else {y = x}
   return(y)
}

trimLeadingPolyT <- function(x,nT){
   polyT_regexp = paste('T{',nT,',}',sep='')
   m <- gregexpr(polyT_regexp,x,perl = TRUE)[[1]]
   m.pos <- m[1]
   m.len <- attr(m,'match.length')[1]
   if (m.pos == 1){
      y <- substr(x, start = m.pos + m.len,stop = nchar(x))
   } else {y = x}
   return(y)
}

flash <- function(seq.df,event,event.tx.dir,out.prefix){
   interlaced_tsv_filename <- file.path(event.tx.dir,'interlaced_RPs.tsv')
   if (file.exists(interlaced_tsv_filename)){
      unlink(interlaced_tsv_filename)
   }
   seq_idx <- which(seq.df$DNA_eventId == event)
   n = length(seq_idx)
   max_overlap_param <- min(unique(nchar(seq.df[seq_idx,'readA'])))
   interlaced.df <- data.frame('tag' = paste(event,1:n,sep='_'),
                               'readA' = seq.df[seq_idx,'readA'],
                               'qualA' = character(n),
                               'readB' = seq.df[seq_idx,'readB'],
                               'qualB' = character(n))
   N = unique(nchar(c(interlaced.df$readA,interlaced.df$readB)))
   if (length(N) == 1){
      interlaced.df$qualA <- paste(rep('.',N),collapse = '')
      interlaced.df$qualB <- paste(rep('.',N),collapse = '')
   } else {
      for (i_seq in 1:nrow(interlaced.df)){
         interlaced.df[i_seq,'qualA'] <- paste(rep('.',nchar(interlaced.df[i_seq,'readA'])),collapse='')
         interlaced.df[i_seq,'qualB'] <- paste(rep('.',nchar(interlaced.df[i_seq,'readB'])),collapse='')
      }
   }
   write_tsv(x = interlaced.df,file = interlaced_tsv_filename,
             append = FALSE,col_names = FALSE,quote = 'none')
   
   
   cmd <- paste(flash_path,'-q -x 0 -t 10 -Ti -To -O -o "flash.out" -M',max_overlap_param,'-d',event.tx.dir,interlaced_tsv_filename)
   suppressMessages(system(cmd))
   tab_outfile <- file.path(event.tx.dir,'flash.out.readsAndPairs.tab')
   new_reads <- strsplit(readLines(tab_outfile),'\t')
   nr_len <- unlist(lapply(X = new_reads,FUN = 'length'))
   RPs.ss <- DNAStringSet(x = sapply(new_reads[nr_len == 5],'[',c(2,4)))
   if (length(RPs.ss) > 0){names(RPs.ss) <- paste(event,c('A','B'),sub(paste(event,'_',sep=''),'',sapply(new_reads[nr_len == 5],'[',c(1,1))),sep='_')}
   
   ol_RPs.ss <- DNAStringSet(x = sapply(new_reads[nr_len == 3],'[',2))
   if (length(ol_RPs.ss) > 0){names(ol_RPs.ss) <- paste(event,'AB',sub(paste(event,'_',sep=''),'',sapply(new_reads[nr_len == 3],'[',1)),sep='_')}
   
   #spanningReads.ss <- c(RPs.ss,ol_RPs.ss)
   
   unlink(file.path(event.tx.dir,grep('flash.out',list.files(event.tx.dir),value = TRUE)))
   unlink(interlaced_tsv_filename)
   flash.out <- list('pairs' = RPs.ss,'combined' = ol_RPs.ss)
   return(flash.out)
}


writeJunctionFa <- function(reads,pairedFastaFileRoot,write_paired) {
   if (missing(write_paired)) {write_paired <- FALSE}
   if (grepl('.fa',pairedFastaFileRoot)){
      pairedFastaFileRoot <- sub('.fa','',pairedFastaFileRoot)
   }
   readsA <- reads[grepl('A',names(reads))]
   readsB <- reads[grepl('B',names(reads))]
   reads.il <- DNAStringSet()
   for (i in 1:length(readsA)) {
      reads.il <- c(reads.il,readsA[i],readsB[i])
   }
   if (write_paired) {
      writeXStringSet(readsA, paste(pairedFastaFileRoot,'_A.fa',sep = ''), append = FALSE, format = "fasta")
      writeXStringSet(readsB, paste(pairedFastaFileRoot,'_B.fa',sep = ''), append = FALSE, format = "fasta")
   } else {
      writeXStringSet(reads.il, paste(pairedFastaFileRoot,'.fa',sep = ''), append = FALSE, format = "fasta")
   }
}



reformatFastaHeader <- function(filename,jct_no,k,use_kmer_in_tx_name = TRUE){
   tx.file <- filename
   tx.file.tmp.1 <- paste(tx.file,'tmp.1',sep='.')
   tx.file.tmp.2 <- paste(tx.file,'tmp.2',sep='.')
   tx.file.tmp.3 <- paste(tx.file,'tmp.3',sep='.')
   tx.file.tmp.4 <- paste(tx.file,'tmp.4',sep='.')
   tx.file.tmp.5 <- paste(tx.file,'tmp.5',sep='.')
   if (missing(k) | !use_kmer_in_tx_name){   
      sed.cmd.1 <- paste('sed "s/>/>J_',jct_no,'_/g" ',tx.file, ' > ', tx.file.tmp.1,sep='')
   }   else {
      sed.cmd.1 <- paste('sed "s/>/>J_',jct_no,'_k_',k,'_/g" ',tx.file, ' > ', tx.file.tmp.1,sep='')
      
   }
   system(sed.cmd.1)
   Sys.sleep(rw.pause.time)
   sed.cmd.2 <- paste('sed "s/Locus/Loc/g" ',tx.file.tmp.1,' > ', tx.file.tmp.2,sep='')
   system(sed.cmd.2)
   Sys.sleep(rw.pause.time)
   sed.cmd.3 <- paste('sed "s/Transcript/TX','/g" ',tx.file.tmp.2,' > ', tx.file.tmp.3,sep='')
   system(sed.cmd.3)
   Sys.sleep(rw.pause.time)
   sed.cmd.4 <- paste('sed "s/Confidence/Conf','/g" ',tx.file.tmp.3,' > ', tx.file.tmp.4,sep='')
   system(sed.cmd.4)
   Sys.sleep(rw.pause.time)
   sed.cmd.5 <- paste('sed "s/Length/L','/g" ',tx.file.tmp.4,' > ', tx.file.tmp.5,sep='')
   system(sed.cmd.5)
   Sys.sleep(rw.pause.time)
   mv.cmd <- paste('cp',tx.file.tmp.5,tx.file)
   system(mv.cmd)
   rm.cmd <- paste('rm',tx.file.tmp.1, tx.file.tmp.2, tx.file.tmp.3, tx.file.tmp.4, tx.file.tmp.5)
   system(rm.cmd)
}

spanningRPsPreProcessing <- function(seq.df,event,event.tx.dir,file.prefix,
                                     flash.threshold = 100,
                                     remove.polyA = FALSE,
                                     prehash_reads = FALSE){
   e_seq.df <- seq.df[seq.df$DNA_eventId == event,]
   event.readPairs.fa <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'_pairs.fa',sep=''))
   event.flashed.reads.fa <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'_flash.fa',sep=''))
   
   
   #out.file <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'.fa',sep=''))
   if (nrow(e_seq.df) > flash.threshold){
      flash.out <- flash(e_seq.df,event,event.tx.dir,out.prefix = out.prefix)
      RP.ss <- flash.out[['pairs']]
      cRP.ss <- flash.out[['combined']]
      
      readsPP.out <- list('n_combined' = length(cRP.ss),'n_pairs' = length(RP.ss)/2,
                          'combined_file' = event.flashed.reads.fa, 'paired_file' = event.readPairs.fa)
      if (readsPP.out[['n_combined']] > 0){writeXStringSet(x = cRP.ss,file = event.flashed.reads.fa)}
      if (readsPP.out[['n_pairs']] > 0){writeXStringSet(x = RP.ss,file = event.readPairs.fa)}
   } else {
      RP.ss <- DNAStringSet(x = c(e_seq.df$readA,e_seq.df$readB))
      names(RP.ss) <- c(paste(event,'A',1:nrow(e_seq.df),sep='_'),paste(event,'B',1:nrow(e_seq.df),sep='_'))
      readsPP.out <- list('n_combined' = 0,'n_pairs' = length(RP.ss)/2,
                          'combined_file' = event.flashed.reads.fa, 'paired_file' = event.readPairs.fa)
      writeJunctionFa(reads = RP.ss,pairedFastaFileRoot = event.readPairs.fa)
      if (file.exists(event.flashed.reads.fa)){unlink(event.flashed.reads.fa)}
   }
   # remove.polyA = FALSE
   # if (remove.polyA){
   #    nA = 7
   #    polyA_idx <- which(endsWith(as.character(reads.ss),paste(rep('A',nA),collapse='')))
   #    reads.ss[polyA_idx] <- sapply(as.character(reads.ss[polyA_idx]),'trimTrailingPolyA',nA)
   #    
   #    polyT_idx <- which(startsWith(as.character(reads.ss),paste(rep('T',nA),collapse='')))
   #    reads.ss[polyT_idx] <- sapply(as.character(reads.ss[polyT_idx]),'trimLeadingPolyT',nA)
   # }
   
   # readsPP.out <- list('n_combined' = length(cRP.ss),'n_pairs' = length(RP.ss)/2,
   #                   'combined_file' = event.flashed.reads.fa, 'paired_file' = event.readPairs.fa)
   # if (readsPP.out[['n_combined']] > 0){writeXStringSet(x = cRP.ss,file = event.flashed.reads.fa)}
   # if (readsPP.out[['n_pairs']] > 0){writeXStringSet(x = cRP.ss,file = event.flashed.reads.fa)}
   # 
   # if (readsPP.out[['n_pairs']] > 0){writeJunctionFa(reads = RP.ss,pairedFastaFileRoot = event.readPairs.fa)}
   
   # Recall: 
   # event.readPairs.fa <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'_pairs.fa',sep=''))
   # event.flashed.reads.fa <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'_flash.fa',sep=''))
   # 
   if (file.exists(event.readPairs.fa) & file.size(event.readPairs.fa) > 0){
      short.paired_param <- paste('-shortPaired',event.readPairs.fa)
   } else {short.paired_param <- ''}
   if (readsPP.out$n_combined > 0 & file.exists(event.flashed.reads.fa) & file.size(event.flashed.reads.fa) > 0){
      short.combined_param <- paste('-short',event.flashed.reads.fa)
   } else {short.combined_param <- ''}
   
   if (prehash_reads){
      velveth.cmd <- paste(paste(velvet_path,'velveth',sep='/'), event.tx.dir, 21,'-fasta',
                           short.paired_param,short.combined_param,'-noHash')
      system(velveth.cmd, ignore.stdout = TRUE)
   }
   
   return(readsPP.out)
}

velvetOptimiser <- function(event.tx.dir,
                            file.prefix = NA,
                            insert_length = NA,
                            flash.flag = TRUE,
                            oases_cov_cutoff = 0){
   event = unlist(strsplit(grep('event',unlist(strsplit(event.tx.dir,'/')),value=TRUE),'_'))[2]
   
   short.paired.file <- file.path(event.tx.dir,grep('pairs.fa',list.files(event.tx.dir),value=TRUE))
   short.file = file.path(event.tx.dir,grep('flash.fa',list.files(event.tx.dir),value=TRUE))
   
   if (length(short.paired.file) > 0 &&
       file.exists(short.paired.file) & file.size(short.paired.file) > 0){
      short.paired_param <- paste('-shortPaired',short.paired.file)
   } else {short.paired_param <- ''}
   if (flash.flag & length(short.file) > 0 ){
      short.combined_param <- paste('-short',short.file)
   } else {short.combined_param <- ''}
   
   
   
   if (oases_cov_cutoff != 'auto' ){
      oases_cov_cutoff = min(3,oases_cov_cutoff)
      cov_cutoff_param <- paste('-cov_cutoff',format(oases_cov_cutoff,digits = 3))
      
   } else {cov_cutoff_param <- ''}
   
   if (!is.na(insert_length)){
      insert_length_param <- paste('-ins_length',insert_length)
   } else {insert_length_param = ''}
   
   #cov_cutoff <- 0
   scaffolding = FALSE
   if (scaffolding){
      scaffolding_param <- ''
   } else {scaffolding_param = '-scaffolding no'}
   
   
   
   assembly_dir <- file.path(event.tx.dir,'assembly/')
   assembly_prefix <- 'vo'
   if (!dir.exists(assembly_dir)){dir.create(assembly_dir)}
   
   vo_outdir <- file.path(event.tx.dir,'vo')
   if (dir.exists(vo_outdir)){unlink(vo_outdir,recursive = TRUE)}
   # log.file <- file.path(vo_outdir,format(Sys.time(),'VO_oases_%Y_%H%M%S.log'))
   
   setwd(file.path(event.tx.dir,'assembly'))
   # setwd(assembly_dir)
   # system(paste('cd',assembly_dir))
   
   vo.cmd <- paste(perl_path, vo_path,
                   '-s 11 -e 31 -x 4',
                   '-p',assembly_prefix,
                   '-f \'',paste(short.paired_param,short.combined_param),'\'',
                   '-d',vo_outdir,'-a','-c tbp',
                   '-o','\'',paste(insert_length_param,scaffolding_param,'-shortMatePaired yes'),'\'')
   
   system(vo.cmd,ignore.stderr = TRUE,ignore.stdout = TRUE)
   
   # Side-stepping VO's finnicky programming ----
   if (!dir.exists(vo_outdir)){
      vo_log <- readLines(file.path(event.tx.dir,'assembly','vo_logfile.txt'))
      assembly_id <- unlist(regmatches(x = grep('Best assembly by assembly',vo_log,value=TRUE),m = gregexpr('(?<=Best assembly by assembly score - assembly id: )\\d',grep('Best assembly by assembly',vo_log,value=TRUE),perl=TRUE)))
      assembly_block <- vo_log[tail(grep(paste('Assembly id: ',assembly_id,sep=''),vo_log),1): 
                                  (tail(grep('\\*',vo_log,perl=TRUE),1)-1)]
      assembly_dir_line <- grep('Assembly directory:',assembly_block,value=TRUE)
      assembly_kmer_line <- grep('Velvet hash value: ',assembly_block,value=TRUE)
      assembly_kmer <- unlist(regmatches(m = gregexpr('(?<=Velvet hash value: )\\d{1,}',assembly_kmer_line,perl=TRUE),x = assembly_kmer_line))
      
      dir_for_oases <- unlist(regmatches(x = assembly_dir_line,m = gregexpr('(?<=Assembly directory: )\\S+',assembly_dir_line,perl=TRUE)))
      
      data <- read.table(file.path(dir_for_oases,'stats.txt'),header = TRUE)
      exp_cov <- weighted.mean(data$short1_cov, data$lgth)
      velvetg_cov_cutoff = exp_cov / 2
      
      tx <- velvetOases(iJct = event,kmer = assembly_kmer,
                        file.prefix = file.prefix,
                        event.tx.dir = event.tx.dir,
                        velvetg_cov_cutoff = velvetg_cov_cutoff,
                        oases_cov_cutoff = oases_cov_cutoff,
                        insert_length = insert_length,
                        exp_cov = exp_cov,
                        flash.flag = flash.flag, 
                        use_kmer_in_tx_name = FALSE)
      message('VO fail :( using vanilla velvetOases')
   } else {
      dir_for_oases <- vo_outdir
      
      oases.cmd <- paste(paste(oases_path,'oases',sep='/'), dir_for_oases,cov_cutoff_param)
      system(oases.cmd,ignore.stderr = TRUE,ignore.stdout = TRUE)
      Sys.sleep(rw.pause.time)
      tx.file <- file.path(vo_outdir,'transcripts.fa')
      
      if (file.exists(tx.file)){
         reformatFastaHeader(filename = tx.file,jct_no = event)
         tx <- readDNAStringSet(filepath = tx.file)
      } else {tx = DNAStringSet()}
   }
   
   
   #    message('No TX assembled.')
   #    tx = DNAStringSet()
   # }
   # if (dir.exists(old_wd)){setwd(old_wd)}
   setwd(event.tx.dir)
   unlink(assembly_dir,recursive = TRUE)
   return(tx)
}

velvetOases <- function(iJct,kmer,event.tx.dir,file.prefix,
                        flash.flag = FALSE,
                        take.rc = FALSE,
                        oases_cov_cutoff = 3,
                        velvetg_cov_cutoff = 0,
                        exp_cov = 'auto',
                        use_kmer_in_tx_name = FALSE,
                        insert_length = NA) {
   
   
   
   #if (missing(verbose)) {verbose = FALSE}
   if (missing(file.prefix)) {file.prefix <- ''}
   # if (missing(kmer)) {kmer <- 13}
   if (missing(take.rc)){take.rc = FALSE}
   kmer_dir <- file.path(event.tx.dir,paste('k',kmer,sep='_'))
   if (!dir.exists(kmer_dir)) {dir.create(kmer_dir)}
   log.file <- paste(kmer_dir,format(Sys.time(),'oases_run_%Y_%H%M%S_R.log'),sep = '/')
   #if (file.exists(log.file)) {system(paste('rm',log.file))}
   # write fasta file
   sequences.file <- file.path(event.tx.dir,'Sequences')
   if (file.exists(sequences.file)){
      seq_cp.cmd <- paste('cp',sequences.file,kmer_dir)
      system(seq_cp.cmd)
      reuse_sequences_param <- '-reuseSequences'
   } else {
      reuse_sequences_param <- ''
   }
   tx.file <- file.path(kmer_dir,'transcripts.fa')
   if (file.exists(tx.file)){unlink(tx.file)}
   
   velvet.out.dir <- file.path(event.tx.dir,paste('k',kmer,sep='_'),paste(file.prefix,'_event_',iJct,'.fa',sep = ''))
   # tx.file <- paste(kmer_dir,'/transcripts.fa',sep = '')
   
   
   # velveth -----
   
   
   
   short.paired.file <- file.path(event.tx.dir,paste(file.prefix,'_event_',iJct,'_pairs.fa',sep=''))
   short.file = file.path(event.tx.dir,paste(file.prefix,'_event_',iJct,'_flash.fa',sep=''))
   
   if (file.exists(short.paired.file) & file.size(short.paired.file) > 0){
      short.paired_param <- paste('-shortPaired',short.paired.file)
   } else {short.paired_param <- ''}
   if (flash.flag & file.exists(short.file) & file.size(short.file) > 0){
      short.combined_param <- paste('-short',short.file)
   } else {short.combined_param <- ''}
   
   
   if (reuse_sequences_param == ''){
      velveth.cmd <- paste(paste(velvet_path,'velveth',sep = '/'),kmer_dir,kmer,'-fasta',
                           short.paired_param,short.combined_param,
                           '>>',log.file)
   } else {
      velveth.cmd <- paste(paste(velvet_path,'velveth',sep = '/'),kmer_dir,kmer,
                           '-reuse_Sequences >>',log.file)
   }
   
   
   
   # system(velveth.cmd)
   
   # velvetg ----
   # use_exp_cov <- FALSE
   if (!missing(exp_cov)){
      # exp_cov <- 2;
      #exp_cov <- unique(x$RNA_RPs)/3
      exp_cov_param <- paste('-exp_cov',exp_cov)
   } else {exp_cov_param <- ''}
   
   if (!is.na(insert_length)){
      insert_length_param <- paste('-ins_length',insert_length)
   } else {insert_length_param = ''}
   
   if (velvetg_cov_cutoff != 'auto'){
      velvetg_cov_cutoff_param <- paste('-cov_cutoff',velvetg_cov_cutoff)
   } else {velvetg_cov_cutoff_param <- ''}
   
   if (oases_cov_cutoff != 'auto'){
      oases_cov_cutoff_param <- paste('-cov_cutoff',oases_cov_cutoff)
   } else {oases_cov_cutoff_param <- ''}
   
   #cov_cutoff <- 0
   scaffolding = FALSE
   if (scaffolding){
      scaffolding_param <- ''
   } else {scaffolding_param = '-scaffolding no'}
   
   velvetg.cmd <- paste(paste(velvet_path,'velvetg',sep = '/'),
                        kmer_dir,
                        insert_length_param,
                        exp_cov_param,
                        velvetg_cov_cutoff_param,
                        scaffolding_param,
                        '-read_trkg yes -shortMatePaired yes','>>',log.file)
   
   system(velveth.cmd)
   system(velvetg.cmd)
   
   oases.cmd <- paste(paste(oases_path,'oases',sep='/'), kmer_dir,insert_length_param,oases_cov_cutoff_param,'-unused_reads yes >>',log.file)
   system(oases.cmd,ignore.stderr = TRUE,ignore.stdout = TRUE)
   
   
   
   
   # oases ----
   
   
   Sys.sleep(rw.pause.time)
   
   attempts = 0
   
   while ((!file.exists(tx.file) | file.mtime(tx.file) < Sys.time() - 5) & attempts < 10){
      
      system(oases.cmd,ignore.stderr = TRUE)
      attempts = attempts + 1
      Sys.sleep(rw.pause.time)
   }
   if (file.exists(tx.file)){
      reformatFastaHeader(filename = tx.file,jct_no = iJct,k = kmer,use_kmer_in_tx_name = use_kmer_in_tx_name)
      tx <- readDNAStringSet(filepath = tx.file)
      
   } else {
      message('Oases failed!!!')
      tx = DNAStringSet()
   }
   
   # tx.file <- file.path(kmer_dir,'transcripts.fa')
   
   #unused_reads <- readDNAStringSet(filepath = paste(pathToFaDir,'/UnusedReads.fa',sep = ''))
   return(tx)
} 


contigParser <- function(tx_out_dir,merged){
   if (missing(merged)){merged <- FALSE}
   contig_dir <- file.path(tx_out_dir,'contigs')
   if (!dir.exists(contig_dir) & merged){
      dir.create(contig_dir)
   } else if (dir.exists(contig_dir)) {
      unlink(contig_dir,recursive = TRUE)
      # dir.create(contig_dir)
   }
   kmer_dir <- grep('k|merged',unlist(strsplit(tx_out_dir,'/')),value=TRUE)
   # Process LastGraph file -----
   lastGraph.file <- file.path(tx_out_dir,'LastGraph')
   lastGraph <- sapply(unlist(strsplit(read_file(file = lastGraph.file),'\n')),FUN = 'strsplit','\t')
   
   lastGraph.nodes <- data.frame(x = t(matrix(unlist(lastGraph[grep('NODE',lastGraph)]),nrow = 7))[,2:6])
   
   if (ncol(lastGraph.nodes) == 1){
      lastGraph.nodes <- t(lastGraph.nodes)
      colnames(lastGraph.nodes) <- c('nodeId','cov_short1','zero_cov_short1','cov_short2','zero_cov_short2')
   } else {
      colnames(lastGraph.nodes) <- c('nodeId','cov_short1','zero_cov_short1','cov_short2','zero_cov_short2')
   }
   # NODE $NODE_ID $COV_SHORT1 $O_COV_SHORT1 $COV_SHORT2 $O_COV_SHORT2
   
   
   if (!merged){
      lastGraph_reads.idx <- grep('NR',lastGraph)
      if (length(lastGraph_reads.idx) == 1 ){
         lastGraph.reads <- data.frame(t(unlist(lastGraph[lastGraph_reads.idx])[2:3]))
         colnames(lastGraph.reads) <- c('nodeId','n_reads')
         lastGraph.reads$lg_idx = lastGraph_reads.idx
      }else {
         lastGraph.reads <- data.frame(t(matrix(unlist(lastGraph[grep('NR',lastGraph)]),nrow = 3))[,2:3],
                                       'lg_idx' = grep('NR',lastGraph))
         colnames(lastGraph.reads)[c(1,2)] <- c('nodeId','n_reads')
      }
      
   } 
   
   # Define structures ----
   seq_names_raw <- names(read.fasta(paste(tx_out_dir,'/Sequences',sep='')))
   seq_names <- regmatches(seq_names_raw,regexpr('[.\\/\\w]{4,}(?=\\t)',seq_names_raw,perl=TRUE))
   seq_ids <- sapply(strsplit(seq_names_raw,'\t'),FUN = '[',2)
   
   reads.df <- data.frame('Read'=seq_names,
                          'readId' = seq_ids,
                          'Contig' = character(length = length(seq_names)), 
                          'nContig' = numeric(length = length(seq_names)),
                          'Transcript' = character(length = length(seq_names)),
                          'nTX' = numeric(length = length(seq_names)),
                          stringsAsFactors = FALSE)
   
   contig_ordering_file <- file.path(tx_out_dir,'contig-ordering.txt')
   contigs.fa <- read.fasta(contig_ordering_file,as.string = TRUE)
   
   nodes.fa <- contigs.fa[grepl('Node',names(contigs.fa))]
   nodes.df <- data.frame('Locus' = regmatches(names(nodes.fa),regexpr('(?<=Locus_)\\d{1,3}(?=_Node)',names(nodes.fa),perl=TRUE)),
                          'Node' = as.integer(regmatches(names(nodes.fa),regexpr('(?<=Node_)\\d{1,3}',names(nodes.fa),perl=TRUE))),
                          'Length' = nchar(as.character(nodes.fa)),
                          'nReads' = numeric(length = length(nodes.fa)),
                          seq = as.character(nodes.fa),
                          'LG_reads' = character(length=length(nodes.fa)),
                          'LG_reads_twin' = character(length=length(nodes.fa)),
                          'LG_nReads' = numeric(length = length(nodes.fa)),
                          stringsAsFactors = FALSE)
   
   transcripts.fa <- contigs.fa[grepl('Transcript',names(contigs.fa))]
   if (grepl('merged',tx_out_dir)){
      kmer_num = 'merged'
      kmer_dir = 'merged'
   } else {
      kmer_dir <- grep('k',unlist(strsplit(tx_out_dir,'/')),value=TRUE)
      kmer_num = unlist(strsplit(kmer_dir,'_'))[2]
   }
   transcripts_names_reformatted <- paste(sub('event','J',grep('event',unlist(strsplit(tx_out_dir,'/')),value=TRUE)),'k',kmer_num,
                                          sub('Transcript','TX',sub('Locus','Loc',regmatches(names(transcripts.fa),regexpr('\\w+(?=\\/)',names(transcripts.fa),perl=TRUE)))),sep='_')
   transcripts.df <- data.frame('Name' = transcripts_names_reformatted,
                                'Locus' = regmatches(names(transcripts.fa),regexpr('(?<=Locus\\_)\\d{1,3}(?=\\_Transcript)',names(transcripts.fa),perl=TRUE)),
                                'Transcript' = regmatches(names(transcripts.fa),regexpr('(?<=Transcript\\_)\\d{1,3}(?=\\/)',names(transcripts.fa),perl=TRUE)),
                                'Confidence' = as.numeric(regmatches(names(transcripts.fa),regexpr('(?<=Confidence\\_)[\\d\\.]+(?=\\_Length)',names(transcripts.fa),perl=TRUE))),
                                'Length' = regmatches(names(transcripts.fa),regexpr('(?<=Length\\_)\\d+',names(transcripts.fa),perl=TRUE)),
                                'Composition' = as.character(transcripts.fa),
                                'nReads' = numeric(length = length(transcripts.fa)),
                                stringsAsFactors = FALSE)
   
   # Prepare hash tables ----
   
   read_node_hash.mat <- matrix( 0,nrow = length(seq_names),ncol = nrow(nodes.df))
   rownames(read_node_hash.mat) <- seq_names; colnames(read_node_hash.mat) <- nodes.df$Node
   
   tx_node_hash.mat <- matrix(0, ncol = nrow(nodes.df),nrow = nrow(transcripts.df))
   rownames(tx_node_hash.mat) <- transcripts.df$Name; colnames(tx_node_hash.mat) <- nodes.df$Node
   
   read_tx_hash.mat <- matrix(0, ncol = nrow(transcripts.df),nrow = length(seq_names))
   rownames(read_tx_hash.mat) <- seq_names; colnames(read_tx_hash.mat) <- transcripts.df$Name
   
   extractContigScript <- file.path(velvet_path,'/contrib/extractContigReads/extractContigReads.pl')
   #print(paste('Running extractContigReads.pl over',length(nodes.fa),'nodes...'))
   # Loop over nodes ----
   for (iNode in 1:nrow(nodes.df)){ # fill in Reads, nReads, and seq fields of nodes.df
      
      # extract from LastGraph
      if (merged) {
         # extract contig reads to fasta
         node_name <- nodes.df$Node[iNode]
         seq_fa_filename <- paste(contig_dir,'/contig_',node_name,'_seq.fa',sep='')
         
         extract.cmd <- paste('perl',extractContigScript,node_name,tx_out_dir,'>>',seq_fa_filename)
         system(extract.cmd)
         
         Sys.sleep(0.005)
         
         contig_reads.fa <- readDNAStringSet(seq_fa_filename)
         
         read_list <- unique(names(contig_reads.fa))
         nodes.df[iNode,'nReads'] <- length(read_list)
         if (length(read_list) > 100){
            nodes.df[iNode,'Reads'] <- paste(paste(read_list[1:100],collapse=' '),'...',sep='')
         } else {
            nodes.df[iNode,'Reads'] <- paste(read_list,collapse=' ')
         }
         
         read_node_hash.mat[read_list,iNode] <- 1
         
      } else {
         i_lg_idx <- lastGraph.reads[lastGraph.reads$nodeId == nodes.df[iNode,'Node'],'lg_idx']
         n_read_ld <- as.numeric(lastGraph.reads[lastGraph.reads$nodeId == nodes.df[iNode,'Node'],'n_reads'])
         reads_from_LG <- sapply(lastGraph[(i_lg_idx+1):(i_lg_idx+n_read_ld)],FUN = '[',1)
         
         seq_names_from_LG <- reads.df[reads.df$readId %in% reads_from_LG,'Read']
         
         
         if (paste('-',nodes.df[iNode,'Node'],sep='') %in% lastGraph.reads$nodeId){
            i_twin_lg_idx <- lastGraph.reads[lastGraph.reads$nodeId == paste('-',nodes.df[iNode,'Node'],sep=''),'lg_idx']
            n_twin_read_ld <- as.numeric(lastGraph.reads[lastGraph.reads$nodeId == paste('-',nodes.df[iNode,'Node'],sep=''),'n_reads'])
            twin_reads_from_LG <- sapply(lastGraph[(i_twin_lg_idx+1):(i_twin_lg_idx+n_twin_read_ld)],FUN = '[',1)
            
            twin_seq_names_from_LG <- reads.df[reads.df$readId %in% twin_reads_from_LG,'Read']
            
            all_LG_seq_names <- unique(c(seq_names_from_LG,twin_seq_names_from_LG))
         } else {
            all_LG_seq_names <- seq_names_from_LG
            n_twin_read_ld <- 0
         }
         
         read_node_hash.mat[all_LG_seq_names,iNode] <- 1
         nodes.df[iNode,'nReads'] <- length(all_LG_seq_names)
         if (length(all_LG_seq_names) < 100){ 
            nodes.df[iNode,'Reads'] <- paste(all_LG_seq_names,collapse = ' ')
         } else {
            nodes.df[iNode,'Reads'] <- paste(all_LG_seq_names[1:100],'...',collapse = ' ')
         }
      }
   }
   
   node_list <- nodes.df$Node
   # Loop over transcripts ----
   if (nrow(transcripts.df)>0){
      for (i_tx in 1:nrow(transcripts.df)){
         contig_comp <- unlist(strsplit(transcripts.df[i_tx,'Composition'],'-\\(0\\)->'))
         read_list <- ''
         for (j_node in 1:length(contig_comp)){
            # tx_node_hash.mat <- matrix(0, ncol = nrow(nodes.df),nrow = nrow(transcripts.df))
            # rownames(tx_node_hash.mat) <- transcripts.df$Name; colnames(tx_node_hash.mat) <- nodes.df$Node
            
            node_no <- as.character(abs(as.numeric(strsplit(contig_comp[j_node],':')[[1]][1])))
            
            tx_node_hash.mat[i_tx,node_no] <- 1
            
            node_associated_reads <- seq_names[read_node_hash.mat[,node_no] == 1]
            
            read_list <- c(read_list,node_associated_reads)
            
            
            read_tx_hash.mat[node_associated_reads,i_tx] <- 1
         }
         
         if (any(read_list == '')){read_list <- read_list[-which(read_list == '')]}
         read_list <- sort(unique(read_list))
         if (length(read_list) > 100){
            transcripts.df[i_tx,'Reads'] <- paste(read_list[1:100],'...',collapse = ' ')
         } else {
            transcripts.df[i_tx,'Reads'] <- paste(read_list,collapse = ' ')
         }
         
         tx_nodes <- as.character(node_list[tx_node_hash.mat[i_tx,] > 0])
         if (length(tx_nodes) == 1){
            read_list_from_hash_tab <- names(read_node_hash.mat[,tx_nodes] > 1)
         } else {
            read_list_from_hash_tab <- names(rowSums(read_node_hash.mat[,tx_nodes]) > 1)
         }
         
         
         if (length(read_list) != sum(read_tx_hash.mat[,i_tx])){
            message(i_tx)
         }
         
         transcripts.df[i_tx,'nReads'] <- length(read_list)
         # transcripts.df[i,'AB_shared'] <- sum(sub('A','B',a_reads) %in% b_reads)
         # transcripts.df[i,'AB_shared_prop'] <- transcripts.df[i,'AB_shared'] / transcripts.df[i,'nReads']
         
         reads.df[reads.df$Read %in% read_list,'nTX'] <- reads.df[reads.df$Read %in% read_list,'nTX'] + 1
         
         #try(reads.df[reads.df$Read %in% read_list,'Transcript'] <- paste(unique(c(reads.df[reads.df$Read %in% read_list,'Transcript'] ,transcripts.df[i_tx,'Name'])),collapse=''))
         #reads.df[reads.df$Read %in% read_list,'nTX'] <- reads.df[reads.df$Read %in% read_list,'nTX'] + rep(1,length(reads.df[reads.df$Read %in% read_list,'nTX']))
      }
   }
   
   if (merged){
      event_tx_dir <- paste(unlist(strsplit(tx_out_dir,'/'))[1:(grep('merged',unlist(strsplit(tx_out_dir,'/')))- 1)],collapse='/')
      read_names_raw <- names(read.fasta(paste(event_tx_dir,'/Sequences',sep='')))
      read_names <- regmatches(read_names_raw,regexpr('[.\\/\\w]{4,}(?=\\t)',read_names_raw,perl=TRUE))
      
      merged_read_tx.hash <- matrix(0,nrow = length(read_names),ncol = nrow(transcripts.df))
      colnames(merged_read_tx.hash) <- transcripts.df$Name; rownames(merged_read_tx.hash) <- read_names
      
      for (i_m_tx in 1:nrow(transcripts.df)){
         
         tx_list <- seq_names[read_tx_hash.mat[,i_m_tx] == 1]
         for (j_tx in tx_list){
            j_tx_kmer <- regmatches(j_tx,regexpr('(?<=_)\\d+(?=_Loc)',j_tx,perl=TRUE))
            hash_rds_file <- file.path(event_tx_dir,paste('k',j_tx_kmer,sep='_'),'read_tx_hash.Rds')
            j_read_tx.hash <- readRDS(hash_rds_file)
            j_tx_reads <- rownames(j_read_tx.hash)[j_read_tx.hash[,j_tx] == 1]
            
            merged_read_tx.hash[j_tx_reads,i_m_tx] <- 1
         }
         
      }
      output <- list(reads.df = reads.df, 
                     nodes.df = nodes.df, 
                     transcripts.df = transcripts.df,
                     read_node.hash = read_node_hash.mat,
                     tx_node.hash = tx_node_hash.mat,
                     read_tx.hash = read_tx_hash.mat,
                     merged_read_tx.hash = merged_read_tx.hash)
      read.tx.hash.filename <- file.path(tx_out_dir,'merged_read_tx_hash.Rds')
      saveRDS(merged_read_tx.hash,read.tx.hash.filename)
      
   } else {
      output <- list(reads.df = reads.df, 
                     nodes.df = nodes.df, 
                     transcripts.df = transcripts.df,
                     read_node.hash = read_node_hash.mat,
                     tx_node.hash = tx_node_hash.mat,
                     read_tx.hash = read_tx_hash.mat)
      read.tx.hash.filename <- file.path(tx_out_dir,'read_tx_hash.Rds')
      saveRDS(read_tx_hash.mat,read.tx.hash.filename)
   }
   
   # for (i_tx in 1:nrow(reads.df)){
   #    #reads.df[i,'nContig'] <- length(unlist(strsplit(reads.df[i,'Contig'], ' ')))
   #    reads.df[i_tx,'nTX'] <- length(unlist(strsplit(reads.df[i_tx,'Transcript'], ' ')))
   # }
   
   
   
   return(output)
   
}

