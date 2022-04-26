
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
#
# spec = matrix(c(
#    'reads.file', 'f',"character","location of file containing RNA spanning reads pairs",NA,
#    'k','k','character','hash lengths for Velveth',paste(seq(15,31,by=4),collapse=','),
#    'reads.threshold','r','integer','minimum number of spanning read pairs required for assembly',0,
#    'insert.length','i','numeric','insert length of RNA library',NA,
#    'cov_cutoff','c','character','coverage cutoff for velvet','dynamic',
#    'flash.threshold','fl','integer','number of read pairs to trigger quick combining of read pairs',50),
#    byrow = TRUE,ncol = 5)
#
# debug = FALSE
# if (debug){
#    readsFilename = '/research/labs/experpath/vasm/shared/NextGen/Projects/ExVivoEndom/RNAseq/EM28502/GRCh38/svar-1/dev/EM102_RNA_supportForDNAseqJcts_reads_EM102_MP.csv'
#    tx.dir = '/research/labs/experpath/vasm/shared/NextGen/Projects/ExVivoEndom/RNAseq/EM28502/GRCh38/svar-1/dev/tx'
#    file.prefix = 'EM102_RNA'
#    k = paste(seq(15,31,by=4),collapse=',')
#    reads.threshold = 0
#    cov_cutoff = 'dynamic'
#    flash.threshold = 50
#    insert_length = NULL
# }
#
# velvet_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/velvet/1.2.10'
# oases_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/oases'
# flash_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/FLASH-1.2.11/flash'
# perl_path <- '/research/labs/experpath/vasm/shared/NextGen/judell/lib/perl5/perls/perl-ju/bin/perl'
# vo_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/VelvetOptimiser/2.2.6/VelvetOptimiser.pl'
# old_path <- Sys.getenv("PATH")
# if (!grepl('velvet',old_path)){Sys.setenv(PATH = paste(old_path, velvet_path, sep = ":"))}
# rw.pause.time <- 0.01

# Functions -----



intruderAssembly <- function(reads.filename,tx.dir,file.prefix,
                             reads1 = NA,reads2 = NA,
                             aligned.bam = NA,
                             chimeric.bam = NA,
                             insert_length = NA,
                             flash.threshold = 50,
                             k = seq(15,31,by=4),
                             reads.threshold = 0,
                             cov_cutoff = 'dynamic',
                             overwrite_tx = TRUE,
                             overwrite_plots = FALSE,
                             overwrite_raw_files = FALSE,
                             # take.rc = FALSE,
                             cleanup = FALSE,
                             numCores = 12,
                             all_reads = FALSE,
                             tpm_threshold = 30,
                             use_preRR = TRUE,
                             method = 'spades') {
   # options for "method" parameter: "velvet_optimizer","transrate"
   # Setup variables ----
   # k = seq(15,31,by=4); COV_CUTOFF = 'dynamic'; reads.threshold = 0;
   #   flash.threshold = 50;  overwrite_tx = TRUE;

   #optimize_velvet <- method == 'velvet_optimizer'


   # COV_CUTOFF <- cov_cutoff
   seq.df <- loadRNAReads(reads.filename)
   # message(paste('Loaded ',nrow(seq.df), ' spanning RNA RPs','...',sep=''))

   if (nrow(seq.df) == 0){
      message('No RNA reads found for any junctions.')
      return(NULL)
   }

   # output[['seq.df']] <- seq.df


   if (!dir.exists(tx.dir)){dir.create(tx.dir)}

   # tx.df.filename <- file.path(tx.dir,paste(file.prefix,'TX.csv',sep='_'))

   event_list <- as.integer(names(table(seq.df$DNA_eventId))[table(seq.df$DNA_eventId) >= reads.threshold])

   if (length(event_list) == 0){message('No events meet reads threshold.')}

   check_prot <- FALSE

   output_tx <- data.frame()#data.frame(unique(inputSeqDF[,!names(inputSeqDF) %in% c('seqIdNum','ReadA','ReadB')]))
   tx.gr <- GRanges()
   df_header <- c('DNA_sampleId','DNA_eventId','RNA_RPs','chrA','chrB','DNA_posA','DNA_posB','geneA','geneB')

   DNA_sampleId = unique(seq.df$DNA_sampleId)

   # Pre-process event info & runSpades() in parallel ----

   # check_spades = TRUE
   use_preRR <- TRUE

   if (!overwrite_tx){
      spanning_dir_list <- file.path(tx.dir,paste0('event_',formatC(x = event_list,format = 'd',flag = 0,width = 3)),'spanning')


      good_tx_fa_list <- file.path(tx.dir,paste0('event_',formatC(x = event_list,format = 'd',flag = 0,width = 3)),'good.transcripts.fa')
      event_list <- event_list[file.exists(good_tx_fa_list)]
   }
   tx_progress_logfile <- file.path(tx.dir,'tx.log')


   cat(paste('runSpades called on',format(Sys.time()),'over',numCores,'threads:','\n'),append = FALSE,file = tx_progress_logfile)
   tic()
   closeAllConnections()
   cl <- makeForkCluster(nnodes = min(numCores,length(event_list)))
   # if (length(event_list) < numCores){register(MulticoreParam(workers = 4),default = TRUE)} else {register(MulticoreParam(workers = numCores),default = TRUE)}
   dna_prep_out <- clusterApply(cl = cl,x = event_list,fun = function(i) eventwiseDeNovoAssemblyPrep(seq.df = seq.df,event = i,tx.dir = tx.dir,overwrite = overwrite_raw_files,chimeric.bam = chimeric.bam,aligned.bam = aligned.bam) )
   runSpades_out <- clusterApply(cl = cl,x = event_list,
                                 fun  = function(event) tryCatch(runSpades(event = event,
                                                                           seq.df = seq.df,
                                                                           tx.dir = tx.dir,
                                                                           overwrite_tx = overwrite_tx,
                                                                           use_preRR = use_preRR,
                                                                           pause = TRUE),
                                                                 error = function(err_msg) {cat(paste('Event',event,'first pass assembly failed with error:',err_msg),sep='\n',file = tx_progress_logfile,append = TRUE); return('error')},
                                                                 warning = function(warn_msg) {cat(paste('Event',event,'first pass assembly failed with warning:',warn_msg),sep='\n',file = tx_progress_logfile,append = TRUE); return('warning')},
                                                                 finally = cat(paste('Event',event,'yielded',
                                                                                     ifelse(test = file.exists(file.path(tx.dir,paste0('event_',formatC(x = event,format = 'd',flag = 0,width = 3)),'good.transcripts.fa')),
                                                                                            yes = countFastaEntries(file.path(tx.dir,paste0('event_',formatC(x = event,format = 'd',flag = 0,width = 3)),'good.transcripts.fa')),
                                                                                            no = 0),'transcripts.\n'),file = tx_progress_logfile,append = TRUE)))
   stopCluster(cl = cl)
   toc()
   spades_out_files <- system(command = paste('ls',file.path(tx.dir,'event_*/good.transcripts.fa')),intern=TRUE)

   remaining_events <- event_list[sapply(X = event_list,FUN = function(i) !any(grepl(formatC(x = i,width = 3,format = 'd',flag = 0),spades_out_files)),USE.NAMES = FALSE)]
   n_cycles = 0
   numCores_in_loop <- numCores
   # cl <- makeForkCluster(nnodes = numCores/2)
   # length(remaining_events)
   # message(length(remaining_events),' events left to be assembled...')
   while (length(remaining_events) > 0 & n_cycles < 5){

      numCores_in_loop <- ifelse(test = numCores_in_loop <= length(remaining_events),
                                 yes = numCores_in_loop,
                                 no = max(numCores_in_loop/2,1))
      cl <- makeForkCluster(nnodes = numCores_in_loop)

      message('Event(s) ',paste(remaining_events,collapse=', '),' events left to be assembled...')

      runSpades_out <- clusterApply(cl = cl,x = event_list,
                                    fun  = function(event) tryCatch(runSpades(event = event,
                                                                              seq.df = seq.df,
                                                                              tx.dir = tx.dir,
                                                                              overwrite_tx = overwrite_tx,
                                                                              use_preRR = use_preRR,
                                                                              pause = TRUE),
                                                                    error = function(err_msg) {cat(paste('Event',event,'second+ pass assembly failed with error:',err_msg),sep='\n',file = tx_progress_logfile,append = TRUE); return('error')},
                                                                    warning = function(warn_msg) {cat(paste('Event',event,'second+ pass assembly failed with warning:',warn_msg),sep='\n',file = tx_progress_logfile,append = TRUE); return('warning')},
                                                                    finally = cat(paste('Event',event,'yielded',
                                                                                        ifelse(test = file.exists(file.path(tx.dir,paste0('event_',formatC(x = event,format = 'd',flag = 0,width = 3)),'good.transcripts.fa')),
                                                                                               yes = countFastaEntries(file.path(tx.dir,paste0('event_',formatC(x = event,format = 'd',flag = 0,width = 3)),'good.transcripts.fa')),
                                                                                               no = 0),'transcripts.\n'),file = tx_progress_logfile,append = TRUE)))


      # try(runSpades_out <- clusterApply(cl = cl,x = remaining_events,fun  = function(event) try(runSpades(event = event,seq.df = seq.df,tx.dir = tx.dir,overwrite_tx = overwrite_tx,use_preRR = use_preRR,pause = TRUE))))
      spades_out_files <- system(command = paste('ls',file.path(tx.dir,'event_*/good.transcripts.fa')),intern=TRUE)
      remaining_events <- event_list[sapply(X = event_list,FUN = function(i) !any(grepl(formatC(x = i,width = 3,format = 'd',flag = 0),spades_out_files)),USE.NAMES = FALSE)]
      n_cycles = n_cycles + 1
      stopCluster(cl = cl)
   }

   # runSpades_out <- clusterApply(cl = cl,x = remaining_events,fun  = function(event) try(runSpades(event = event,seq.df = seq.df,tx.dir = tx.dir,overwrite_tx = overwrite_tx,check_spades = check_spades,use_preRR = use_preRR,pause = TRUE)))

   # stopCluster(cl = cl)
   toc()

   # Cat results ----
   ## TX.csv
   all_TX.csv <- file.path(tx.dir,paste0(file.prefix,'_TX.csv'))
   TX.csv_file_list <- sapply(X = event_list,FUN = function(i) file.path(tx.dir,paste0('event_',formatC(x = i,width = 3,flag = 0,format = 'd')),'TX.csv')) # system(command = paste('ls',file.path(tx.dir,'event_*/TX.csv')),intern=TRUE)
   TX.csv_file_list <- TX.csv_file_list[file.exists(TX.csv_file_list)]
   # cat tx.dir/event_*/TX.csv[1] > all_TX.csv && cat $(head --lines=-1 tx.dir_event/TX.csv[all-but-1]) >> all_TX.csv

   if (length(TX.csv_file_list) > 0){
      system(command = paste('cat',TX.csv_file_list[1],'>',all_TX.csv,'&& awk FNR-1',paste(TX.csv_file_list[-1],collapse=' '),' >>',all_TX.csv))
      tx.df <- read.table(all_TX.csv,header = TRUE,sep=',',stringsAsFactors = FALSE)#; tx.df <- tx.df[-which(tx.df$DNA_sampleId == 'DNA_sampleId'),]

      ## good.transcripts.fa
      all_assembled_tx.fa <- file.path(tx.dir,'transcripts.fa')
      system(command = paste('cat',file.path(tx.dir,'event_*/good.transcripts.fa'),'>',all_assembled_tx.fa))

      tx.df$goodTranscript <- tx.df$tx_num %in% getFastaHeader(fasta_file = all_assembled_tx.fa)
   } else {
      message('No TX loaded!')
      tx.df <- data.frame()
   }



   return(tx.df)

}

getFastaHeader <- function(fasta_file){
   # grep "^>" <filename> | sed 's/;/,/g'
   header_list <- system(command = paste('grep "^[>@]"',fasta_file, '| sed "s/[>@]//g"'),intern=TRUE)
   return(header_list)
}

countFastaEntries <- function(fasta_file){
   # grep "^>" <filename> | sed 's/;/,/g'
   if (file.exists(fasta_file) && file.size(fasta_file) > 0){
      header_count <- as.integer(system(command = paste('grep -c "[>@]"',fasta_file),intern=TRUE))
   } else {
      header_count <- 0
   }

   return(header_count)
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
      if (!any(grepl('RNA_RPs',names(seq.df)))){
         seq.df$RNA_RPs <- numeric(length = nrow(seq.df))
      }
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


            seq.df$RNA_RPs[which(seq.df$DNA_eventId == event_list[j])] <- rep(length(unique(jct$seqID)),length(which(seq.df$DNA_eventId == event_list[j])))
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

spanningRPsPreProcessing <- function(seq.df,event,event.tx.dir,
                                     flash.threshold = 100,
                                     remove.polyA = FALSE,
                                     prehash_reads = FALSE,
                                     write.fq = TRUE){

   e_seq.df <- seq.df[seq.df$DNA_eventId == event,]
   # event.readPairs.fa <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'_pairs.fa',sep=''))
   # event.flashed.reads.fa <- file.path(event.tx.dir,paste(file.prefix,'_event_',event,'_flash.fa',sep=''))
   if (write.fq){
      span1.fq <- file.path(event.tx.dir,'span.1.fq')
      span1.ss <- DNAStringSet(x = e_seq.df$readA); names(span1.ss) <- e_seq.df$seqID
      writeXStringSet(x = span1.ss,filepath = span1.fq,format = 'fastq')

      span2.fq <- file.path(event.tx.dir,'span.2.fq')
      span2.ss <- DNAStringSet(x = e_seq.df$readB); names(span2.ss) <- e_seq.df$seqID
      writeXStringSet(x = span2.ss,filepath = span2.fq,format = 'fastq')
   } else {
      span1.fq = ''; span2.fq = ''
   }

}

runVelvet <- function(seq.df,tx.dir,file.prefix){
   event_list <- sort(unique(seq.df$DNA_eventId[seq.df$RNA_RPs >= reads.threshold]))

   e_seq.df <- seq.df[seq.df$DNA_eventId == event,]
   event.tx.dir <- file.path(tx.dir,paste('event',event,sep='_'))
   if (!dir.exists(event.tx.dir)){dir.create(event.tx.dir)}


   n_spann_RPs <- nrow(e_seq.df)

   # Using VO ----
   # flash.flag = readsPP.out$n_combined > 0
   ## test assembly ----
   tx <- velvetOases(iJct = event,kmer = 15,
                     event.tx.dir = event.tx.dir,
                     file.prefix = file.prefix,
                     flash.flag = FALSE)

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
   } else {cov_cutoff = 'auto'}

   if (velvet_check & optimize_velvet){
      #tic()
      ### Call Velvet Optimiser ----
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
      ## else, run regular velvet ------
      message(paste('Assembling ',n_spann_RPs,' read pairs at event ',event,'...',sep=''))
      ### Loop over kmers ----
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

         ### Remove transcripts that only use reference sequence ----
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

      ### Write assembled novel TX to file in event dir for assembly ----
      tx.ss <- DNAStringSet(x = output_tx$TX[output_tx$DNA_eventId == event ],use.names = TRUE)
      names(tx.ss) <- output_tx$tx_num[output_tx$DNA_eventId == event & output_tx$kmer_size != 'merged']
      tx.fa_file <- file.path(event.tx.dir,'raw_assembled_transcripts.fa')
      writeXStringSet(x = tx.ss,filepath = tx.fa_file)

      ### Merge assembled TX ----
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
      ## Parsing contigs of merged TX ----
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


            ### format output data ----
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



   if (optimize_velvet){tx.df <- output_tx} else {tx.df <- merged_tx}


   ## Cleanup -----
   if (cleanup){
      merged.tx.dir <- paste(tx.dir,'merged',sep='/')
      unlink(merged.tx.dir,recursive = TRUE)
      ##kmer_dirs <- paste(tx.dir,'/k_',k,sep='')
      for (kmer in k){
         kmer_dir <- paste(tx.dir,'/k_',k,sep='')
         unlink(kmer_dir,recursive = TRUE)
      }
   }

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
         nodes.df[iNode,'nReads'] <- length(unique(all_LG_seq_names))
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

runSpades <- function(seq.df,tx.dir,event,overwrite_tx,use_preRR = FALSE,pause = FALSE,CPU = 1){
   if (pause){Sys.sleep(time = abs(rnorm(1,mean = .5)))}
   event0 <- formatC(x = event,width =3,flag = 0,format = 'd')
   e_seq.df <- seq.df[seq.df$DNA_eventId == event,]
   DNA_sampleId <- unique(e_seq.df$DNA_sampleId)
   event.tx.dir <- file.path(tx.dir,paste('event',event0,sep='_'))
   if (!dir.exists(event.tx.dir)){dir.create(event.tx.dir,recursive = TRUE)}

   n_spann_RPs <- nrow(e_seq.df)

   jct.gr <- getGRangesFromInputDF(seq.df = seq.df,event = event)

   event.dir = event.tx.dir
   tx.fa <- file.path(event.dir,'good.transcripts.fa')
   tx.stats <- file.path(event.dir,'good.transcripts.stats.csv')
   tx.out.csv <- file.path(event.dir,'TX.csv')

   event_input.dir <- file.path(event.dir,'raw_input')

   event.bam <- file.path(event_input.dir,'event.bam')
   span1.fq <- file.path(event_input.dir,'span.1.fq')
   span2.fq <- file.path(event_input.dir,'span.2.fq')

   run_spades <- !file.exists(tx.fa) | !file.exists(tx.stats) | !file.exists(tx.out.csv) | overwrite_tx

   if (run_spades){

      min.tx.length = 200

      n_reads <- system(command = paste('grep -c "@"',file.path(event_input.dir,'all.1.fq')),intern=TRUE)

      unlink(x = c(tx.fa,tx.stats,tx.out.csv,file.path(event.dir,'spades'),file.path(event.dir,'spanning'),file.path(event.dir,'transrate')))

      spades.out <- spades(event = event,event.bam = event.bam,jct.gr = jct.gr,event_input.dir = event_input.dir,
                           event.dir = event.dir,use_preRR = use_preRR,CPU = CPU,
                           overwrite_tx = overwrite_tx)
   }

   if (file.exists(tx.fa) && file.size(tx.fa) > 0){
      n_tx <- as.integer(system(paste('grep -c ">"',tx.fa),intern = TRUE))
   } else {
      unlink(x = file.path(event.dir,'TX.csv'))
      n_tx = 0
   }

   if (n_tx > 0){
      e_stats.df <- read.table(file = file.path(event.dir,'good.transcripts.stats.csv'),sep=' ',header=TRUE)

      e_tx.df <- data.frame('DNA_sampleId' = rep(DNA_sampleId,times=n_tx),
                            'DNA_eventId' = rep(event,n_tx),
                            'RNA_RPs' = rep(n_spann_RPs,n_tx),
                            'chrA' = rep(seqnames(jct.gr[1]),n_tx),
                            'chrB' = rep(seqnames(jct.gr[2]),n_tx),
                            'posAmin' = rep(start(jct.gr[1]),n_tx),
                            'posAmax' = rep(end(jct.gr[1]),n_tx),
                            'posBmin' = rep(start(jct.gr[2]),n_tx),
                            'posBmax' = rep(end(jct.gr[2]),n_tx),
                            'DNA_posA' = rep(unique(e_seq.df$DNA_posA),n_tx),
                            'DNA_posB' = rep(unique(e_seq.df$DNA_posB),n_tx),
                            'geneA' = rep(unique(e_seq.df$geneA),n_tx),
                            'geneB' = rep(unique(e_seq.df$geneB),n_tx))
      e_tx.df <- cbind(e_tx.df,e_stats.df)
      e_tx.df$tx_num <- e_tx.df$contig_name
      write.table(x = e_tx.df, file = file.path(event.dir,'TX.csv'),quote = FALSE,row.names = FALSE,sep=',')

   }

   # output event summary ----

   n_RNA_RPs <- countFastaEntries(file.path(event_input.dir,'all.1.fq'))
   bam.stats <- capture.output(quickBamFlagSummary(file.path(event_input.dir,'event.bam')))
   n_contigs <- countFastaEntries(file.path(event.dir,'spanning/good_spanning_contigs.fa'))
   n_before_rr <- countFastaEntries(file.path(event.dir,'spanning/spades/before_rr.fasta'))

   n_raw_spades_out <- countFastaEntries(file.path(event.dir,'spanning/spades/transcripts.fasta'))


   n_hf_spades_out <- countFastaEntries(file.path(event.dir,'spanning/spades/hard_filtered_transcripts.fasta'))


   summary_str <- c(paste('Event:\t',event),
                    paste0('Region A:\t',as.character(jct.gr[1])),
                    paste0('Region B:\t',as.character(jct.gr[2])),
                    paste0('nLocalReads:\t',n_RNA_RPs),
                    '\n-------- BAM stats --------',
                    bam.stats,
                    '\n---- Spanning Assembly ----',
                    paste0('Num spanning RPs:\t',n_spann_RPs),
                    paste0('Assembled spanning contigs:\t',n_raw_spades_out,' (',n_raw_spades_out - n_hf_spades_out,' low quality, ',n_before_rr - n_raw_spades_out,' before RR)'))



   if (file.exists(file.path(event.dir,'spanning/transrate/assemblies.csv'))){
      spanning_spades_out <- read.table(file = file.path(event.dir,'spanning/transrate/assemblies.csv'),header = TRUE,sep=',')

      n_good_spanning_contigs <- countFastaEntries(file.path(event.dir,'spanning/good_spanning_contigs.fa'))
      n_fragments_mapped <- spanning_spades_out[1,'fragments_mapped']

      spanning_transrate_line <- paste0('Transrate metrics:\n  o Good contigs: ',n_good_spanning_contigs,'\n  o Fragments mapped: ',spanning_spades_out[1,'fragments_mapped'],' / ',spanning_spades_out[1,'fragments'],' (',format(x = 100*spanning_spades_out[1,'p_fragments_mapped'],digits = 3),'%)',
                                        '\n  o Good mappings: ',spanning_spades_out[1,'good_mappings'],'\n  o Assembly score: ',spanning_spades_out[1,'score'],'; optimal score: ',spanning_spades_out[1,'optimal_score'])
      summary_str <- c(summary_str,spanning_transrate_line)
   } else {n_good_spanning_contigs <- 0}

   n_contigs <- countFastaEntries(file.path(event.dir,'good.transcripts.fa'))
   n_raw_spades_out <- countFastaEntries(file.path(event.dir,'spades/transcripts.fasta'))
   n_hf_spades_out <- countFastaEntries(file.path(event.dir,'spades/hard_filtered_transcripts.fasta'))
   n_before_rr <- countFastaEntries(file.path(event.dir,'spades/before_rr.fa'))

   summary_str <- c(summary_str,
                    '\n---- General Assembly ----',
                    paste0('Assembled contigs:\t',n_raw_spades_out,' (',n_raw_spades_out - n_hf_spades_out,' low quality, ',n_before_rr - n_before_rr,' removed via RR)'))


   if (n_contigs > 0){
      transrate_spades_out <- read.table(file = file.path(event.dir,'transrate/assemblies.csv'),header = TRUE,sep=',')

      n_fragments_mapped <- transrate_spades_out[1,'fragments_mapped']

      transrate_str <- paste0('Transrate metrics:\n  o Good contigs: ',n_good_spanning_contigs,'\n  o Fragments mapped: ',transrate_spades_out[1,'fragments_mapped'],' / ',transrate_spades_out[1,'fragments'],' (',format(x = 100*transrate_spades_out[1,'p_fragments_mapped'],digits = 3),'%)',
                              '\n  o Good mappings: ',transrate_spades_out[1,'good_mappings'],'\n  o Assembly score: ',transrate_spades_out[1,'score'],'; optimal score: ',transrate_spades_out[1,'optimal_score'])
      summary_str <- c(summary_str,transrate_str)
   } else {
      n_fragments_mapped <- 0
   }



   tx.log <- file.path(event.dir,'tx.log')
   writeLines(text = summary_str,con = tx.log)
   if (n_contigs > 0){
      message('Spades assembly successful for event ',event)
      tx.df <- read.table(file =  file.path(event.dir,'good.transcripts.stats.csv'),header = TRUE)
      return(tx.df)
   } else {
      message('Spades assembly failed for ',event)
      return()
   }
}

spades <- function(event,event.bam,event.dir,jct.gr,
                   event_input.dir,
                   overwrite_tx = TRUE,
                   min.tx.length = 200,
                   use_preRR = FALSE,
                   filterByCDS = FALSE,
                   CPU = 1,
                   use_conda = FALSE){
   # numCores = 1 # because the whole process is run on multiple cores
   # use_conda = FALSE
   if (use_conda){
      spades_path <- 'conda run -n spades_sandbox rnaspades.py'
   } else {
      spades_path <- 'python /research/labs/experpath/vasm/shared/NextGen/judell/lib/SPAdes-3.15.4-Linux/bin/rnaspades.py'
   }

   event0 <- formatC(x = event,width = 3,flag = 0,format = 'd')

   all1.fq <- file.path(event_input.dir,'all.1.fq')
   all2.fq <- file.path(event_input.dir,'all.2.fq')
   span1.fq <- file.path(event_input.dir,'span.1.fq')
   span2.fq <- file.path(event_input.dir,'span.2.fq')

   dn.out <- list('tx.ss' = DNAStringSet(),'tx.gr' = GRanges(),'tx.df' = data.frame())

   # Spanning read assembly ----
   ## Prep ----

   spanning.spades.dir <- file.path(event.dir,'spanning/spades')
   if (!dir.exists(spanning.spades.dir)){dir.create(spanning.spades.dir,recursive = TRUE)}
   spanning_fa <- file.path(spanning.spades.dir,'transcripts.fasta')

   ## Assembly ----

   if (!file.exists(spanning_fa) | overwrite_tx){
      unlink(spanning.spades.dir,recursive = TRUE)

      dir.create(spanning.spades.dir,recursive = TRUE)
      spades.log <- file.path(spanning.spades.dir,'stdout.txt')

      spades.cmd <- paste(spades_path,'-o',spanning.spades.dir,'-t',CPU,
                          '-1', span1.fq,'-2',span2.fq,' > ',spades.log)
      suppressMessages(system(spades.cmd))

      if (any(grepl(pattern='terminate called without an active exception|finished abnormally',x = readLines(spades.log)))){
         spades.cmd <- paste(spades_path,'-o',spanning.spades.dir,
                             '--continue',' &> ',spades.log)
         system(spades.cmd)
      }



      spades_preRR.fa  <- file.path(spanning.spades.dir,'before_rr.fasta')
      if (!file.exists(spanning_fa) & file.exists(spades_preRR.fa) & use_preRR){
         spanning_fa = spades_preRR.fa
      }

   }

   if (file.exists(spanning_fa) && file.size(spanning_fa) > 0){
      n_spanning_contigs <- as.integer(system(paste('grep ">" -c',spanning_fa),intern=TRUE))
   } else {
      n_spanning_contigs = 0
   }

   filtered_spanning_fa <- file.path(event.dir,'spanning/good_spanning_contigs.fa')
   transrate.dir <- file.path(event.dir,'spanning/transrate')
   run_transrate <- !file.exists(filtered_spanning_fa) | overwrite_tx
   if (n_spanning_contigs > 0 & run_transrate ){
      ## Transrate ----
      unlink(filtered_spanning_fa)

      unlink(transrate.dir,recursive = TRUE); dir.create(transrate.dir,recursive = TRUE)
      transrate_assemblies_out <- file.path(transrate.dir,'assemblies.csv')

      unlink(transrate_assemblies_out)

      transrate_input.fa <- file.path(transrate.dir,'raw_spanning_contigs.fa')
      unlink(transrate_input.fa)
      file.copy(from = spanning_fa,to = transrate_input.fa)

      suppressMessages(spanning.contigs.df <- transrate(fasta.list = transrate_input.fa,
                                                        f1.fq = span1.fq,f2.fq = span2.fq,CPU=CPU,transrate.dir = transrate.dir,outfile = filtered_spanning_fa))

   } else if (n_spanning_contigs == 0 ) {
      unlink(filtered_spanning_fa)
      unlink(transrate.dir,recursive = TRUE)
   }

   if (file.exists(filtered_spanning_fa) && file.size(filtered_spanning_fa) > 0){
      n_filtered_spanning_contigs = as.integer(system(paste('grep ">" -c',filtered_spanning_fa),intern=TRUE))
   } else {
      n_filtered_spanning_contigs = 0
   }

   # Spades assembly ----
   spades.dir <- file.path(event.dir,'spades')
   if (!dir.exists(spades.dir)){dir.create(spades.dir)}
   spades_fa <- file.path(event.dir,'spades/spades.fa')

   if (!file.exists(spades_fa) | overwrite_tx){
      unlink(spades.dir,recursive = TRUE)
      dir.create(spades.dir)
      spades_output <- file.path(spades.dir,'transcripts.fasta')
      #unlink(spades_output)
      if (n_filtered_spanning_contigs == 0){
         trusted_contigs_param <- ''
      } else {
         trusted_contigs_param <- paste('--trusted-contigs',filtered_spanning_fa)
      }

      spades.log <- file.path(spades.dir,'stdout.txt')
      spades.cmd <- paste(spades_path,'-o',spades.dir,
                          '-1',all1.fq,'-2',all2.fq,trusted_contigs_param,'-t',CPU,
                          '>',spades.log)
      system(spades.cmd)

      if (any(grepl(pattern='terminate called without an active exception|finished abnormally',x = readLines(spades.log)))){
         spades.cmd <- paste(spades_path,'-o',spades.dir,
                             '--continue',
                             '>',spades.log)
         system(spades.cmd)
      }

      if (file.exists(spades_output) && file.size(spades_output) > 0){
         tx <- read.fasta(file = spades_output,
                          as.string = TRUE,set.attributes = FALSE)
         names.df <- data.frame(strsplit(names(tx),'_'))
         spades_tx.df <- data.frame('names' = names(tx),
                                    'node' = as.character(names.df[2,]),
                                    'length' = as.integer(names.df[4,]),
                                    'cov' = as.numeric(names.df[6,]),
                                    'gene' = as.integer(sub('g','',names.df[7,])),
                                    'Loc' = sub('i','',names.df[8,]),
                                    'seq' = as.character(tx),
                                    stringsAsFactors = FALSE)
         spades_tx.df$gene_obvs <- sapply(X = spades_tx.df$gene,FUN = function(i) sum(spades_tx.df$gene == i))
         spades_tx.df$tx_num <- paste('J',event0,'Loc',spades_tx.df$gene,'Iso',spades_tx.df$Loc,sep='_')
         spades.ss <- DNAStringSet(x = spades_tx.df$seq); names(spades.ss) <- spades_tx.df$tx_num
         writeXStringSet(x = spades.ss,filepath = spades_fa)
      } else {

         n_assembled_contigs = 0
      }
   } else {
      spades.ss <- readDNAStringSet(filepath = spades_fa,use.names = TRUE)
   }

   if (file.exists(spades_fa) && file.size(spades_fa) > 0){
      n_assembled_contigs = length(spades.ss)
   } else {
      n_assembled_contigs = 0
   }

   if (n_assembled_contigs == 0){
      return(dn.out)
   }

   ## Transrate ----

   transrate.input.list <- c(spades_fa)
   transrate.dir <- file.path(event.dir,'transrate'); if (!dir.exists(transrate.dir)){dir.create(transrate.dir)}
   assembled_tx_out.fa <- file.path(transrate.dir,'transcripts.fa')
   assemblies.file <- file.path(transrate.dir,'assemblies.csv')
   if (n_assembled_contigs > 0 & (!file.exists(assembled_tx_out.fa) | overwrite_tx)){
      unlink(transrate.dir,recursive = TRUE)
      unlink(assemblies.file)
      fasta.list = transrate.input.list;f1.fq = all1.fq;f2.fq = all2.fq; outfile = assembled_tx_out.fa
      suppressWarnings(tx.out <- transrate(fasta.list = fasta.list,f1.fq = f1.fq,f2.fq = f2.fq,CPU=CPU,transrate.dir=transrate.dir,outfile = outfile))
   } else {
      tx.out <- read.csv(file = file.path(transrate.dir,'good_contigs.csv'),header = TRUE,stringsAsFactors = FALSE)
   }

   tx.ss <- readDNAStringSet(assembled_tx_out.fa)
   n_transrated_contigs <- length(tx.ss)

   # Wrangle output ----

   tx.df.file <- file.path(event.dir,'good.transcripts.stats.csv')
   out.fa <- file.path(event.dir,'good.transcripts.fa')
   if (n_transrated_contigs > 0){
      writeXStringSet(x = tx.ss,filepath = out.fa)
      tx.df = tx.out
      tx.df$nSpanningContigs = n_filtered_spanning_contigs
      write.table(x = tx.df,file = tx.df.file,quote = FALSE, row.names = FALSE)
   } else {
      unlink(x = c(tx.df.file,out.fa))
      # tx.gr <- GRanges()
      tx.df <- data.frame()
      tx.ss <- DNAStringSet()
   }

   dn.out <- list('tx.ss' = tx.ss,'tx.df' = tx.df)
   return(dn.out)
}



transrate <- function(fasta.list,f1.fq,f2.fq,CPU,transrate.dir,outfile = NULL){
   commands.file <- file.path(transrate.dir,'transrate_commands.txt')
   if (!dir.exists(transrate.dir)){dir.create(transrate.dir)}
   write(paste('# Transrate Code'),file= commands.file,append  = FALSE)
   if (any (!basename(fasta.list) %in% list.files(transrate.dir))){
      cp.cmd <- paste('cp',paste(fasta.list[!basename(fasta.list) %in% list.files(transrate.dir)],collapse=' '),transrate.dir)
      write(cp.cmd,file = commands.file,append = TRUE)
      system(cp.cmd)
   }


   tr.fasta.list <- file.path(transrate.dir,basename(fasta.list))


   f1.base <- strsplit(basename(f1.fq),'\\.')[[1]][1]
   f2.base <- strsplit(basename(f2.fq),'\\.')[[1]][1]
   snap.fasta.list <- c(tr.fasta.list)
   #tic()
   for (i.fa in snap.fasta.list){

      input.base <- strsplit(basename(i.fa),'\\.')[[1]][1]
      snap.input.fa <- i.fa

      snap.dir <- file.path(transrate.dir,input.base)
      snap.out.dir <- file.path(snap.dir,input.base); if (!dir.exists(snap.out.dir)){dir.create(snap.out.dir,recursive=TRUE)}
      snap.bam.file <- file.path(snap.dir,paste(basename(f1.fq),basename(f2.fq),input.base,'bam',sep='.'))
      snap.path <- '/research/labs/experpath/vasm/shared/NextGen/judell/home/miniconda3/envs/orp/bin/snap-aligner'
      #tr.input.fa <- file.path(mmseq.dir,'stv_rep_seq.fasta')
      snap.index.cmd <- paste(snap.path,'index', snap.input.fa,snap.out.dir,paste0('-t',CPU),'-s 23 -bSpace -locationSize 7')
      write(snap.index.cmd,file = commands.file,append = TRUE)
      system(snap.index.cmd,ignore.stderr = TRUE,ignore.stdout = TRUE)
      snap.align.cmd <- paste(snap.path,'paired',snap.out.dir,f1.fq,f2.fq,'-o',snap.bam.file,'-t',CPU,'-s 0 1000 -H 300000 -h 2000 -d 15 -b -M -D 5 -om 5 -omax 10')
      write(snap.align.cmd,file = commands.file,append = TRUE)
      system(snap.align.cmd,ignore.stderr = TRUE,ignore.stdout = TRUE)


   }
   #snap.dt <- toc()



   # Analysis via transrate (8.2 sec) ----

   tr.input.fa <- paste(snap.fasta.list,collapse=',')


   transrate.log <- file.path(transrate.dir,'transrate.log')
   transrate.cmd <- paste('conda run -n trans_test',
                          'transrate',
                          '--assembly',tr.input.fa,
                          '--left', f1.fq,'--right',f2.fq,
                          '--threads',CPU,
                          '--output',transrate.dir,
                          '--loglevel debug')# -o',transrate.dir)
   # conda create -n trans_test -y -c conda-forge -c bioconda transrate transrate-tools snap-aligner=1.0dev.96 salmon=0.6.0 boost=1.60.0
   #use_condaenv(condaenv = 'trans_test')
   #tic()
   write(transrate.cmd,file = commands.file,append = TRUE)
   loginfo <- system(transrate.cmd,intern = TRUE,ignore.stderr = TRUE)
   write(x = loginfo,file = transrate.log)
   #tr.dt <- toc()

   if (!is.null(outfile)) {
      assemblies.df <- read.table(file = file.path(transrate.dir,'assemblies.csv'),sep=',',header=TRUE)
      good_list <- sapply(X = fasta.list,FUN = function(i) file.path(transrate.dir,strsplit(basename(i),'\\.')[[1]][1],paste('good',basename(i),sep='.')), USE.NAMES = FALSE)

      #good_list <- file.path(transrate.dir,'*',paste0('good.',basename(fasta.list)))

      # out.fa <- file.path(transrate.dir,'good_contigs.fa')
      cat.cmd <- paste('cat',paste(good_list,collapse=' '),'>',outfile)
      write(cat.cmd,file = commands.file,append = TRUE)
      system(cat.cmd)

      good.ss <- read.fasta(file = outfile)
      contigs.df <- data.frame()

      for (i in 1:length(fasta.list)){
         contigs.filename <- file.path(transrate.dir,unlist(strsplit(basename(fasta.list[i]),'\\.'))[1],'contigs.csv')
         contigs.seq <- read.fasta(file = fasta.list[i],as.string = TRUE)
         i_contigs.df <- read.table(file = contigs.filename,header = TRUE,sep=',')
         i_contigs.df$tx = sapply(X = i_contigs.df$contig_name, FUN = function(i) toupper(as.character(contigs.seq[names(contigs.seq) == i])),USE.NAMES = FALSE)
         contigs.df <- rbind(contigs.df,i_contigs.df[i_contigs.df$contig_name %in% names(good.ss),])

      }
      if (nrow(contigs.df) > 0){
         write.csv(x = contigs.df,file = file.path(transrate.dir,'good_contigs.csv'),quote = FALSE,row.names = FALSE)
      }
   } else {contigs.df <- data.frame()}
   return(contigs.df)
}





deNovoAssemblyPrep <- function(event.dir,jct.gr,aligned.bam,chimeric.bam,overwrite = FALSE){

   #if (!dir.exists(event.dir)){dir.create(event.dir)}
   event.bam <- file.path(event.dir,'event.bam')
   event.bam.bai <- file.path(event.dir,'event.bam.bai')
   # eventSamFile <- file.path(raw.dir,'event.sam')
   event_gr_info <- file.path(event.dir,'event.txt')
   #old_jct.gr <- makeGRangesFromDataFrame(read.table(event_gr_info))

   #all1.fa <- file.path(event.dir,'all.1.fa')
   #all2.fa <- file.path(event.dir,'all.2.fa')
   all1.fq <- file.path(event.dir,'all.1.fq')
   all2.fq <- file.path(event.dir,'all.2.fq')


   # Check for chimeric.bam file ----


   bamFileWriteFlag <- !file.exists(event.bam) | !file.exists(event.bam.bai) | overwrite #| sum(ranges(old_jct.gr) != ranges(jct.gr)) > 0
   if (bamFileWriteFlag){
      # message('Writing event.bam')
      # Filter bam files according to event region (jct.gr)
      write.table(x = jct.gr,file = event_gr_info, quote = FALSE)
      # aligned.bamfile <- BamFile(file = aligned.bam,index = aligned.bam,asMates = TRUE)
      #event.aligned.unsorted.bamfile <- file.path(event.dir,'aligned.unsorted.bam')
      #sort.check <- system(command = paste('samtools stats',aligned.bam,'| grep "is sorted:"'))

      aligned.A.bam <-  file.path(event.dir,'aligned.A.bam')
      aligned.B.bam <-  file.path(event.dir,'aligned.B.bam')

      eventAlignedReadsA <- filterBam(file = aligned.bam,index = aligned.bam,
                                      destination = aligned.A.bam,
                                      #indexDestination = FALSE,
                                      param = ScanBamParam(which = jct.gr[1],
                                                           flag = scanBamFlag(isPaired = TRUE,
                                                                              isDuplicate = FALSE,
                                                                              isProperPair = TRUE,
                                                                              isUnmappedQuery = FALSE,
                                                                              hasUnmappedMate = FALSE),
                                                           what = c('flag','qname','rname','strand','pos','seq')))

      eventAlignedReadsB <- filterBam(file = aligned.bam,index = aligned.bam,
                                      destination = aligned.B.bam,
                                      #indexDestination = FALSE,
                                      param = ScanBamParam(which = jct.gr[2],
                                                           flag = scanBamFlag(isPaired = TRUE,
                                                                              isDuplicate = FALSE,
                                                                              isProperPair = TRUE,
                                                                              isUnmappedQuery = FALSE,
                                                                              hasUnmappedMate = FALSE),
                                                           what = c('flag','qname','rname','strand','pos','seq')))


      # aligned_reads <- scanBam(file = event.aligned.bamfile,index = event.aligned.bamfile,
      #                          param = ScanBamParam(what = c('flag','mapq','qname','rname','strand','pos','seq')))

      event.chimeric.A.bam <- file.path(event.dir,'chimeric.A.bam')
      event.chimeric.B.bam <- file.path(event.dir,'chimeric.B.bam')
      eventChimericReadsA <- filterBam(file = chimeric.bam,index = chimeric.bam,# region = jct.gr[1],
                                       indexDestinaton = FALSE, overwrite = TRUE,
                                       destination = event.chimeric.A.bam,
                                       param = ScanBamParam(which = jct.gr[1],
                                                            flag = scanBamFlag(isPaired = TRUE,
                                                                               isDuplicate = FALSE,
                                                                               isProperPair = TRUE,
                                                                               isUnmappedQuery = FALSE,
                                                                               hasUnmappedMate = FALSE),
                                                            what = c('mate_status','flag','mapq','qname','rname','strand','pos','seq')))
      eventChimericSortedReadsA <- sortBam(file = event.chimeric.A.bam,byQname = FALSE, destination = file.path(event.dir,'chimeric.A.sorted'))

      eventChimericReadsB <- filterBam(file = chimeric.bam,index = chimeric.bam,# region = jct.gr[2],
                                       destination = event.chimeric.B.bam,
                                       indexDestination = FALSE, overwrite = TRUE,
                                       param = ScanBamParam(which = jct.gr[2],
                                                            flag = scanBamFlag(isPaired = TRUE,
                                                                               isDuplicate = FALSE,
                                                                               isProperPair = TRUE,
                                                                               isUnmappedQuery = FALSE,
                                                                               hasUnmappedMate = FALSE),
                                                            what = c('mate_status','flag','mapq','qname','rname','strand','pos','seq')))
      eventChimericSortedReadsB <- sortBam(file = event.chimeric.B.bam,byQname = FALSE, destination = file.path(event.dir,'chimeric.B.sorted'))


      # # Save BAM files ----


      foo <- mergeBam(files = c(aligned.A.bam,aligned.B.bam,eventChimericSortedReadsA,eventChimericSortedReadsB),destination = event.bam,
                      indexDestination = TRUE,
                      overwrite = TRUE)

      # cleanup ----
      to_rm <- c(aligned.A.bam,aligned.B.bam,event.chimeric.A.bam,event.chimeric.B.bam,eventChimericSortedReadsB,eventChimericSortedReadsA)
      unlink(c(to_rm,paste0(to_rm,'.bai')))
   }
   fa.write.flag <-  !file.exists(all1.fq) | !file.exists(all2.fq) | overwrite
   if (fa.write.flag) {
      # fasta.cmd <- paste0('samtools collate -u -O ',event.bam,' | samtools fasta -1 ',all1.fa,' -2 ',all2.fa,' -0 /dev/null -s /dev/null - -n &> /dev/null')
      # system(fasta.cmd)


      fastq.cmd <- paste0('samtools collate -u -O ',event.bam,' | samtools fastq -1 ',all1.fq,' -2 ',all2.fq,' -0 /dev/null -s /dev/null - -n &> /dev/null')
      system(fastq.cmd)

      # chimeric reads only : event.chimeric.bamfile



      # fasta.cmd <- paste0('samtools collate -u -O ',event.chimeric.bamfile,' | samtools fasta -1 ',chi1.fa,' -2 ',chi2.fa,' -0 /dev/null -s /dev/null - -n')
      # system(fasta.cmd)
      #
      #
      # fastq.cmd <- paste0('samtools collate -u -O ',event.chimeric.bamfile,' | samtools fastq -1 ',chi1.fq,' -2 ',chi2.fq,' -0 /dev/null -s /dev/null - -n')
      # system(fastq.cmd)


   }

   # Write chimeric.bam file if necessary ----

   # Write jct info ----

   # Write event.bam ----

   # Write all1.fa, all2.fa files ----

   # Trimmomatic? ----
   run.trimmomatic <- FALSE
   if (run.trimmomatic){
      trim.summary.file <- file.path(raw.dir,'trimmomatic_summary.txt')
      trim.log.file <- file.path(raw.dir,'trimmomatic.log')
      out.basename <- file.path(raw.dir,'reads')
      trim.cmd <- paste('java -jar /usr/local/biotools/trimmomatic/0.39/trimmomatic-0.39.jar PE',
                        '-threads',CPU,'-summary',trim.summary.file,'-trimlog',trim.log.file,'-validatePairs',
                        '-baseout',out.basename,all1.fq,all2.fq)
      system(trim.cmd)
   }
   # Write local reference genome ----

   # Return created file names ----
}


eventwiseDeNovoAssemblyPrep <- function(seq.df,event,tx.dir,aligned.bam = NA,chimeric.bam = NA,
                                        overwrite = TRUE){
   # message(paste('Processing event',event,'files...'))
   event0 <- formatC(x = event,width =3,flag = 0,format = 'd')
   event.dir <- file.path(tx.dir,paste0('event_',event0))
   raw_event.dir <- file.path(event.dir,'raw_input')
   if (!dir.exists(raw_event.dir)){dir.create(raw_event.dir,recursive = TRUE)}


   # if (method == 'transrate'){flash.threshold = 10000000}
   span1.fq <- file.path(raw_event.dir,'span.1.fq')
   span2.fq <- file.path(raw_event.dir,'span.2.fq')
   write_spanning_frag_flag <- !file.exists(span1.fq) | !file.exists(span1.fq)
   if (overwrite | write_spanning_frag_flag){
      readsPP.out <- spanningRPsPreProcessing(seq.df = seq.df,event = event,
                                              event.tx.dir = raw_event.dir,
                                              #file.prefix = file.prefix,
                                              # flash.threshold = 1000000000,
                                              # remove.polyA = FALSE,
                                              # prehash_reads = FALSE,
                                              write.fq = TRUE)
   }

   jct.gr <- getGRangesFromInputDF(seq.df = seq.df,event = event)




   deNovoAssemblyPrep(event.dir = raw_event.dir,jct.gr = jct.gr,
                      aligned.bam = aligned.bam,
                      chimeric.bam = chimeric.bam,overwrite = overwrite)

   span.bam <- file.path(raw_event.dir,'span.bam')
   if (!file.exists(span.bam) | overwrite){
      blep <- system(paste('bwa mem',gen_ref_path,span1.fq,span2.fq,' 2> /dev/null | samtools view -b - | samtools sort - -o',span.bam,'&& samtools index',span.bam,'&> /dev/null'),
                     ignore.stdout = TRUE,ignore.stderr = TRUE,intern = FALSE)
   }

}

getGRangesFromInputDF <- function(seq.df = seq.df,event = event){
   jct.gr <- unique(GRanges(seqnames = c(paste0('chr',seq.df[seq.df$DNA_eventId == event,'chrA']),
                                         paste0('chr',seq.df[seq.df$DNA_eventId == event,'chrB'])),
                            ranges = IRanges(start = c(seq.df[seq.df$DNA_eventId == event,'posAmin'],
                                                       seq.df[seq.df$DNA_eventId == event,'posBmin']),
                                             end = c(seq.df[seq.df$DNA_eventId == event,'posAmax'],
                                                     seq.df[seq.df$DNA_eventId == event,'posBmax']))))

   return(jct.gr)
}

salmon <- function(tx.dir,salmon_input.fa,reads1,reads2,tpm_threshold = 30,overwrite = FALSE,numCores = 12){
   if (!file.exists(salmon_input.fa)){return(data.frame())}

   salmon.dir <- file.path(tx.dir,'salmon'); if (!dir.exists(salmon.dir)){dir.create(salmon.dir)}

   salmon_out <- file.path(salmon.dir,'quant.sf')
   if (!file.exists(salmon_out) | overwrite){
      input.file <- file.path(salmon.dir,'salmon.input.fa')
      index.file <- file.path(salmon.dir,'salmon.idx')
      salmon.log <- file.path(salmon.dir,'salmon.log')
      tx_ref_file <- '/research/labs/experpath/vasm/shared/Genome/Human/referenceFiles/transcriptome/gencode.v39.transcripts.fa'

      decoy.file <- '/research/labs/experpath/vasm/shared/Genome/Human/referenceFiles/transcriptome/gencode.v39.transcripts.decoys.txt'

      tx.cat.cmd <- paste('cat', salmon_input.fa,tx_ref_file,'>',input.file); system(tx.cat.cmd)

      # salmon index -t $SALMON_IN -d $TX_DECOY -i $SALMON_IDX -p 24 --gencode
      # salmon quant -l A -i $SALMON_IDX -1 $READS1 -2 $READS2 -o $SALMON_DIR -p 16
      tic()
      salmon.idx.cmd <- paste('conda run -n salmon_env salmon index --gencode -t', input.file, '-d',decoy.file,'-i',index.file,'-p',numCores,'>',salmon.log)
      foo <- system(salmon.idx.cmd,ignore.stderr = TRUE,intern = TRUE)

      salmon.quant.cmd <- paste('conda run -n salmon_env salmon quant --validateMappings -l A -i',index.file,'-1',reads1,'-2',reads2,'-p',numCores,'-o',salmon.dir)
      foo <- system(salmon.quant.cmd,intern=TRUE)
      salmon.dt <- toc()



   }


   quant.df <- read.table(file = salmon_out,header=TRUE)

   #quant_threshold = 30
   # quant_pass.df <- quant.df[quant.df$TPM > tpm_threshold,]
   # write.table(x = quant_pass.df,file = file.path(salmon.dir,'quant_pass.csv'),quote = FALSE,row.names = FALSE)
   return(quant.df)
}
