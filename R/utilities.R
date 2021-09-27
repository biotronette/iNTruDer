cohortAnalysis <- function(out.dir.list,case.prefix.list){
   n = length(out.dir.list)
   log.df <- data.frame('Case' = case.prefix.list,
                        'Directory' = out.dir.list,
                        'RNAsupportRunDate' = character(length = n),
                        'hasSpanningRNAfile' = logical(length = n),
                        'nSpanningRPs' = integer(length = n),
                        'nEvents' = integer(length = n),
                        'assemblyDate' = character(length = n),
                        'hasTX' = logical(length = n),
                        'nTX' = integer(length = n),
                        'map_rate' = numeric(length = n),
                        'ref_rate' = numeric(length = n),
                        'neopepRunDate' = character(length = n),
                        'nORFs' = integer(length = n),
                        'nNeopeps' = integer(length = n))
   
   for (i_case in 1:n){
      reads_file <- file.path(out.dir.list[i_case],grep('supportForDNAseqJcts_reads',list.files(out.dir.list[i_case]),value = TRUE))
      log.df[i_case,'hasSpanningRNAfile'] <- any(file.exists(reads_file))
      if (any(file.exists(reads_file))){
         
         if (length(reads_file) > 1){
            n_SR <- integer(length(reads_file))
            for (j_rf in 1:length(reads_file)){
               n_SR[j_rf] <- length(readLines(reads_file[j_rf])) - grep('DNA_eventId',readLines(reads_file[j_rf]))
               
            }
            reads_file <- reads_file[n_SR == max(n_SR)][1]
            
            
         }
         log.df[i_case,'RNAsupportRunDate'] <- format(file.mtime(reads_file))
         log.df[i_case,'nSpanningRPs'] <- length(readLines(reads_file)) - grep('DNA_eventId',readLines(reads_file))
         log.df[i_case,'nEvents'] <- ifelse(log.df[i_case,'nSpanningRPs'] > 0,
                                            length(unique(sapply(strsplit(readLines(reads_file)[(grep('DNA_eventId',readLines(reads_file)) + 1):length(readLines(reads_file))],','),FUN='[',1))),
                                            0)
      }
      tx_dir <- file.path(out.dir.list[i_case],'tx')
      tx_file <- file.path(out.dir.list[i_case],'tx',grep('TX_out.csv',list.files(tx_dir),value = TRUE))
      
      if ( length(tx_file) > 0 && any(file.exists(tx_file) & length(readLines(tx_file)) > 1)){
         log.df[i_case,'hasTX'] <- length(readLines(tx_file)) > 1
         log.df[i_case,'assemblyDate'] <- format(file.mtime(tx_file))
         #tx <- data.frame()
         tx = read.csv(tx_file,row.names = NULL)
         log.df[i_case,'nTX'] <- nrow(tx)
         # log.df[i_case,'map_rate'] <- sum(tx$blat) / nrow(tx)
      }
      
      blat_file <- file.path(out.dir.list[i_case],'blat',grep('junction_tx.rds',list.files(file.path(out.dir.list[i_case],'blat')),value = TRUE))
      if (length(blat_file) > 0){
         blat.gr <- readRDS(blat_file)
         log.df[i_case,'map_rate'] <- sum(tx$tx_num %in% blat.gr$query) / nrow(tx)
      }
      if (log.df[i_case,'nEvents'] > 0){
         ref_files <- file.path(out.dir.list[i_case],'prot/ref',list.files(file.path(out.dir.list[i_case],'prot/ref')))
         log.df[i_case,'ref_rate'] <- length(ref_files) / length(unique(tx$DNA_eventId))
      }
      
         event_file <- file.path(out.dir.list[i_case],grep('event_summary',list.files(out.dir.list[i_case]),value = TRUE))
      if (length(event_file) > 0 && any(file.size(event_file) > 1)){
         event <- read.csv(event_file)
         
         log.df[i_case,'ref_rate'] <- sum(event$nRefProt > 0) / nrow(event)
        
         log.df[i_case,'nNeopeps'] <- sum(event$nNeopep)
      }
      np_file <- file.path(out.dir.list[i_case],'prot',grep('neopeptides.fa',list.files(file.path(out.dir.list[i_case],'prot')),value=TRUE))
      if (length(np_file) > 0){
         orf_file <- file.path(out.dir.list[i_case],'prot',grep('orfs.csv',list.files(file.path(out.dir.list[i_case],'prot')),value=TRUE))
         
         log.df[i_case,'nORFs'] <- length(readLines(orf_file)) - 1
         log.df[i_case,'nNeopeps'] <- sum(grepl('>',readLines(np_file)))
         log.df[i_case,'neopepRunDate'] <- format(file.mtime(np_file))
      }
   }
   
   return(log.df)
}

captainsLog <- function(sid,eid,tid,msg,
                        log.file = '/research/labs/experpath/vasm/shared/NextGen/judell/intruder.log'){
   newLine <- paste(sid,eid,tid,format.Date(Sys.time()),msg,sep='\t')
   fileConn<-file(log.file)
   writeLines(newLine, fileConn)
   close(fileConn)
}