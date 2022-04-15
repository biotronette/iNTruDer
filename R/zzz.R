.onLoad <- function(libname,pkgname){

   velvet_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/velvet/1.2.10'
   oases_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/oases'
   flash_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/FLASH-1.2.11/flash'
   perl_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/perlbrew/perls/perl-5.34.0/bin/perl'
   vo_path <- '/research/labs/experpath/vasm/shared/NextGen/tools/VelvetOptimiser/2.2.6/VelvetOptimiser.pl'
   old_path <- Sys.getenv("PATH")
   if (!grepl('velvet',old_path)){Sys.setenv(PATH = paste(old_path, velvet_path, sep = ":"))}
   if (!grepl('perl_path',old_path)){Sys.setenv(PATH = paste(perl_path, old_path, sep = ":"))}
   rw.pause.time <- 0.01

   hg38 <- BSgenome.Hsapiens.UCSC.hg38

   edb <- EnsDb.Hsapiens.v86
   edb.ucsc <- edb
   ensembldb::seqlevelsStyle(edb.ucsc) <- 'UCSC'

   blat_path <- '/research/bsi/tools/biotools/blat/36.5/bin/blat'
   ooc_file <- '/research/labs/experpath/vasm/shared/Genome/Human/GRCh38/GCRh38.11.ooc'
   gen_ref_path <- '/research/labs/experpath/vasm/shared/Genome/Human/GRCh38/GCRh38_chromosomesOnly_noMask.fna'

   # source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/utilities.R", echo=FALSE)
   # source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/intruderPep.R", echo=FALSE)
   # source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/intruderAssembly.R", echo=FALSE)
   if (!exists('edb.ucsc')){
      edb <- EnsDb.Hsapiens.v86
      edb.ucsc <- edb
      ensembldb::seqlevelsStyle(edb.ucsc) <- 'UCSC'}
   # source("/research/labs/experpath/vasm/shared/NextGen/judell/Rprojects/vasm_bmd/julia/utilities.R", echo=FALSE)


}
