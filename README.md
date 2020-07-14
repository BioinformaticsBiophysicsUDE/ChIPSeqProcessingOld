# ChIPSeqProcessingOld

## requirements

the script requires the following programs
bwa

prinseq-lite

fastqc

samtools

qualimap

picard tools

and an indexed mouse reference genome

## parameters
all adjustable parameters are set in the config file
WD + the working directory, directory that contains the input file, one fastq file
FILE + name of the input fastq file
THREADS + number og threads, some steps are multithreaded 
RG + Read Group
TRIM + 1 if the reads should be trimmed, else 0
TRIM_CO + trim cutoff, if trim = 1, the llength reads are be trimmed to
MEAN_QUAL + reads with an average Phread score below this threshold are removed
DIFF + used in bwa, aligned reads with a this difference to the reference are removed
PRINSEQ + path to the folder that contains prinseq-lite
BWA + path to the folder that contains bwa
PICARD + path to the folder that contains picard tools
REF_GENOME + <reference chromosome index file with fullpath>
SAMTOOLS + path to the folder that contains samtools
QUALIMAP + full path of qualimap
BAM_QUALI + aligend reads with a mapping quality below this treshold are removed
TMP + temporary file
FASTQC + full path of fastqc
TMP_DIR + path of a temporary directory, used by picard tools
  
## run the script
./chipSeqPreProcessingV2.sh example.config.txt
