#!/bin/bash

printf "ChipSeq processing starts ..."


##  import parameter from config file ------------------------------------
	input_file=$1
	param_array_key=()
	eval param_array_key=($(awk -F+ '{printf($1 "\n")}' $input_file))
	
	param_array_value=()
	eval param_array_value=($(awk -F+ '{printf($2 "\n")}' $input_file))
	
	length=`expr ${#param_array_key[@]} - 1`
	for i in `seq 0 $length`
	do
	  case ${param_array_key[$i]} in
		WD)
		  WD=${param_array_value[$i]}  
		  ;;
		PRINSEQ)
		  PRINSEQ=${param_array_value[$i]}
		  ;;	
		BWA)
		  BWA=${param_array_value[$i]}
		  ;;
		PICARD)
		  PICARD=${param_array_value[$i]}
		  ;;
		TRIM)
		  TRIM=${param_array_value[$i]}
		  ;;
                QUALIMAP)
                  QUALIMAP=${param_array_value[$i]}
                  ;;
		FILE)
		  INPUT_FILE=${param_array_value[$i]}
		  ;;
		REF_GENOME)
		  REF_GENOME=${param_array_value[$i]}
		  ;;
		THREADS)
		  THREADS=${param_array_value[$i]}
		  ;;
		RG)
		  RG=${param_array_value[$i]}
		  ;;
		MEAN_QUAL)
		  MEAN_QUAL=${param_array_value[$i]}
		  ;;
		TRIM_CO)
		  TRIM_CO=${param_array_value[$i]}
		  ;;
		DIFF)
		  DIFF=${param_array_value[$i]}
		  ;;
                BAM_QUALI)
                  BAM_QUALI=${param_array_value[$i]}
                  ;;
                FASTQC)
                  FASTQC=${param_array_value[$i]}
                  ;;
                TMP_DIR)
                  TMP_DIR=${param_array_value[$i]}
                  ;;
	  esac
	done
#---------------------------------------------------------------------------------
GENOME=MOUSE

## change directory to working directory, folder with fastq raw data

cd $WD

## create data and log folder

LOG_F=${WD}/logs
RESULTS_F=${WD}/results/
mkdir $LOG_F
mkdir $RESULTS_F

# create main log file with important infos
LOG=${LOG_F}/log.txt
touch $LOG

# tmp log folder
TMP_LOG=${TMP_DIR}log.txt

## write config file settings to  log file
date 2>&1>>$LOG
cat $1 >>${LOG}

		  echo $WD >>${LOG}
		  echo $INPUT_FILE >>$LOG
		  echo $BWA >>$LOG
		  echo $PICARD >>$LOG
		  echo $PRINSEQ >>$LOG
		  echo $SAMTOOLS >>$LOG
		  echo $REF_GENOME >>$LOG
		  echo $RG >>$LOG
                  echo "Mapping quality filter: $BAM_QUALI" >>$LOG
                  echo  Trimming : $TRIM  >>$LOG
		  echo  Trimming cutoff: $TRIM_CO  >>$LOG
		  echo "Diff: $DIFF " >>$LOG
		  echo "Mean Qual: $MEAN_QUAL" >>$LOG

printf "\n" >> $LOG

##
set -x

## write tool versions in Log file
${SAMTOOLS}samtools --version &>>${LOG}
${PRINSEQ}prinseq-lite --version &>>${LOG}
${BWA}bwa | awk '1; NR == 5 { exit }' &>>${LOG}
printf "${PICARD} \n" &>>${LOG}
${FASTQC} --version &>>$LOG

printf "\n" >> $LOG

printf "The working directory is $WD"
########################################################################################################
##----------------  start processing ------------------------------------------------


### fastqc ------------------------------------------
if [ -e ${INPUT_FILE} ]; 
then
echo fastqc ...
        $FASTQC ${INPUT_FILE} -o ${LOG_F} -t ${THREADS} -d $TMP_DIR -f fastq 
fi

## filtering of raw data
## trimming if Trim parameter set to 1

touch $TMP_LOG

if [ "$TRIM" = "1" ];
then
        TRIM_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim
        if [ ! -e ${TRIM_OUT}.fastq ]
        then
                echo start trimming ...
                ${PRINSEQ}prinseq-lite -fastq $INPUT_FILE -out_good $TRIM_OUT -trim_to_len $TRIM_CO -log $TMP_LOG
        printf "\n" >> $LOG
        fi
  
else
 TRIM_OUT=${WD}${INPUT_FILE/.fastq}
fi

### fastqc ------------------------------------------
if [ -e ${TRIM_OUT}.fastq ]; 
then
echo fastqc ...
        $FASTQC ${TRIM_OUT}.fastq -o ${LOG_F} -t ${THREADS} -d $TMP_DIR -f fastq 
fi


## quality filtering ---------------------------------------
QUAL_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual
QUAL_OUT_BAD=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.bad
if [ ! -e ${QUAL_OUT}.fastq ];
then
        echo start quality filtering ...
        ${PRINSEQ}prinseq-lite -fastq ${TRIM_OUT}.fastq -out_good $QUAL_OUT -out_bad $QUAL_OUT_BAD -min_qual_mean $MEAN_QUAL -log $TMP_LOG
fi

cat $TMP_LOG >> ${LOG_F}/prinseq.log

### fastqc ------------------------------------------
if [  -e ${QUAL_OUT}.fastq ]; 
        then
        $FASTQC ${QUAL_OUT}.fastq -o ${LOG_F} -t ${THREADS} -d $TMP_DIR -f fastq 
fi


## alignment to reference genome ------------------------------
#-q INT quality threshold for read trimming down to 35bp [0], ist better not to use it
ALN_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.sai
if [ ! -e "${ALN_OUT}" ];
then
        echo start aligning reads ...
        ${BWA}bwa aln -n $DIFF -t $THREADS -f $ALN_OUT $REF_GENOME ${QUAL_OUT}.fastq &>>${LOG_F}/aln.log
fi

SAMSE_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.sam
if [ ! -e ${SAMSE_OUT} ];
then
        echo generate .sam ...
        ${BWA}bwa samse  -f $SAMSE_OUT -r $RG $REF_GENOME $ALN_OUT ${QUAL_OUT}.fastq &>>${LOG_F}/samse.log
fi

## convert sam to bam and bai and sort accoding to coordinate on ref genome -------------
SORTSAM_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.bam
if [ ! -e ${SORTSAM_OUT} ];
then
        echo sorting and binarizing reads 
java -jar -Xmx4g -Dsamjdk.try_use_intel_deflater=false ${PICARD}picard.jar SortSam SO=coordinate TMP_DIR=$TMP_DIR INPUT=${SAMSE_OUT} OUTPUT=$SORTSAM_OUT CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT &>>${LOG_F}/sortsam.log
fi

### qualimap -------------------------------------
OUTFILE=bamqc.sortSam.html
if [ ! -e ${LOG_F}/$OUTFILE ];
then
        echo generate bamqc statistic ...
        ${QUALIMAP} bamqc  --java-mem-size=12G -bam ${SORTSAM_OUT} -c -nt $THREADS -gd $GENOME -outformat HTML -outdir $LOG_F/SortSamQC -outfile $OUTFILE
fi

## remove reads with mapping quality below BAM_Quali and unmapped reads (-F 4) ------------------
MQ_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.uniqueM.bam
if [ ! -e $MQ_OUT ];
then
        echo start filtering mapped reads
${SAMTOOLS}samtools view -b -F 4 -q $BAM_QUALI $SORTSAM_OUT > $MQ_OUT
fi

# create index (not sure if necessary)
if [ ! -e ${MQ_OUT}.bai ];
then
        echo create index 
        ${SAMTOOLS}samtools index ${MQ_OUT}
fi

### qualimap -----------------------------------------------
OUTFILE=bamqc.mqOUT.html
if [ ! -e ${LOG_F}/$OUTFILE ];
then
        echo generate bamqc statistic ...
        ${QUALIMAP} bamqc  --java-mem-size=12G -bam ${MQ_OUT} -c -nt $THREADS -gd $GENOME -outformat HTML -outdir $LOG_F/mqOutQC -outfile $OUTFILEduplication level
fi

## mark duplicate reads ------------------------------------
DEDUP_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.uniqueM.markdups.bam
DEDUP_MET=${LOG_F}/${INPUT_FILE/.fastq}.trim.qual.uniqueM.markdups.metrics
if [ ! -e $DEDUP_OUT ];
then
        echo mark duplicate reads ...
java -jar -Xmx4g -Dsamjdk.try_use_intel_deflater=false ${PICARD}picard.jar MarkDuplicates INPUT=${MQ_OUT} OUTPUT=$DEDUP_OUT TMP_DIR=$TMP_DIR METRICS_FILE=$DEDUP_MET REMOVE_DUPLICATES=FALSE CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT &>>${LOG_F}/markDups.log
fi

### remove non standard chroms
STDCHROM_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.uniqueM.markdups.stdChrom.bam
${SAMTOOLS}samtools view -b -o $STDCHROM_OUT $DEDUP_OUT chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY

## remove duplicate reads -------------------------------------
DEDUP_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.uniqueM.rmdups.bam
DEDUP_MET=${LOG_F}/${INPUT_FILE/.fastq}.trim.qual.uniqueM.rmdups.metrics
if [ ! -e $DEDUP_OUT ];
then
        echo remove duplicate reads
java -jar -Xmx4g -Dsamjdk.try_use_intel_deflater=false ${PICARD}picard.jar MarkDuplicates INPUT=${MQ_OUT} OUTPUT=$DEDUP_OUT TMP_DIR=$TMP_DIR METRICS_FILE=$DEDUP_MET REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT &>>${LOG_F}/rmDups.log
fi

### qualimap ---------------------------------------------------
OUTFILE=bamqc.rmDups.html
if [ ! -e ${LOG_F}/$OUTFILE ];
then
        echo generate bamqc statistic ...
        ${QUALIMAP} bamqc  --java-mem-size=12G -bam ${DEDUP_OUT} -c -nt $THREADS -gd $GENOME -outformat HTML -outdir $LOG_F/rmDupsQC -outfile $OUTFILE
fi

### remove non standard chroms
STDCHROM_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.uniqueM.rmdups.stdChrom.bam
${SAMTOOLS}samtools view -b -o $STDCHROM_OUT $DEDUP_OUT chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY

### make index file
${SAMTOOLS}samtools index $STDCHROM_OUT

### sort bam file
SORT_OUT=${RESULTS_F}${INPUT_FILE/.fastq}.trim.qual.uniqueM.rmdups.stdChrom.sorted.bam
${SAMTOOLS}samtools sort -o $SORT_OUT $STDCHROM_OUT

### make index file
${SAMTOOLS}samtools index $SORT_OUT


### qualimap ---------------------------------------------------
OUTFILE=bamqc.rmNonStdChrom.html
if [ ! -e ${LOG_F}/$OUTFILE ];
then
        echo generate bamqc statistic ...
        ${QUALIMAP} bamqc  --java-mem-size=12G -bam ${SORT_OUT} -c -nt $THREADS -gd $GENOME -outformat HTML -outdir $LOG_F/rmNonStdChrom -outfile $OUTFILE
fi

#### remove big intermediate files sam files --------------------
#rm $ALN_OUT
#rm $SAMSE_OUT
#rm $INPUT_FILE
#if [ "$TRIM" = "1" ]
#        then
#        rm ${TRIM_OUT}.fastq
#fi
#rm ${QUAL_OUT}.fastq


#rm ${QUAL_OUT}.bai
#rm $MQ_OUT
#rm ${MQ_OUT}.bai
#-----------------------------------------------------------------
rm TMP_DIR=$TMP_LOG

echo $ALN_OUT
echo $SAMSE_OUT
echo $INPUT_FILE
if [ "$TRIM" = "1" ]
        then
        echo ${TRIM_OUT}.fastq
fi
echo ${QUAL_OUT}.fastq
echo ${QUAL_OUT}.bai
echo $MQ_OUT
echo ${MQ_OUT}.bai

date 2>&1>>$LOG

