#!/bin/bash
#
#SBATCH --job-name analyse_paired_reads
#SBATCH --cpus-per-task=8
#SBATCH --mem=2GB
#SBATCH --account=rch-biol337
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err



# usage

# sbatch ../scripts/analyse_paired_reads data_directory

# data directory should contain all paired read files (.fastq.gz), an adapter file (adapters.fa),
# and a genome directory containing .gtf, .gff, and .fna files.



# Setup global environment

MODULEPATH=$MODULEPATH:/scratch/projects/sbs-apps/modulefiles
module load FastQC/0.12.1-Java-11
module load cutadapt/5.0
module load prefix/2025
module load BWA/0.7.18-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load htseq/htseq-2.1.2

READ_DATA_DIR=$1
GENOME_DATA_DIR="${READ_DATA_DIR}/genome"
ADAPTER_FILE="${READ_DATA_DIR}/adapters.fa"


# Functions

do_quality_check () {
	local READFILE=$1
	local OUTDIR=../results/00_fastqc/
	mkdir -p $OUTDIR
	printf "Doing fastqc quality check on ${READFILE} ...\n"
	fastqc $READFILE -o $OUTDIR > /dev/null
	if [[ $? -eq 0 ]]; then
		printf "Quality check results available at ${OUTDIR}\n"
	else
		printf "Failed to quality check ${READFILE}\n"
	fi
}


remove_adapters () {
	local READFILE1=$1
	local READFILE2=$2
	local OUTDIR=../results/01_cutadapt
	mkdir -p $OUTDIR
	local TEMP1="${OUTDIR}/$(basename $READFILE1 fastq.gz)_temp.fastq"
	local TEMP2="${OUTDIR}/$(basename $READFILE2 fastq.gz)_temp.fastq"
	local OUTFILE1="${OUTDIR}/$(basename $READFILE1 .gz)"
	local OUTFILE2="${OUTDIR}/$(basename $READFILE2 .gz)"

	printf "Removing adapter sequences using cutadapt ...\n"
 	cutadapt -a "file:$ADAPTER_FILE" --cores $SLURM_CPUS_PER_TASK -o $TEMP1 -p $TEMP2 $READFILE1 $READFILE2 > /dev/null
	cutadapt -q 15 -m 50 --cores $SLURM_CPUS_PER_TASK -o $OUTFILE1 -p $OUTFILE2 $TEMP1 $TEMP2 > /dev/null
	rm $TEMP1
	rm $TEMP2
	gzip -f $OUTFILE1
	gzip -f $OUTFILE2
	if [[ $? -eq 0 ]]; then
		printf "Adapter sequences removed from ${READFILE}\n"
	else
		printf "Failed to remove adapter sequences from ${READFILE}\n"
		exit 1
	fi
}


align () {
	local READFILE1="../results/01_cutadapt/$(basename $1)"
	local READFILE2="../results/01_cutadapt/$(basename $2)"

	# Get genome file ending in .fna extension
	local GENOMEFILE="$(find $GENOME_DATA_DIR -type f -name "*.fna" -print -quit)"

	# Index genome if required
	printf "Indexing genome file ${GENOMEFILE} using bwa index ...\n"
	local INDEXGENOMEFILE1="${GENOMEFILE}.amb"
	if [ -f $INDEXGENOMEFILE1 ]; then
		printf "Genome file already indexed\n"
	else
		bwa index $GENOMEFILE
		if [[ $? -eq 0 ]]; then
			printf "Genome file indexed\n"
		else
			printf "Failed to index genome file ${GENOMEFILE}\n"
			exit 1
		fi
	fi

	# Align reads against genome
	printf "Aligning using bwa mem ...\n"
	local OUTDIR=../results/02_bwa
	mkdir -p $OUTDIR
	local OUTFILE="${OUTDIR}/$(basename $READFILE1 _1.fastq.gz).bam"
	local TEMPFILE="${OUTDIR}/$(basename $READFILE1 _1.fastq.gz)_temp.bam"
	bwa mem -p -t $SLURM_CPUS_PER_TASK $GENOMEFILE $READFILE1 $READFILE2 > $TEMPFILE
	samtools sort $TEMPFILE > $OUTFILE
	rm $TEMPFILE
	samtools index $OUTFILE 
	if [[ $? -eq 0 ]]; then
		printf "Successfully aligned reads\n"
	else
		printf "Failed to align reads\n"
	fi
}


count () {
	local BASE_FILE_NAME="$(basename $1 _1.fastq.gz)"
	local ALIGN_FILE="../results/02_bwa/${BASE_FILE_NAME}.bam"
	local FEATURE_FILE="$(find $GENOME_DATA_DIR -type f -name "*.gtf" -print -quit)"
	mkdir -p ../results/03_htseq
	local OUTFILE="../results/03_htseq/${BASE_FILE_NAME}.tsv"

	printf "Counting reads using htseq-count ...\n"
	printf "DEBUG : ALIGN_FILE=${ALIGN_FILE}\n"
	printf "DEBUG : FEATURE_FILE=${FEATURE_FILE}\n"
	printf "DEBUG : OUTFILE=${OUTFILE}\n"
	htseq-count --type gene $ALIGN_FILE $FEATURE_FILE > $OUTFILE 
	if [[ $? -eq 0 ]]; then
		printf "Reads counted successfully\n"
	else
		printf "Failed to count reads\n"
		exit 1
	fi
}


# Main script
printf "\n\n============ RUNNING PAIRED READ ANALYSIS ============\n\n"
for READFILE1 in "${READ_DATA_DIR}"/*1.fastq.gz;
do
	#if [[ "$(basename $READFILE1)" == "ERR3140233_1.fastq.gz" ]]; then
	#	continue
	#fi

	READFILE2=${READFILE1::-10}2.fastq.gz
	printf "\nProcessing file pair: ${READFILE1} ${READFILE2}\n"

	# Generate quality check reports
	do_quality_check $READFILE1
	do_quality_check $READFILE2

	# Remove adapter sequences from both files
	remove_adapters $READFILE1 $READFILE2

	# Generate a sorted alignment file
	align $READFILE1 $READFILE2

	# Generate a .tsv counts file
	count $READFILE1

done

printf "\nReached end of script!\n\n\n"

# Run R analysis
