#!/bin/bash
#
#SBATCH --job-name process_single_read_rrbs_files
#SBATCH --cpus-per-task=8
#SBATCH --mem=2GB
#SBATCH --account=rch-biol337
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err



# usage

# sbatch ../scripts/process_single_read_rrbs_files data_directory

# data directory should contain:
# all single read illumina-generated files (.fastq.gz),
# an adapter file (adapters.fa)


# Setup global environment

MODULEPATH=$MODULEPATH:/scratch/projects/sbs-apps/modulefiles
module load FastQC/0.12.1-Java-11
module load cutadapt/5.0
module load prefix/latest
module load bismark/0.25.1
module load SAMtools/1.21-GCC-13.3.0

READ_DATA_DIR=$1
ADAPTER_FILE="${READ_DATA_DIR}/adapters.fa"
GENOME_FOLDER="${READ_DATA_DIR}/genome"

JOB_START="$(date +'%H-%M-%S')"
ERROR_LOG="job_${JOB_START}_error_log.txt"
printf "$(date +"%Y-%m-%d %H-%M-%S") - ERROR LOG START\n\n" > $ERROR_LOG

# Functions

do_quality_check () {
	local READFILE=$1
	local FILE_ID="$(basename $READFILE .fastq.gz)"
	local OUTDIR=../results/00_fastqc/
	mkdir -p $OUTDIR
	printf "Starting fastqc quality check on ${FILE_ID} ...\n"
	fastqc $READFILE -o $OUTDIR > /dev/null
	if [[ $? -eq 0 ]]; then
		printf "${FILE_ID} quality check complete\n"
	else
		printf "Failed to quality check ${FILE_ID}\n"
		printf "$(date +"%Y-%m-%d %H-%M-%S") - FAIL: QUALITY CHECK ${FILE_ID}\n" > $ERROR_LOG
	fi
}


remove_adapters () {
	local READFILE=$1
	local FILE_ID="$(basename $READFILE .fastq.gz)"
	local OUTDIR=../results/01_cutadapt
	mkdir -p $OUTDIR
	local TEMP="${OUTDIR}/${FILE_ID}_temp.fastq"
	local OUTFILE="${OUTDIR}/${FILE_ID}.fastq.gz"

	printf "Removing adapter sequences from ${FILE_ID} using cutadapt ...\n"
 	cutadapt -a "file:$ADAPTER_FILE" --cores $SLURM_CPUS_PER_TASK -o $TEMP $READFILE > /dev/null
	cutadapt -q 15 -m 50 --cores $SLURM_CPUS_PER_TASK -o $OUTFILE $TEMP > /dev/null
	rm $TEMP
	if [[ $? -eq 0 ]]; then
		printf "Adapter sequences removed from ${FILE_ID}\n"
	else
		printf "Failed to remove adapter sequences from ${FILE_ID}\n"
		printf "$(date +"%Y-%m-%d %H-%M-%S") - FAIL: ADAPTER REMOVAL ${FILE_ID}\n" > $ERROR_LOG
	fi
}


align () {
	local READFILE=$1
	local FILE_ID="$(basename ${READFILE} .fastq.gz)"
	local OUTDIR="../results/02_bismark"
	local TEMPDIR="${OUTDIR}/temp"
	mkdir -p $OUTDIR
	mkdir -p $TEMPDIR

	printf "Aligning ${FILE_ID} using bismark ...\n"
	bismark -o $OUTDIR -p $SLURM_CPUS_PER_TASK --temp_dir $TEMPDIR --parallel $SLURM_CPUS_PER_TASK $GENOME_FOLDER $READFILE > /dev/null
	if [[ $? -eq 0 ]]; then
		printf "Aligned ${FILE_ID} successfully\n"
	else
		printf "Failed to align ${FILE_ID} using bismark\n"
		printf "$(date +"%Y-%m-%d %H-%M-%S") - FAIL: ALIGNMENT ${FILE_ID}\n" > $ERROR_LOG
	fi
}



generate_coverage_file () {
	local READFILE=$1
	local FILE_ID="$(basename ${READFILE} .bam)"
	local OUTDIR="../results/03_bismark_coverage"
	mkdir -p $OUTDIR

	printf "Generating ${FILE_ID} coverage file ...\n"
	bismark_methylation_extractor -s -o $OUTDIR --bedGraph --ucsc --cytosine_report --parallel $(($SLURM_CPUS_PER_TASK/3)) --genome_folder $GENOME_FOLDER $READFILE > /dev/null
	if [[ $? -eq 0 ]]; then
		printf "Generated ${FILE_ID} coverage file successfully\n"
	else
		printf "Failed to generate ${FILE_ID} coverage report\n"
		printf "$(date +"%Y-%m-%d %H-%M-%S") - FAIL: COVERAGE REPORT ${FILE_ID}\n" > $ERROR_LOG
	fi
}



##### Main script -----------------------------------------------------------------------------------


printf "\n\n============ PROCESSING SINGLE READ RRBS FILES ============\n\n"


printf "$(date +"%Y-%m-%d %H-%M-%S") - Cleaning RRBS files ...\n\n"
#for READFILE in "${READ_DATA_DIR}"/*.fastq.gz;
#do
	# Generate quality check report
#	do_quality_check $READFILE &

	# Remove adapter sequences from file
#	remove_adapters $READFILE &

#done
#wait
printf "$(date +"%Y-%m-%d %H-%M-%S") - RRBS files cleaned\n\n"


# Check if genome is already indexed - don't need to repeat if already indexed
printf "$(date +"%Y-%m-%d %H-%M-%S") - Indexing genome file in ${GENOME_FOLDER} ...\n"
#INDEXGENOMEFILE1="${GENOME_FOLDER}/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2"
#if [ -f $INDEXGENOMEFILE1 ]; then
#	printf "$(date +"%Y-%m-%d %H-%M-%S") - Genome already indexed\n"
#else
#	bismark_genome_preparation $GENOME_FOLDER > /dev/null
#	if [[ $? -eq 0 ]]; then
#		printf "$(date +"%Y-%m-%d %H-%M-%S") - Successfully indexed genome\n\n"
#	else
#		printf "$(date +"%Y-%m-%d %H-%M-%S") - Failed to index genome\n"
#		printf "$(date +"%Y-%m-%d %H-%M-%S") - FAIL: GENOME INDEX $(basename ${READFILE} .fastq.gz)\n\n" > $ERROR_LOG
#	fi
#fi


# Run bismark on all adapterless files
#printf "$(date +"%Y-%m-%d %H-%M-%S") - Aligning adapterless files ...\n\n"
#for READFILE in "../results/01_cutadapt"/*fastq.gz;
#do
#	align $READFILE &
#done
#wait
#printf "\n\n$(date +"%Y-%m-%d %H-%M-%S") - Completed file alignments\n\n"



printf "$(date +"%Y-%m-%d %H-%M-%S") - Generating coverage reports ...\n\n"
for READFILE in "../results/02_bismark"/*.bam;
do
	generate_coverage_file $READFILE &
done
wait
printf "\n\n$(date +"%Y-%m-%d %H-%M-%S") - Generated coverage reports\n\n"



printf "$(date +"%Y-%m-%d %H-%M-%S") - ERROR LOG END\n\n" > $ERROR_LOG
printf "\n############# END OF SCRIPT ##############\n\n\n"
