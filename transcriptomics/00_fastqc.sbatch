#!/bin/bash
#
#SBATCH --job-name fastqc
#SBATCH --cpus-per-task=8
#SBATCH --mem=2GB
#SBATCH --account=rch-biol337
#SBATCH --time=1-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#usage
#sbatch ../scripts/01_fastqc.sbatch <fastqFile>

# add additional MODULEPATH location
MODULEPATH=$MODULEPATH:/scratch/projects/sbs-apps/modulefiles

# load modules
module load FastQC/0.12.1-Java-11



# run a script



OUTDIR=../results/fastqc
# below command makes an output folder for fastqc results if it does not already exist
mkdir -p $OUTDIR
# above command makes an output folder for fastqc results if it does not already exist


# The below command will analyse one fastq.gz file
#fastqc /scratch/projects/RCH-BIOL337/BIOL337_2025/Module_1_transcriptomics/DATA/ERR3140233_1.fastq.gz  -o $OUTDIR

# Like above but with variables
#DATADIR=/scratch/projects/RCH-BIOL337/BIOL337_2025/Module_1_transcriptomics/DATA/
#READFILE=ERR3140233_1.fastq.gz
#fastqc $DATADIR/$READFILE -o $OUTDIR

# Even better - iterate through all fastq.gz files
DATADIR=/scratch/projects/RCH-BIOL337/BIOL337_2025/Module_1_transcriptomics/DATA/
for READFILE in $DATADIR/*fastq.gz
do 
	fastqc $READFILE -o $OUTDIR
done
