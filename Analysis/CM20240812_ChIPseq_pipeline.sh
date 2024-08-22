# CM20240812
# User: Mathias, Christine
# Original Author: Tom van Schaik

# Process ChIP-seq
# This is a wrapper script with function calls to process the ChIP-seq data.
# The original wrapper script (/DATA/usr/m.eder/projects/ChIP_seq/me20220218_E1846_CTCF_F1mESCs/TvS_analysis/me20220222_ChIPseq_pipeline.sh) was written by Tom to align the ChIP-seq data 
# to an in silico genome with the transposon as an extra chromosome.

# Here I (Christine) adapt it to just align to mm10, for the cell line that does not contain any transposon.


# To do: explain screens
#   Initiate: screen -S name (ChipAnalysis)
#   Resume: screen -r name
#   List screens: screen -ls
#   Go back to main session: control+A+D


# To do: explain condas
#   Activate: conda activate name
#   Deactivate: conda deactivate


# I did not run step 2, 3 (fastqc and trimming), but copied them from the previous analysis, since they are independent of what genome we align to. 

#############################
### 1) Prepare output
cd /DATA/projects/Sox2/CTCF_ChIP_analysis

dir_chipseq="results"
mkdir $dir_chipseq

# List files
# samples_dir="/shared/gcf/m.eder/6764/fastq_files"
# samples_dir="data/"
# samples=$(ls $samples_dir/*.fastq.gz)
samples="/shared/gcf/m.eder/6764/fastq_files/6764_2_ME_E2_CGATGT_S2_R1_001.fastq.gz"

# In Unix to print the samples: `echo $samples`

#############################
### 2) FastQC reports
# conda env create -f "/DATA/scratch/usr/t.v.schaik/proj/tests/results/ts220218_chipseq_processing/conda_environments/fastqc.yaml"
conda activate fastqc
mkdir $dir_chipseq/fastqc
fastqc $samples -o $dir_chipseq/fastqc -t 12

# Reports show good reads, but with lots of over-represented sequences. I will
# need some adapter trimming

#############################
### 3) FastP adapter trimming
# conda env create -f /home/t.v.schaik/mydata/proj/tests/results/ts220218_chipseq_processing/conda_environments/fastp.yaml
conda activate fastp
mkdir $dir_chipseq/fastp

# Process them separately
# Note that this is single-end instead of the paired-end replication timing
for p1 in $samples; do
base=$(basename ${p1%.fa*})
# Run FastP
echo "processing $base"
fastp \
-i ${samples_dir}/${base}.fastq.gz \
-o $dir_chipseq/fastp/${base}.fastq.gz \
--html $dir_chipseq/fastp/${base}_fastp.html \
--json $dir_chipseq/fastp/${base}_fastp.json \
-w 12 &> $dir_chipseq/fastp/${base}.log
done

conda activate fastqc
mkdir $dir_chipseq/fastp_fastqc
fastqc $dir_chipseq/fastp/*fastq.gz -o $dir_chipseq/fastp_fastqc -t 12

#############################
### 4) BWA Alignment
conda activate /DATA/scratch/usr/t.v.schaik/miniconda3/envs/4DN_mapper
mkdir $dir_chipseq/mapping

#load genome index
genome_index="/DATA/scratch/usr/t.v.schaik/data/genomes/mm10/mm10"

# Process the samples - all together - single-end datacon
for p1 in $samples; do
base=$(basename ${p1%.fa*})
# Run FastP
echo "processing $base"
/DATA/scratch/usr/t.v.schaik/proj/3D_nucleus/results/ts200921_LaminaNucleolus_AdditionalExperiments/bin/mapping/mapping_bwa.sh \
-r $dir_chipseq/fastp/${base}.fastq.gz \
-i $genome_index \
-o $dir_chipseq/mapping \
-d -c 12 -f $dir_chipseq/mapping/${base}_bwa.log
done

#############################
### 5) MultiQC - combined quality report
# conda env create -f /DATA/scratch/usr/t.v.schaik/proj/tests/results/ts220218_chipseq_processing/conda_environments/multiqc.yaml
conda activate /DATA/scratch/usr/t.v.schaik/miniconda3/envs/multiqc
multiqc $dir_chipseq -f -o $dir_chipseq

#############################
### 6) Make BigWigs with bamCoverage (deeptools suite)
# conda env create -f /DATA/scratch/usr/t.v.schaik/proj/tests/results/ts220218_chipseq_processing/conda_environments/deeptools.yaml
conda activate /DATA/scratch/usr/t.v.schaik/miniconda3/envs/deeptools
mkdir $dir_chipseq/bigwig

bam_files=$(ls $dir_chipseq/mapping/ | grep bam | grep -v bai)

for bam in $bam_files; do
base=$(basename ${bam%.bam})
# Run deeptools - bamCoverage
bamCoverage \
-b $dir_chipseq/mapping/$bam \
-o $dir_chipseq/bigwig/$base.bw \
--binSize 5 \
--numberOfProcessors 12 \
--effectiveGenomeSize 2652783500 \
--normalizeUsing RPKM \
--ignoreDuplicates \
--minMappingQuality 50
done

#############################
### 7) Call peaks with MACS2
# conda env create -f /DATA/scratch/usr/t.v.schaik/proj/tests/results/ts220218_chipseq_processing/conda_environments/macs2.yaml
conda activate /DATA/scratch/usr/t.v.schaik/miniconda3/envs/macs2
mkdir $dir_chipseq/peaks

bam_files="$(ls $dir_chipseq/mapping/ | grep bam | grep -v bai)"

for bam in $bam_files; do
base=$(basename ${bam%.bam})
# Run macs2 callpeak
macs2 callpeak \
-t $dir_chipseq/mapping/$bam \
--name $base \
-g 2652783500 \
--outdir $dir_chipseq/peaks 2> \
$dir_chipseq/peaks/${base}.log
done
