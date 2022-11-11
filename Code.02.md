# Epigenetic marks on hybrid species

## Tools

https://epigenie.com/epigenetic-tools-and-databases/

https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

## Mapping sequences & variant calling

### Downloading software

#### conda 

https://docs.conda.io/en/latest/

```bash
# Download the lastest version
# Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Conda activate
source ~/miniconda3/etc/profile.d/conda.sh
conda activate

vim ~/.bashrc
export LS_OPTIONS='--color=auto'
eval "$(dircolors -b)"
alias ls='ls $LS_OPTIONS'

source ~/.bashrc
#vim ~/.profile
#source ~/.bashrc

# Setup channel
conda config --add channels defaults 
conda config --add channels bioconda 
conda config --add channels conda-forge
```

#### Homebrew

https://brew.sh/

```bash
# On CeBiTec, must download on ThinLinc interface
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
echo "eval \"\$($(brew --prefix)/bin/brew shellenv)\"" >>~/.profile
```

#### BWA

https://github.com/bwa-mem2/bwa-mem2

Working dir: /opt/bwa/

```bash
# Setup directory
sudo chmod -R 777 /opt/
mkdir -p /opt/bwa/
cd /opt/bwa/
# Download
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 | sudo tar jxf -
# Test-1
bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index ref.fa
    ##Error message
    Please verify that both the operating system and the processor support Intel(R) X87, CMOV, MMX, 	FXSAVE, SSE, SSE2, SSE3, SSSE3, SSE4_1, SSE4_2, MOVBE, POPCNT, F16C, AVX, FMA, BMI, LZCNT, AVX2, AVX512DQ, AVX512F, ADX, AVX512CD, AVX512BW, AVX512VL and CLWB instructions.

# Test-2
/opt/bwa/bwa-mem2-2.0pre2_x64-linux/bwa-mem2.sse41 index ref.fa

# On CeBiTec
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 | tar jxf -
# Test-1
/prj/ycc-backup/software/bwa-mem2-2.0pre2_x64-linux/bwa-mem2
```

#### Samtools

https://github.com/samtools/samtools

```bash
#conda install -c bioconda samtools
#conda install -c conda-forge -c bioconda samtools bzip2
conda search -c bioconda -f samtools

samtools                        1.10      h9402c20_2  bioconda
samtools                        1.11      h6270b1f_0  bioconda
samtools                        1.12      h9aed4be_1  bioconda
samtools                        1.12      hd5e65b6_0  bioconda
samtools                        1.13      h8c37831_0  bioconda
samtools                        1.14      hb421002_0  bioconda

conda create -n samtools -c bioconda samtools=1.13
	# version
```

#### sra-tools

https://github.com/ncbi/sra-tools

https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

```bash
# Denbi
sh setup-apt.sh
source /etc/profile.d/sra-tools.sh
	#export PATH=$HOME/tools/sra-tools/sratoolkit.2.9.6-ubuntu64/bin:$PATH
# CebiTec
# Download pre-built binaries to /prj/ycc-backup/software
```

#### bcftools

https://github.com/samtools/bcftools

```bash
conda search -c bioconda -f bcftools
conda install -c bioconda bcftools=1.14
```

#### bedtools 

https://github.com/arq5x/bedtools2

```bash
# On CeBiTec, must be installed on ThinLinc interface
cd /prj/ycc-backup/software
# Upload the latest release
tar zxvf bedtools-2.30.0.tar.gz
cd /prj/ycc-backup/software/bedtools
make
# Alternative way: home brew
brew install gcc
brew install bedtools
```

#### Trimmomatic

https://github.com/usadellab/Trimmomatic

```bash
conda install -c bioconda trimmomatic
```

#### AutoTrim

https://github.com/schellt/autotrim

```bash

```

#### Picard

https://github.com/broadinstitute/picard

```bash
conda activate gatk
conda install -c bioconda picard
```

#### GATK

https://github.com/broadinstitute/gatk

```bash
#git clone https://github.com/broadinstitute/gatk.git
#conda install -c bioconda java-jdk
cd /prj/ycc-backup/software/gatk-4.2.5.0
conda env create -f gatkcondaenv.yml
conda activate gatk

#On Denbi
conda create -n gatk4 -c bioconda gatk4
```

#### R

```bash
conda create -n r_env r-essentials r-base
conda activate r_env


```



### Download great tit and related species genomes

**On CeBiTec**

File dir: /prj/ycc-backup/ParusMA/pma.01.01/

```bash
# Great tit
cd /prj/ycc-backup/ParusMA/pma.01.01/
wget -r -nH --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.3_Parus_major1.1/

# Blue tit
cd /prj/ycc-backup/ParusMA/pma.01.01/
wget -r -nH --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/901/205/GCF_002901205.1_cyaCae2/

# Ground tit
cd /prj/ycc-backup/ParusMA/pma.01.01/
wget -r -nH --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/331/425/GCF_000331425.1_PseHum1.0/

# Collared flycatcher
cd /prj/ycc-backup/ParusMA/pma.01.01/
wget -r -nH --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/247/815/GCF_000247815.1_FicAlb1.5/

# Zebra finch
cd /prj/ycc-backup/ParusMA/pma.01.01/
wget -r -nH --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/
```

The project id of other great tit species

Sample list: D:\Dropbox\VHF274\ParusMA_hybrid\pma\pma.01.01\Sample_list.xlsx



### Test bwa with the highest coverage sample

**On CeBiTec (ThinLinc)** 

Generic data pre-processing pipeline: https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl

gatk: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651

https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

https://github.com/gencorefacility/variant-calling-pipeline-gatk4

https://github.com/nickhir/GermlineMutationCalling

https://evodify.com/gatk-in-non-model-organism/ (old version)

File dir: /prj/ycc-backup/ParusMA/pma.01.01/

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/

Test with SRX2714718 (run: SRR5423297)



#### Download sequence file

File dir: /prj/ycc-backup/ParusMA/pma.01.01/

```bash
# Setup directory
mkdir -p /prj/ycc-backup/ParusMA/pma.01.01/sra

# Download sra file
cd /prj/ycc-backup/ParusMA/pma.01.01/sra
prefetch --verbose --force all SRR5423297 --max-size 50G

#####Future reference#####
# Get all srr id from project id
esearch -db biosample -query SRX2714718 | efetch -mode xml
esearch -db sra -query PRJNA208335 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR
#esearch -db sra -query PRJNA208335 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump -X 10 --split-files

# Extract SRR/ERR:
grep -E 'SRR|ERR' XDR_169_ids.txt > downloads.txt

# Find SRAs from SRS:
grep 'SRS' XDR_169_ids.txt | parallel "esearch -db sra -query {} | efetch --format runinfo | cut -d ',' -f 1 | grep SRR" >> downloads.txt

# Now make sure there are no duplicates, then download using GNU parallel to have 4 (or as many your disk can handle) streams in parallel:
sort -u downloads.txt | parallel -j 4 "prefetch {}"
#####Future reference#####

# Transfer to fastq
fasterq-dump -v -p --split-3 /vol/storage/ParusMA/pma.01.01/sra/SRR5423297/SRR5423297.sra 

```



#### QC check

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/fastqc

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/trimmed

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/fastqc_trim

```bash
. ~/.bashrc

# QC check with fastqc
mkdir /prj/ycc-backup/ParusMA/pma.01.02/fastqc
cd /prj/ycc-backup/ParusMA/pma.01.02/fastqc

nohup fastqc -o /prj/ycc-backup/ParusMA/pma.01.02/fastqc -f fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR5423297_1.fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR5423297_2.fastq &

fastqc -o /prj/ycc-backup/ParusMA/pma.01.02/fastqc -f fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR1793430_1.fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR1793430_2.fastq

# Trim seq with Trimmomatic
mkdir /prj/ycc-backup/ParusMA/pma.01.02/trimmed
cd /prj/ycc-backup/ParusMA/pma.01.02/trimmed

nohup trimmomatic PE -threads 28 -phred33 -trimlog SRR5423297.txt /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR5423297_1.fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR5423297_2.fastq -baseout trim_SRR5423297.fastq.gz ILLUMINACLIP:/prj/ycc-backup/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &

nohup trimmomatic PE -threads 28 -phred33 -trimlog SRR1793430.txt /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR1793430_1.fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR1793430_2.fastq -baseout trim_SRR1793430.fastq.gz ILLUMINACLIP:/prj/ycc-backup/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
    # Trimmomatic usage
    # ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
    # phred33: based on the encoding score in the Illumina platform
    # TruSeq3: as used by HiSeq and MiSeq machines
    # ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
    # Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
    # Remove leading low quality or N bases (below quality 3) (LEADING:3)
    # Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
    # Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
    # Drop reads below the 36 bases long (MINLEN:36)
    # SLIDINGWINDOW:4:15, perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. By considering multiple bases, a single poor quality base will not cause the removal of high quality data later in the read.
    
##### if use qsub #####
vim /prj/ycc-backup/ParusMA/pma.01.02/fastqc/trim.sh

#!/bin/bash
trimmomatic PE -threads 28 -phred33 -trimlog SRR5423297.txt /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR5423297_1.fastq /prj/ycc-backup/ParusMA/pma.01.01/sra/SRR5423297_2.fastq -baseout trim_SRR5423297.fastq.gz ILLUMINACLIP:/prj/ycc-backup/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

nohup qsub -P fair_share -l vf=1G -pe multislot 28 -cwd /prj/ycc-backup/ParusMA/pma.01.02/fastqc/trim.sh &
	#cwd: use the current directory
# Check status
qstat -pe multislot
qstat -f | grep chnyuch
qstat -f -q 16647127
qdel -f 16647127
	#https://bioinformatics.mdc-berlin.de/intro2UnixandSGE/sun_grid_engine_for_beginners/how_to_submit_a_job_using_qsub.html
	
# qsub example
qsub -pe multi-slot 5 myjobscript.sh
qsub -P fair_share -pe multi-slot 5 myjobscript.sh
qsub -P fair_share -pe -o out.txt -e error.txt multi-slot 5 -cwd myjobscript.sh
qsub -P fair_share -pe multi-slot 16 -l vf=4G -o out.txt -e error.txt  -cwd myjobscript.sh
qsub -P fair_share -l vf=64G -o out.txt -e error.txt -cwd myjobscript.sh
##### use qsub end #####

# QC check after trimming with fastqc
mkdir /prj/ycc-backup/ParusMA/pma.01.02/fastqc_trim
cd /prj/ycc-backup/ParusMA/pma.01.02/fastqc_trim

nohup fastqc -o /prj/ycc-backup/ParusMA/pma.01.02/fastqc_trim -f fastq /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR5423297_1P.fastq.gz /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR5423297_2P.fastq.gz &

nohup fastqc -o /prj/ycc-backup/ParusMA/pma.01.02/fastqc_trim/ -d /prj/ycc-backup/ParusMA/pma.01.02/fastqc_trim/tmp -f fastq /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR1793430_1P.fastq.gz /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR1793430_2P.fastq.gz &
	#Error-1
	#OpenJDK 64-Bit Server VM warning: Insufficient space for shared memory file:
    #739885
	#Try using the -Djava.io.tmpdir= option to select an alternate temp location.
	#Solution-1
	# -d: tmp dir 
```



#### Mapping

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/mapping

```bash
# Index reference fa
mkdir -p /prj/ycc-backup/ParusMA/pma.01.02/ref
cd /prj/ycc-backup/ParusMA/pma.01.02/ref
cp /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna.gz /prj/ycc-backup/ParusMA/pma.01.02/ref/

gzip -d /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna.gz

/prj/ycc-backup/software/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index GCF_001522545.3_Parus_major1.1_genomic.fna

# Mapping using BWA
mkdir -p /prj/ycc-backup/ParusMA/pma.01.02/mapping

/prj/ycc-backup/software/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 26 -R '@RG\tID:SRR5423297\tSM:SRR5423297\tLB:sample01\tPL:ILLUMINA' -M /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR5423297_1P.fastq.gz /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR5423297_2P.fastq.gz > /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297.sam 
	# Parameter setting
	# -M: mark shorter split hits as secondary
	# -Y: use soft clipping for supplementary alignments
	# By default, BWA-MEM uses soft clipping for the primary alignment and hard clipping for supplementary alignments (complementary to primary alignment).
	# Soft-clipped: bases in 5' and 3' of the read are NOT part of the alignment; hard clipped: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file. The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
	# https://www.biostars.org/p/109333/
	# https://www.biostars.org/p/119537/
	# Read group setting
	# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
	# In fastq header: @SRR1793435.1 HWI-ST1307:84:C1W60ACXX:3:1101:1374:2173
	# instrument name: HWI-ST1307
	# run id: 84
	# flowcell id: C1W60ACXX
	# flowcell lane: 3
	# tile number within the flowcell lane: 1101
	# 'x'-coordinate of the cluster within the tile: 1374
	# 'y'-coordinate of the cluster within the tile: 2173
	# https://www.gdc-docs.ethz.ch/MDA/site/getdata/
	
/prj/ycc-backup/software/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 10 -R '@RG\tID:SRR1793430\tSM:SRR1793430\tLB:sample02\tPL:ILLUMINA' -M /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR1793430_1P.fastq.gz /prj/ycc-backup/ParusMA/pma.01.02/trimmed/trim_SRR1793430_2P.fastq.gz > /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430.sam 

#####For MarkDuplicats but not MarkDuplicatesSpark#####
# Conversion to bam file
conda activate samtools
nohup samtools sort -@ 24 -o /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297.bam /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297.sam &
	# Error-1
	samtools: error while loading shared libraries: libncurses.so.5: cannot open shared object file: No such file or directory
	# Solve-1-1 (failed)
	conda install -c conda-forge ncurses
	# Solve-1-2 (cause error-2)
	# Reinstall samtools
	# Error-2
	samtools: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
	# Solve-2-1
	# Reinstall samtools

nohup samtools sort -@ 24 -o /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430.bam /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430.sam &

# index sample
samtools index /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297.bam
samtools index /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430.bam 
#####MarkDuplicats but not MarkDuplicatesSpark End#####

# index reference
samtools faidx /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna

# Back to base conda env
conda deactivate

```



#### Realign around putative indels

https://github.com/broadinstitute/gatk-docs/blob/master/blog-2012-to-2019/2016-06-21-Changing_workflows_around_calling_SNPs_and_indels.md?id=7847



#### MarkDuplicates

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/mapping

```bash
conda activate gatk
chmod 700 /prj/ycc-backup/software/gatk-4.2.5.0/gatk
mkdir -p /prj/ycc-backup/ParusMA/pma.01.02/mapping/tmp/

nohup /prj/ycc-backup/software/gatk-4.2.5.0/gatk --java-options "-Xmx60G -Xms60G -Djava.io.tmpdir=/prj/ycc-backup/ParusMA/pma.01.02/mapping/tmp/" MarkDuplicatesSpark -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297.bam -M /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_dedup_metrics.txt -O /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam --tmp-dir /prj/ycc-backup/ParusMA/pma.01.02/mapping/tmp/ &
	#Error-1 java too new
	#Solution-1
	# conda install -c conda-forge openjdk=8.0.312
	#Error-2 No space left on device
	#22/02/22 17:22:05 ERROR Executor: Exception in task 10986.0 in stage 15.0 (TID 150502)
	#java.io.IOException: No space left on device
	#Solution-2
	# --tmp-dir 
	# --java-options "-XmxNG -XmsMG -Djava.io.tmpdir=/path/to/tmpdir"
	# https://gatk.broadinstitute.org/hc/en-us/community/posts/360067258451-MarkDuplicatesSpark-doesn-t-work-for-large-bam-files

# Sort bam file
#samtools sort -@ 20 -n -O bam -o /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted.bam /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297.bam
	#Sort reads by name

##### clear ineffective flag from BWA #####
#samtools fixmate -@ 20 -n -O bam -o fixmate_252A_01_L001.bam bwa_s_252A_01_L001.bam 
## sort the bam file by coordinate
#samtools sort -@ 20 -O bam -o fixmate_s_252A_01_L001.bam fixmate_252A_01_L001.bam 
##### clear end #####

#nohup /prj/ycc-backup/software/gatk-4.2.5.0/gatk MarkDuplicates -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted.bam -M /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_dedup_metrics.txt -O /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam &

/prj/ycc-backup/software/gatk-4.2.5.0/gatk --java-options "-Xmx60G -Xms60G -Djava.io.tmpdir=/prj/ycc-backup/ParusMA/pma.01.02/mapping/tmp/" MarkDuplicatesSpark -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430.sam -M /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_dedup_metrics.txt -O /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam --tmp-dir /prj/ycc-backup/ParusMA/pma.01.02/mapping/tmp/

# Remove sam file
rm /prj/ycc-backup/ParusMA/pma.01.02/mapping/*.sam

conda deactivate
```



#### Collect Alignment & Insert Size Metrics

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/mapping

```bash
conda activate gatk

#Alignment metrics
picard CollectAlignmentSummaryMetrics R=/prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna I=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam O=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_alignment_metrics.txt

picard CollectAlignmentSummaryMetrics R=/prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna I=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam O=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_alignment_metrics.txt

#Insert metrics
picard CollectInsertSizeMetrics INPUT=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam OUTPUT=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_insert_metrics.txt HISTOGRAM_FILE=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_insert_size_histogram.pdf

picard CollectInsertSizeMetrics INPUT=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam OUTPUT=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_insert_metrics.txt HISTOGRAM_FILE=/prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_insert_size_histogram.pdf

conda deactivate

conda activate samtools
samtools depth -a /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam > /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_depth_out.txt

samtools depth -a /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam > /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_depth_out.txt

```



#### Call Variants and extract SNP and INDEL #1

Working dir: /prj/ycc-backup/ParusMA/pma.01.02/callVariant

```bash
conda activate gatk
mkdir -p /prj/ycc-backup/ParusMA/pma.01.02/callVariant

# Create dict file for reference
/prj/ycc-backup/software/gatk-4.2.5.0/gatk CreateSequenceDictionary -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna

# Call variants
/prj/ycc-backup/software/gatk-4.2.5.0/gatk HaplotypeCaller -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_variants.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk HaplotypeCaller -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_variants.vcf

# SNPs
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_variants.vcf --select-type-to-include SNP -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_snps.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_variants.vcf --select-type-to-include SNP -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_snps.vcf

# INDELS
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_variants.vcf --select-type-to-include INDEL -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_indels.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_variants.vcf --select-type-to-include INDEL -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_indels.vcf

```



#### Filter and exclude SNPs & Indels

!!!check each parameters mean

https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

```bash
# Mark SNP and indels that should be filter out
# Filter snp (default setting)
/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_snps.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

    #Warning message
    #09:35:25.744 WARN  JexlEngine - ![0,14]: 'ReadPosRankSum < -8.0;' undefined variable ReadPosRankSum
    #09:35:25.744 WARN  JexlEngine - ![0,9]: 'MQRankSum < -12.5;' undefined variable MQRankSum

/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_snps.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Filter indel (default setting)
/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

# Exclude filter SNPs and Indels
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants --exclude-filtered -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_filtered_snps.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_bqsr_snps.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants --exclude-filtered -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_filtered_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_bqsr_indels.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants --exclude-filtered -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_filtered_snps.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_bqsr_snps.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants --exclude-filtered --exclude-filtered -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_filtered_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_bqsr_indels.vcf

```



#### Base Quality Score Recalibration (BQSR) 

```bash
# BQSR #1 
/prj/ycc-backup/software/gatk-4.2.5.0/gatk BaseRecalibrator -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_bqsr_snps.vcf --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_bqsr_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recal_data.table

/prj/ycc-backup/software/gatk-4.2.5.0/gatk BaseRecalibrator -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_bqsr_snps.vcf --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_bqsr_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recal_data.table

# Apply BQSR #1
# *_recal_reads.bam is analysis-ready
/prj/ycc-backup/software/gatk-4.2.5.0/gatk ApplyBQSR -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR5423297_sorted_dedup.bam -bqsr /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recal_data.table -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recal_reads.bam

/prj/ycc-backup/software/gatk-4.2.5.0/gatk ApplyBQSR -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/mapping/SRR1793430_sorted_dedup.bam -bqsr /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recal_data.table -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recal_reads.bam

# BQSR #2 for report
/prj/ycc-backup/software/gatk-4.2.5.0/gatk BaseRecalibrator -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recal_reads.bam --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_bqsr_snps.vcf --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_bqsr_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_post_recal_data.table

/prj/ycc-backup/software/gatk-4.2.5.0/gatk BaseRecalibrator -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recal_reads.bam --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_bqsr_snps.vcf --known-sites /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_bqsr_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_post_recal_data.table

# Analyze Covariates
/prj/ycc-backup/software/gatk-4.2.5.0/gatk AnalyzeCovariates -before /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recal_data.table -after /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_post_recal_data.table -plots /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recalibration_plots.pdf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk AnalyzeCovariates -before /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recal_data.table -after /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_post_recal_data.table -plots /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recalibration_plots.pdf

```



#### Call Variants and extract SNP and INDEL #2

```bash
# Call variants (Genotype likelihoods)
/prj/ycc-backup/software/gatk-4.2.5.0/gatk HaplotypeCaller -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_recal_reads.bam -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_variants_recal.vcf -ERC GVCF

/prj/ycc-backup/software/gatk-4.2.5.0/gatk HaplotypeCaller -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_recal_reads.bam -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_variants_recal.vcf -ERC GVCF

# Extract SNPs and Indels
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_variants_recal.vcf --select-type-to-include SNP -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_snps_recal.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_variants_recal.vcf --select-type-to-include INDEL -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_indels_recal.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_variants_recal.vcf --select-type-to-include SNP -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_snps_recal.vcf

/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_variants_recal.vcf --select-type-to-include INDEL -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_indels_recal.vcf
```



#### Filter and exclude SNPs & Indels

!!!check each parameters mean

```bash
/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_snps_recal.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_filtered_snps_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_raw_indels_recal.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR5423297_filtered_indels_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"

/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_snps_recal.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_filtered_snps_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_raw_indels_recal.vcf -O /prj/ycc-backup/ParusMA/pma.01.02/callVariant/SRR1793430_filtered_indels_final.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0"
```



#### Annotate SNPs and compile statistics

```bash
# combine the files in mapping and callvariant folder
ln -snf /prj/ycc-backup/ParusMA/pma.01.02/mapping/* /prj/ycc-backup/ParusMA/pma.01.02/callVariant/

# annotation
java -jar snpEff.jar -v <snpeff_db> filtered_snps_final.vcf > $filtered_snps_final.ann.vcf

# compile statistic
/prj/ycc-backup/script/parse_metrics.sh SRR5423297 > SRR5423297_report.csv
/prj/ycc-backup/script/parse_metrics.sh SRR1793430 > SRR1793430_report.csv
```



#### Joint genotyping

https://github.com/nickhir/GermlineMutationCalling

https://hpc.nih.gov/training/gatk_tutorial/

```bash
# Prepare inputs
# Generate sample map
# set up directory
mkdir -p /prj/ycc-backup/ParusMA/pma.01.02/vqsr
# Generate chromosome list
grep '>' /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna | sed -e 's/ .*$//g; s/>//g' > ../interval.2.list

/prj/ycc-backup/software/gatk-4.2.5.0/gatk --java-options "-Djava.io.tmpdir=/prj/ycc-backup/ParusMA/pma.01.02/ -Xms10G -Xmx10G -XX:ParallelGCThreads=2" GenomicsDBImport --genomicsdb-workspace-path /prj/ycc-backup/ParusMA/pma.01.02/pma_db -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna --batch-size 2 --sample-name-map /prj/ycc-backup/ParusMA/pma.01.02/test.map --tmp-dir /prj/ycc-backup/ParusMA/pma.01.02/ --max-num-intervals-to-import-in-parallel 2 -L /prj/ycc-backup/ParusMA/pma.01.02/interval.2.list 
	# --genomicsdb-workspace-path must be an new or empty folder

# Joint genotyping
/prj/ycc-backup/software/gatk-4.2.5.0/gatk --java-options "-Djava.io.tmpdir=/prj/ycc-backup/ParusMA/pma.01.02/ -Xms50G -Xmx50G -XX:ParallelGCThreads=1" GenotypeGVCFs -R /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V gendb://NC_031768.1_gdb -O NC_031768.1.vcf.gz
gatk GenotypeGVCFs

# MergeCohortVCFs
gatk MergeVcfs
gatk MakeSitesOnlyVcf

```



#### VQSR





### Pre-analysis of Parus major genomes

**On CeBiTec (ThinLinc)** 

File dir: /prj/ycc-backup/ParusMA/pma.01.01/

Working dir: /prj/ycc-backup/ParusMA/pma.01.03/

script: /prj/ycc-backup/script/

#### Download SRA files

```bash
# Date: 2022/03/01
# script: /prj/ycc-backup/script/get_sra.01.py
python /prj/ycc-backup/script/get_sra.01.py

#Execute download shell
#chmod -R 755 /prj/ycc-backup/ParusMA/pma.01.01/sra/
#nohup sh /prj/ycc-backup/ParusMA/pma.01.01/sra/sh/dwn_sra_exe.sh &

# Execute fasterq shell
# nohup sh /prj/ycc-backup/ParusMA/pma.01.01/sra/sh/fasterq_exe.sh &
```



#### QC check

```bash
# script: /prj/ycc-backup/script/get_sra.01.py
# output dir: /prj/ycc-backup/ParusMA/pma.01.03/fastqc/out/
conda activate
python /prj/ycc-backup/script/readsQC.01.py

# Execute shell script separately
# nohup sh /prj/ycc-backup/ParusMA/pma.01.03/fastqc/sh/qc1_exe.sh &
# nohup sh /prj/ycc-backup/ParusMA/pma.01.03/trim/sh/trim_exe.sh &
# nohup sh /prj/ycc-backup/ParusMA/pma.01.03/fastqc_trim/sh/qc2_exe.sh &
```

##### Result

Per base sequence quality improved, but some still failed in Per base sequence content and/or Per sequence GC content and/or Per tile sequence quality.

https://www.biostars.org/p/365748/



#### Mapping 

```bash
# script: /prj/ycc-backup/script/mapping.01.py
# output dir: /prj/ycc-backup/ParusMA/pma.01.03/ready_reads/
conda activate samtoools
python /prj/ycc-backup/script/mapping.01.py

# Execute shell script separately
# nohup sh /prj/ycc-backup/ParusMA/pma.01.03/mapping/sh/bwa_exe.sh &
```



#### MarkDuplicate and BQSR

```bash
# script: /prj/ycc-backup/script/bqsr.01.py
# output dir: /prj/ycc-backup/ParusMA/pma.01.03/ready_reads/
conda activate gatk
ulimit -n 9000
	# https://gatk.broadinstitute.org/hc/en-us/community/posts/360076989831-MarkDuplicatesSpark-running-but-not-sorting-creating-deduped-bam-files
python /prj/ycc-backup/script/bqsr.04.py

###qsub usage
# https://www.cebitec.uni-bielefeld.de/intranet/en/computer-systems/unix/using-rcinfo
# https://www.cebitec.uni-bielefeld.de/intranet/en/computer-systems/compute-cluster/running-batch-jobs
# https://docs.hpc.shef.ac.uk/en/latest/parallel/JobArray.html


#  user configurable .profile for bash
#
#  define the packages you want to include, e.g.
#  RCINFO_ILIST="A B-1.1"
#
RCINFO_ILIST=""

#
#  define the packages you want to exclude, e.g.
#  RCINFO_XLIST="B-1.0"
#
RCINFO_XLIST=""
 
export RCINFO_ILIST RCINFO_XLIST
 
#
#  call rcinfo to create the environment
#
if [ -x /vol/local/bin/rcinfo ]; then 
   eval "`/vol/local/bin/rcinfo bash`"
fi

qsub -P fair_share -l vf=64G -o /prj/ycc-backup/ParusMA/pma.01.03/ready_reads/sh/gatk/job_gatk0.out -e /prj/ycc-backup/ParusMA/pma.01.03/ready_reads/sh/gatk/job_gatk0.err -cwd /prj/ycc-backup/ParusMA/pma.01.03/ready_reads/sh/gatk/job_gatk0.sh

qstat -pe multislot
qstat -f | grep chnyuch
qstat -f -q 17454825
qdel -f 17457922
qselect -U chnyuch | xargs qdel
qdel -f -s s -u chnyuch

# Summmary result
cat /prj/ycc-backup/ParusMA/pma.01.03/ready_reads/*/*.csv | sort -u > /prj/ycc-backup/ParusMA/pma.01.03/bqsr_stat.csv

# delelte subfolder
cd /prj/ycc-backup/ParusMA/pma.01.03/ready_reads
find . -maxdepth 2 -type d -name *_tmp -exec rm -r {} + 
find . -maxdepth 2 -type f -name *_report.log -exec rm {} + 
```



#### Sanity check

##### ANGSD

Dir: /prj/ycc-backup/ParusMA/pma.01.05

http://www.popgen.dk/angsd/index.php/ANGSD

https://github.com/ANGSD/angsd

http://www.popgen.dk/software/index.php/NgsAdmix

Filtering

https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial2_genotype_snp_calling/markdowns/snp_calling.md

Angsd is for calculating allele frequencies. 

```bash
-doMaf  0 (Calculate persite frequencies '.mafs.gz')
        1: Frequency (fixed major and minor)
        2: Frequency (fixed major unknown minor)
        4: Frequency from genotype probabilities
        8: AlleleCounts based method (known major minor)
        NB. Filedumping is supressed if value is negative
        -doMaf 7 (1+2+4) will use the first three estimators. If the allele frequencies are estimated from the genotype likelihoods then you need to infer the major and minor allele (-doMajorMinor)
-doPost 0 (Calculate posterior prob 3xgprob)
        1: Using frequency as prior
        2: Using uniform prior
        3: Using SFS as prior (still in development)
Filters:
        -minMaf         -1.000000       (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
        -rmTriallelic   0.000000        (Remove sites with a pvalue lower)
Extras:
        -ref    (null)  (Filename for fasta reference)
        -anc    (null)  (Filename for fasta ancestral)
        -eps    0.001000 [Only used for -doMaf &8]
        -beagleProb     0 (Dump beagle style postprobs)
        -indFname       (null) (file containing individual inbreedcoeficients)
NB These frequency estimators requires major/minor -doMajorMinor

-GL 0 
	1: SAMtools
	2: GATK
	3: SOAPsnp
	4: SYK
	-trim		0		(zero means no trimming)
	-tmpdir		angsd_tmpdir/	(used by SOAPsnp)
	-errors		(null)		(used by SYK)
	-minInd		0		(0 indicates no filtering)

Filedumping:
	-doGlf	0
	1: binary glf (10 log likes)	.glf.gz
	2: beagle likelihood file	.beagle.gz
	3: binary 3 times likelihood	.glf.gz
	4: text version (10 log likes)	.glf.gz
```



##### Prepare inputs

```bash
# ANGSD install
conda create -n angsd -c bioconda python=3.9 angsd

# pcangsd install
conda activate angsd
cd /prj/ycc-backup/software/
git clone https://github.com/Rosemeis/pcangsd.git
cd /prj/ycc-backup/software/pcangsd/
pip install --user -r requirements.txt
python setup.py build_ext --inplace
pip3 install -e .

# make chromosome list
mkdir -p /prj/ycc-backup/ParusMA/pma.01.05/
grep '>' /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna | sed -e 's/,.*$//g; s/>//g; s/P.*Abel //g; s/ /\t/1; s/ /_/g' | (head -n 32 && tail -n 1) > /prj/ycc-backup/ParusMA/pma.01.05/interval.lst

# make bam file list 
find /prj/ycc-backup/ParusMA/pma.01.03/ready_reads -name '*_recal_reads.bam' | sort  > /prj/ycc-backup/ParusMA/pma.01.05/bam_file.lst

# label list
sed -e 's/\/prj\/ycc-backup\/ParusMA\/pma.01.03\/ready_reads\///g; s/\/.*//g' /prj/ycc-backup/ParusMA/pma.01.05/bam_file.lst > /prj/ycc-backup/ParusMA/pma.01.05/label.lst

```



##### PCA 

```bash
# PCA
# whole genome
# Use pcangsd now, maybe we can use PCA_MDS (single base approach)
mkdir -p /prj/ycc-backup/ParusMA/pma.01.05/pca.1

# Calculate allele frequencies
# shell file: /prj/ycc-backup/ParusMA/pma.01.05/pca.2/pma_pca.1.sh
#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate angsd
angsd -GL 2 -out /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam /prj/ycc-backup/ParusMA/pma.01.05/bam_file.lst
	# doMajorMinor 1: estimate allele frequencies from genotype likelihoods
	# GL 2: GATK likelihood
	# doGlf 2: beagle genotype likelihood format
	# doMaf 1: Frequency (fixed major and minor)
	# SNP_pval: filtering out the sites that are very likely to be polymorphic with a p-value less than 10^-6

qsub -P fair_share -l idle=1 -l vf=4G -pe multislot 10 -o /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca_sh.1.out -e /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca_sh.1.err -cwd /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.sh
	# Your job 17501798 ("pma_pca.1.sh") has been submitted
	
# PCA plot
cd /prj/ycc-backup/ParusMA/pma.01.05/pca.1/
pcangsd -b /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out.beagle.gz -o /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.plot.out -t 10
# Plot in R
# R script: pca_plot.1.R

# PCA per chromosome
python script: angsd_pca.01.py
mkdir -p /prj/ycc-backup/ParusMA/pma.01.05/pca.2/cov_out
ln -snf /prj/ycc-backup/ParusMA/pma.01.05/pca.2/*/*.cov /prj/ycc-backup/ParusMA/pma.01.05/pca.2/cov_out/

# Plot in R
# R script: pca_plot.2.R
```



##### Admixture

https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html

```bash
# Allele frequencies were calculated in PCA session
conda activate angsd
mkdir -p /prj/ycc-backup/ParusMA/pma.01.05/admix
cd /prj/ycc-backup/ParusMA/pma.01.05/
# shell file: /prj/ycc-backup/ParusMA/pma.01.05/admix/admix.1.sh
#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate angsd
NGSadmix -likes /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out.beagle.gz -K 3 -P 10 -o /prj/ycc-backup/ParusMA/pma.01.05/admix/pma_wg_admix3 -minMaf 0.05 
NGSadmix -likes /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out.beagle.gz -K 2 -P 10 -o /prj/ycc-backup/ParusMA/pma.01.05/admix/pma_wg_admix2 -minMaf 0.05 
NGSadmix -likes /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out.beagle.gz -K 1 -P 10 -o /prj/ycc-backup/ParusMA/pma.01.05/admix/pma_wg_admix1 -minMaf 0.05 
NGSadmix -likes /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out.beagle.gz -K 4 -P 10 -o /prj/ycc-backup/ParusMA/pma.01.05/admix/pma_wg_admix4 -minMaf 0.05
NGSadmix -likes /prj/ycc-backup/ParusMA/pma.01.05/pca.1/pma_pca.1.out.beagle.gz -K 5 -P 10 -o /prj/ycc-backup/ParusMA/pma.01.05/admix/pma_wg_admix5 -minMaf 0.05

qsub -P fair_share -l idle=1 -l vf=30G -pe multislot 10 -o /prj/ycc-backup/ParusMA/pma.01.05/admix/admix_sh_1.out -e /prj/ycc-backup/ParusMA/pma.01.05/admix/admix_sh_1.err -cwd /prj/ycc-backup/ParusMA/pma.01.05/admix/admix.1.sh
	#Your job 17509236 ("admix.1.sh") has been submitted

#Python srcipt: angsd_admx.01.py

# Plot in R 
# R script: angsd_admx.01.R
           1            2            3            4            5 
         Inf 1.403503e+09 5.146001e+01 5.910322e+01 8.683183e+02
	# Choose K = 2
	# Plot K=2 and K=3
```



##### Call heterozygosity, Tajima's D

```bash
# Call heterozygosity (only chromosome)
# Filter out low quality sites due to 
mkdir -p /prj/ycc-backup/ParusMA/pma.01.05/htz/
# shell file: /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.1.sh
#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate angsd
angsd -bam /prj/ycc-backup/ParusMA/pma.01.05/bam_file.lst -ref /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -out /prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1 -anc /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -rf /prj/ycc-backup/ParusMA/pma.01.05/interval.lst -dosaf 1 -GL 1 -P 25 -C 50 -minQ 20 -minmapq 20
	# remove low quality bases: -C 50 -ref ref.fa -minQ 20 -minmapq 30
	# -minMapQ    0   Discard reads with mapping quality below
    # -minQ       13  Discard bases with base quality below 
    # -C          0   adjust mapQ for excessive mismatches (as SAMtools), supply -ref
	
qsub -P fair_share -l idle=1 -l vf=20G -pe multislot 25 -o /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz.1.sh.out -e /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz.1.sh.err -cwd /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.1.sh
	# Your job 17571281 ("pma.1.htz.sh") has been submitted

# realSFDS
# shell file: /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.2.sh
#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate angsd
realSFS /prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1.saf.idx -fold 1 -P 25 > /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.1.sfs

qsub -P fair_share -l idle=1 -l vf=20G -pe multislot 25 -o /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz.2.sh.out -e /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz.2.sh.err -cwd /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.2.sh

# realSFS limited nSites
<<'###BLOCK-COMMENT'
# Change to per chromosome
# shell file: /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.2.sh
#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate angsd
realSFS /prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1.saf.idx -fold 1 -P 10 -nSites 800000000 > /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.1.sfs
	# Error-1
	-> Version of fname:/prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1.saf.idx is:2
	-> Assuming .saf.gz file: /prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1.saf.gz
	-> Assuming .saf.pos.gz: /prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1.saf.pos.gz
	-> args: tole:0.000000 nthreads:10 maxiter:100 nsites:0 start:(null) chr:(null) start:-1 stop:-1 fstout:(null) oldout:0 seed:-1 bootstrap:0 resample_chr:0 whichFst:0 fold:1 ref:(null) anc:(null)
	-> generating offset remapper lookup
	-> Looks like you will allocate too much memory, consider starting the program with a lower -nSites argument
	-> nSites: 984857046
	-> The choice of -nSites will require atleast: 604865.937500 megabyte memory, that is at least: 117.30% of total memory
terminate called after throwing an instance of 'std::bad_alloc'
  what():  std::bad_alloc
/var/spool/codine/laxbach/job_scripts/17571275: line 4: 1190555 Aborted                 (core dumped) realSFS /prj/ycc-backup/ParusMA/pma.01.05/htz/htz.1.saf.idx -fold 1 -P 10 > /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.1.sfs
	# Default nSites: 984857047
	# -nSites 800000000, failed due to further setting needed

# qsub -P fair_share -l idle=1 -l vf=30G -l mem_free=50G,s_vmem=50G,h_vmem=70G -pe multislot 10 -o /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz_sh.2.out -e /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz_sh.2.err -cwd /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.2.sh

qsub -P fair_share -l idle=1 -l vf=30G -pe multislot 10 -o /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz_sh.2.out -e /prj/ycc-backup/ParusMA/pma.01.05/htz/pma_htz_sh.2.err -cwd /prj/ycc-backup/ParusMA/pma.01.05/htz/pma.htz.2.sh
###BLOCK-COMMENT

# realSFS per chromosome
# Fail due to the same output in each chromosome
python scrtipt: angsd_sfs.01.py

# plot prep
mkdir -p /prj/ycc-backup/ParusMA/pma.01.05/htz/sfs
ln -snf /prj/ycc-backup/ParusMA/pma.01.05/htz/*/*.sfs /prj/ycc-backup/ParusMA/pma.01.05/htz/sfs/
cd /prj/ycc-backup/ParusMA/pma.01.05/htz/sfs

#### Plot in R ####
s<-scan('out.sfs')
s<-s[-c(1,length(s))]
s<-s/sum(s)
barplot(s,names=1:length(s),main='SFS')

R script: angsd_sfs.01.R

# Calculate theta
python script: angsd_sfs.01.py
realSFS saf2theta htz.1.saf.idx -sfs htz.1.sfs -outname htz.1.sfs.out -fold 1
# print out theata
conda activate angsd
cd /prj/ycc-backup/ParusMA/pma.01.05/htz
for i in /prj/ycc-backup/ParusMA/pma.01.05/htz/*.out.thetas.idx; do thetaStat print $i 2>/dev/null 

thetaStat print out.thetas.idx 2>/dev/null |head

# Estimate Tajima's D
#thetaStat do_stat out.thetas.idx
#cat out.thetas.idx.pestPG
# Estimate stats with sliding window
thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
```



## Great tit reference genome

### Software installation

#### gffread

https://github.com/gpertea/gffread

```bash
conda install -c bioconda gffread
```



#### EMBOSS CpG island prediction

https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Sequence+Statistics

##### EMBOSS Cpgplot

https://emboss.bioinformatics.nl/cgi-bin/emboss/help/cpgplot

##### EMBOSS Cpgreport

https://www.bioinformatics.nl/cgi-bin/emboss/help/newcpgreport

```bash
#cd /prj/ycc-backup/software/
#wget ftp://ftp.ebi.ac.uk/pub/software/unix/EMBOSS/EMBOSS-6.6.0.tar.gz 
#tar zxvf /prj/ycc-backup/software/EMBOSS-6.6.0.tar.gz
#cd /prj/ycc-backup/software/EMBOSS-6.6.0
#./configure
#make
#make install

conda install -c bioconda emboss
```

##### Bedops

```bash
conda install -c bioconda bedops
```



### Check gff file problem

Working dir: /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1

```bash
# Uncompress fna and gff files
gzip -kd /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.gff.gz

gzip -kd /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna.gz

# Extract cds test
conda activate
gffread -g /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna -x ref_cds.fa -w ref_rna.fa /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.gff

# Brief gff summary 
cat GCF_001522545.3_Parus_major1.1_genomic.gff | grep -v "^#" | cut -f3 | sort -n | uniq -c
```



### Create bed file for CpG island

Working dir: /prj/ycc-backup/ParusMA/pma.01.03/ref/

```bash
conda activate
# Gernerate CpG island report
newcpgreport /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -window 100 -shift 1 -minlen 200 -minoe 0.6 -minpc 50. -outfile /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.cpgreport

cpgplot /prj/ycc-backup/ParusMA/pma.01.02/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -window 100 -minlen 200 -minoe 0.6 -minpc 50. -outfile /prj/ycc-backup/ParusMA/pma.01.03/cpg/GCF_001522545.3_Parus_major1.1_genomic.cpgplot -graph png -outfeat /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.cpgplot.gff

# Transfrom to bed file
grep -v '#' /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.cpgplot.gff | awk -F'\t' '{print $1 "\t" $4 "\t" $5}' > /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.cpgplot.bed

# TpG

# TpA

# Gene body
cat /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.gff | grep -v "^#" | cut -f3 | sort -n | uniq -c
 531222 CDS
      3 C_gene_segment
      8 V_gene_segment
   1539 cDNA_match
 610836 exon
  19015 gene
     15 guide_RNA
   5298 lnc_RNA
  38606 mRNA
     46 pseudogene
      3 rRNA
   1675 region
      1 sequence_feature
     32 snRNA
    171 snoRNA
    214 tRNA
   2078 transcript

gff2bed < /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.gff | cut -f1-6 > /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.bed

grep 'gene' /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.bed > /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic_gbd.bed
```



## 3- and 5- way whole genome alignment

https://genome.cshlp.org/content/early/2020/03/19/gr.255752.119

### Software installation

#### MULTIZ

```bash
conda install -c bioconda multiz
```

#### LASTZ

```bash
conda install -c bioconda lastz
```

#### axtchain

```bash
conda install -c bioconda ucsc-axtchain
```

#### faidx

```bash
pip install pyfaidx
```







### Sample preparation

```bash
# Uncompress gz file
gzip -kd /prj/ycc-backup/ParusMA/pma.01.01/GCF_000247815.1_FicAlb1.5/GCF_000247815.1_FicAlb1.5_genomic.fna.gz
gzip -kd /prj/ycc-backup/ParusMA/pma.01.01/GCF_000331425.1_PseHum1.0/GCF_000331425.1_PseHum1.0_genomic.fna.gz
gzip -kd /prj/ycc-backup/ParusMA/pma.01.01/GCF_002901205.1_cyaCae2/GCF_002901205.1_cyaCae2_genomic.fna.gz
gzip -kd /prj/ycc-backup/ParusMA/pma.01.01/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz
```



### Whole genome alignment

Workind dir: /prj/ycc-backup/ParusMA/pma.01.04/

```bash
# name switching table 
grep '>NC' /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna | sed -e 's/,.*//g; s/Parus major isolate Abel //g; s/ /\t/1; s/>//g' | sed -e 's/ /_/g' > /prj/ycc-backup/ParusMA/pma.01.04/pma_scaff_lst

python /prj/ycc-backup/script/multiw_aln.01.py


# lastz alignment
# Great tit vs zebra finch
lastz /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna[multiple] /prj/ycc-backup/ParusMA/pma.01.01/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna --ambiguous=iupac --allocate:traceback=1.99G --notransition --step=20 --format=maf ‑‑rdotplot+score=/prj/ycc-backup/ParusMA/pma.01.04/lastz/pma_tgu_rplot_1.txt > /prj/ycc-backup/ParusMA/pma.01.04/lastz/pma_tgu_1.maf

lastz /prj/ycc-backup/ParusMA/pma.01.01/GCF_001522545.3_Parus_major1.1/GCF_001522545.3_Parus_major1.1_genomic.fna[multiple] /prj/ycc-backup/ParusMA/pma.01.01/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna --ambiguous=iupac --allocate:traceback=1.99G --notransition --step=20 --format=blastn ‑‑rdotplot+score=/prj/ycc-backup/ParusMA/pma.01.04/lastz/pma_tgu_rplot_1.txt > /prj/ycc-backup/ParusMA/pma.01.04/lastz/pma_tgu_1.blastn.txt


```



We used a whole genome alignment approach to obtain divergence estimates for the UTRs. We downloaded the reference genomes for chicken (v5.0; [Hillier et al. 2004](javascript:;)), zebra finch (v3.2.4; [Warren et al. 2010](javascript:;)) and great tit (v1.0.4; [Laine et al. 2016](javascript:;)). First we created pairwise alignments with the zebra finch as a reference using LASTZ ([Harris 2007](javascript:;)), following the procedures described in previous analyses of avian genomes ([Jarvis et al. 2014](javascript:;); [Zhang et al. 2014](javascript:;)). This was followed by chaining and netting using axtChain and chainNet, respectively ([Kent et al. 2003](javascript:;)). Finally, single coverage was ensured for the reference genome using single_cov2.v11 from the MULTIZ package and the pairwise alignments were aligned with MULTIZ ([Blanchette et al. 2004](javascript:;)). Coordinates for 5′ UTRs and 3′ UTRs in the zebra finch and great tit genomes were obtained from their respective annotation databases, and were analyzed together in each species. Only the UTRs of genes in the orthologous gene set described above were analyzed. This resulted in UTR alignments for 4,524 genes.
