# Downloading software

## conda 

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

## BWA

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

## Samtools

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

## sra-tools

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

## bcftools

https://github.com/samtools/bcftools

```bash
conda search -c bioconda -f bcftools
conda install -c bioconda bcftools=1.14
```

## bedtools 

https://github.com/arq5x/bedtools2

```bash
conda install bioconda::bedtools

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

## Trimmomatic

https://github.com/usadellab/Trimmomatic

```bash
conda install -c bioconda trimmomatic
```

## Picard

https://github.com/broadinstitute/picard

```bash
conda activate gatk
conda install -c bioconda picard
```

## GATK

https://github.com/broadinstitute/gatk

```bash
#git clone https://github.com/broadinstitute/gatk.git
#conda install -c bioconda java-jdk
cd /prj/ycc-backup/software/gatk-4.2.5.0
conda env create -f gatkcondaenv.yml
conda activate gatk

#On Denbi and LiDo
conda create -n gatk -c bioconda gatk4
```

## R

```bash
conda create -n r_env r-essentials r-base
conda activate r_env
```
