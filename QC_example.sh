#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate

mkdir -p /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc1/
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc1/
fastqc -o /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc1/ -f fastq /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0014_EKDN230004373-1A_HTV2CDSX5_L1_1.fq.gz /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0014_EKDN230004373-1A_HTV2CDSX5_L1_2.fq.gz > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc1/Z0014_EKDN230004373-1A_HTV2CDSX5_L1.log 2>&1 
mkdir -p /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/
trimmomatic PE -threads 5 -phred33 -trimlog /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/Z0014_EKDN230004373-1A_HTV2CDSX5_L1.txt /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0014_EKDN230004373-1A_HTV2CDSX5_L1_1.fq.gz /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0014_EKDN230004373-1A_HTV2CDSX5_L1_2.fq.gz -baseout /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/trim_Z0014_EKDN230004373-1A_HTV2CDSX5_L1.fq.gz ILLUMINACLIP:/prj/ycc-backup/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/trim_Z0014_EKDN230004373-1A_HTV2CDSX5_L1.log 2>&1 
rm /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0014_EKDN230004373-1A_HTV2CDSX5_L1_.*
mkdir -p /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc2/
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc2/
fastqc -o /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc2/ -f fastq /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/trim_Z0014_EKDN230004373-1A_HTV2CDSX5_L1_1P.fq.gz /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/trim/trim_Z0014_EKDN230004373-1A_HTV2CDSX5_L1_2P.fq.gz > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0014_EKDN230004373-1A_HTV2CDSX5_L1/qc2/Z0014_EKDN230004373-1A_HTV2CDSX5_L1.log 2>&1 
