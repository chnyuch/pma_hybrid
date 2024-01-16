#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate
mkdir -p /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc1/
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc1/
fastqc -o /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc1/ -f fastq /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_1.fq.gz /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_2.fq.gz > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc1/Z0007_EKDN230004371-1A_HNHMFDSX5_L4.log 2>&1 
mkdir -p /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/
trimmomatic PE -threads 5 -phred33 -trimlog /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/Z0007_EKDN230004371-1A_HNHMFDSX5_L4.txt /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_1.fq.gz /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_2.fq.gz -baseout /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/trim_Z0007_EKDN230004371-1A_HNHMFDSX5_L4.fq.gz ILLUMINACLIP:/prj/ycc-backup/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:2:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/trim_Z0007_EKDN230004371-1A_HNHMFDSX5_L4.log 2>&1 
rm /prj/ycc-backup/ParusMA/pma.01.09/fastq/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_.*
mkdir -p /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc2/
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc2/
fastqc -o /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc2/ -f fastq /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/trim_Z0007_EKDN230004371-1A_HNHMFDSX5_L4_1P.fq.gz /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/trim/trim_Z0007_EKDN230004371-1A_HNHMFDSX5_L4_2P.fq.gz > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/qc2/Z0007_EKDN230004371-1A_HNHMFDSX5_L4.log 2>&1 
