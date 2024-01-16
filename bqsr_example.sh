#!/bin/bash
export PATH=/prj/ycc-backup/miniconda3/bin:$PATH
source activate gatk
cd /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4
/prj/ycc-backup/software/gatk-4.2.5.0/gatk HaplotypeCaller -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_sorted_dedup.bam -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_variants.vcf --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ > /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_variants.vcf -select-type-to-include SNP -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_snps.vcf --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_variants.vcf -select-type-to-include INDEL -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_indels.vcf --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_snps.vcf -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_filtered_snps.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 3.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk VariantFiltration -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -V /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_filtered_indels.vcf -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 200.0" -filter-name "SOR_filter" -filter "SOR > 10.0" --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants --exclude-filtered -V /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_filtered_snps.vcf -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr_snps.vcf --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk SelectVariants --exclude-filtered -V /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_filtered_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr_indels.vcf --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk BaseRecalibrator -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_sorted_dedup.bam --known-sites /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr_snps.vcf --known-sites /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr_indels.vcf -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_recal_data.table --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
/prj/ycc-backup/software/gatk-4.2.5.0/gatk ApplyBQSR -R /prj/ycc-backup/ParusMA/pma.01.03/ref/GCF_001522545.3_Parus_major1.1_genomic.fna -I /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_sorted_dedup.bam -bqsr /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_recal_data.table -O /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_recal_reads.bam --tmp-dir /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_tmp/ >> /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_bqsr.log 2>&1
rm /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_variants.vcf
rm /prj/ycc-backup/ParusMA/pma.01.09/ready_reads/Z0007_EKDN230004371-1A_HNHMFDSX5_L4/Z0007_EKDN230004371-1A_HNHMFDSX5_L4_raw_indels.vcf
