#!/bin/sh
#SBATCH --job-name=star_indices # Job name
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=d.ence@mail.ufl.edu  # Where to send mail
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=60gb                  # Total memory limit
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --output=star_indices_%j.out     # Standard output and error log
#SBATCH --qos=kirst-b
#SBATCH --account=kirst
#SBATCH --partition=hpg1-compute

#export OMP_NUM_THREADS=4

# Load required modules; for example, if your program was
# compiled with Intel compiler, use the following

module load star/2.7.5c

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ./v2.01_genome_annots_index/ --genomeFastaFiles /blue/kirst/share/Genomes/Ptaeda_v201/genome_annotation/annotation_2.01_Feb_2020/Pita.2_01.cds.fa --limitGenomeGenerateRAM=36176571317 --genomeSAindexNbases 11

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ./clustered_longer_than_500.est_0.99/ --genomeFastaFiles /home/d.ence/d.ence_peter_share/pinus_taeda_L/miami_ohio_ESTs_ptaeda_rbh/combining_PacBio_miami_of_ohio/adding_V2.01_t3352/t3352.cdna.v1.contig.Pita.2_01.Feb_2020.cds.longer_than_500.est_0.99  --limitGenomeGenerateRAM=58220435040 --genomeSAindexNbases 11

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir t3352.cdna.v1.contig --genomeFastaFiles /home/d.ence/d.ence_peter_share/pinus_taeda_L/miami_ohio_ESTs_ptaeda_rbh/t3352.cdna.v1.contig.fasta  --limitGenomeGenerateRAM=58220435040 --genomeSAindexNbases 11

STAR --runThreadN 1 --runMode genomeGenerate --genomeDir ./PacBio --genomeFastaFiles /home/d.ence/d.ence_peter_share/pinus_taeda_L//miami_ohio_ESTs_ptaeda_rbh/combining_PacBio_miami_of_ohio/PacBio_transcriptome/Pita.IU.TranscriptomeMainsV1.102013.fasta  --limitGenomeGenerateRAM=58220435040 --genomeSAindexNbases 11
# Run your program with correct path and command line options
#bowtie2-build /blue/kirst/share/Genomes/Ptaeda_v201/genome_annotation/annotation_2.01_Feb_2020/Pita.2_01.cds.fa Pita.2_01.cds
#bowtie2-build ../../miami_ohio_ESTs_ptaeda_rbh/combining_PacBio_miami_of_ohio/adding_V2.01_t3352/t3352.cdna.v1.contig.Pita.2_01.Feb_2020.cds.longer_than_500.est_0.99 clustered_longer_than_500.est_0.99
#bowtie2-build ../../miami_ohio_ESTs_ptaeda_rbh/combining_PacBio_miami_of_ohio/PacBio_transcriptome/Pita.IU.TranscriptomeMainsV1.102013.fasta PacBio
#bowtie2-build ../../miami_ohio_ESTs_ptaeda_rbh/t3352.cdna.v1.contig.fasta t3352.cdna.v1.contig

#rsem-calculate-expression --bowtie2 --num-threads 10 --paired-end /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L007_R1.fastq.gz /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L008_R2.fastq.gz Pita.2_01.cds Pita.2_01.cds_test
#
#rsem-calculate-expression --bowtie2 --num-threads 10 --paired-end /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L007_R1.fastq.gz /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L008_R2.fastq.gz clustered_longer_than_500.est_0.99 clustered_longer_than_500.est_0.99_test

#rsem-calculate-expression --bowtie2 --num-threads 10 --paired-end /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L007_R1.fastq.gz /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L008_R2.fastq.gz PacBio PacBio_test

#rsem-calculate-expression --bowtie2 --num-threads 10 --paired-end /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L007_R1.fastq.gz /home/d.ence/d.ence_peter_share/UF_PETER_Terpene/concatted__fastqs/Sample_UF_Derv_001_ATCACG_L008_R2.fastq.gz t3352.cdna.v1.contig t3352.cdna.v1.contig_test
