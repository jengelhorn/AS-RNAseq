# AS-RNAseq
Scripts for allele-specific analysis of strand-specific RNAseq data

## Description
These scripts allow counting of RNAseq reads based on SNP positions between the two parental genomes in an F1 hybrid. 

./halliftover_and_depuplicate_of_SNPs_for_RNA_new_lines.sh

- Needs SNP file in tsv format from Cactus (in our case called bed but it is not a bed and is 0 based)

/biodata/dep_psl/grp_frommer/Julia/scripts/RNA_seq/SNP_Gene_files_for_RNA_seq.sh

- Needs B73.gff NAM.gff SNP_file(from previous script) chromfile(chromosome length) Gene_translate_file(containing  genotype DIR

/biodata/dep_psl/grp_frommer/Julia/scripts/RNA_seq/Count_SNP_reads_per_gene_for_RNA_seq_new10_22.sh

/Count_SNP_reads_per_gene_for_RNA_seq_new10_22.sh bam1 bam2 bam3 SNPs_genes_B73 SNPs_genes_NAM genotype DIR out

- needs bams for three replicates and SNP files containing gene Names fro the previous script: ${g}.B73.SNPs.genes.${g}_strand.bed and B73.${g}.SNPs.genes.B73_strand.bed
