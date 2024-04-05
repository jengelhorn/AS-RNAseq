# AS-RNAseq
Scripts for allele-specific analysis of strand-specific RNAseq data

## Description
These scripts allow counting of strand-specific RNAseq reads based on SNP positions between the two parental genomes in an F1 hybrid. 
They assume omni hybrids with B73 as a common mother, so for the moment B73 as one allele is hard coded but this will be flexibilised in the future. 

## halliftover_and_depuplicate_of_SNPs_for_RNA_git.sh

Example command:

`./halliftover_and_depuplicate_of_SNPs_for_RNA_git.sh genotype WorkDIR SNPs.tsv hal`

This script generates two files containing the positions of all 1:1 mappable SNPs in an F1 hybrid with unique chromosome names for each parent (where one parent is B73 and one is [genotype]).

- expects SNP file in tsv format from Cactus (or any other 0-based! SNP file in the format CHR POS)
- expects a hal file where one genome is called B73 and the other has the same name as specified in "genotype", e.g. Mo17
- returns two files named B73.[genotype].ID.hallifted.SNPs.bed  [genotype].B73.ID.hallifted.SNPs.bed

To do: bedops path needs to become a variable

# SNP_Gene_files_for_RNA_seq_git.sh

Example command:

`./SNP_Gene_files_for_RNA_seq_git.sh B73.gff parent2.gff SNP_file chromfile Gene_translate_file genotype DIR`

This script generates XXX

- expects, for each parental genome, a gff file with exon annotation. Chromosomes must be named B37-[chromosomename] and [genotype]-[chromosomename]
- expects translated SNP files from halliftover_and_depuplicate_of_SNPs_for_RNA_new_lines.sh, chromosomename being the original names used in the hal and bed files in the previous steps
- expects a genome file size for a diploid genome in the form CHROM LENGTH for each Chromosome, where the two parental chromosomes are named B73-[chromosomename] and [genotype]-[chromosomename]
- expects a gene name translation file in the format XXX, if this is not available, use 
- generates two files B73.[genotype].SNPs.genes.B73_strand.bed and [genotype].B73.SNPs.genes.B73_strand.bed in format XXX

  To do: description, better annotation, flexibilisation of gene name length if possible

# SNP_Gene_files_for_RNA_seq_without_translate_file_git.sh

Example command:

`./SNP_Gene_files_for_RNA_seq_without_translate_file_git.sh B73.gff B73.gff parent2.gff SNP_file chromfile genotype DIR`

This script is employed when no translation file for gene names is available. Such a translation file is generated within this script based on SNP positions overlapping with the genes and only genes with on the same strand in both genotypes are retained. Otherwise, the same input as for SNP_Gene_files_for_RNA_seq_git.sh is required.


# Count_SNP_reads_per_gene_for_RNA_seq_new10_22.sh

Example command:

`./Count_SNP_reads_per_gene_for_RNA_seq_new10_22.sh bam1 bam2 bam3 SNPs_genes_B73 SNPs_genes_NAM genotype DIR out`

This script generates XXX

- expects bam files of RNA-seq reads mapped to a diploid genome, three replicates are assumed
- expects the output files of SNP_Gene_files_for_RNA_seq_git.sh
- generates [genotype].all_genes_with_reads.counts_[out].txt in format XXX
- counts can be used for testing of allele-specific expression in e.g. DEseq2 using the parent as treatment

To do: explain how genes with 0 reads are excluded

# normalise_allele_counts_git.sh

Example command:

`./normalise_allele_counts.sh infile outfile`

This script normalises allele-specific RNAseq count data and returns counts per million counts for each replicate

-expects an infile of the format GENE Counts_Rep1_parent1  Counts_Rep2_parent1  Counts_Rep3_parent1  Counts_Rep1_parent2 Counts_Rep2_parent2 Counts_Rep3_parent2 (assuming 3 reps with columns 2 and 5, 3 and 6 and 4 and 7 to be normalised together)


