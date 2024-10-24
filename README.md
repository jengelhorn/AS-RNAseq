# AS-RNAseq
Scripts for allele-specific analysis of strand-specific RNAseq data

## Description
These scripts allow counting of strand-specific RNAseq reads based on SNP positions between the two parental genomes in an F1 hybrid. 
They assume omni hybrids with B73 as a common mother, so B73 as one allele is hardcoded for the moment. 

## halliftover_and_depuplicate_of_SNPs_for_RNA_git.sh

Example command:

`./halliftover_and_depuplicate_of_SNPs_for_RNA_git.sh genotype WorkDIR SNPs.tsv hal`

This script generates two files containing the positions of all 1:1 mappable SNPs in an F1 hybrid with unique chromosome names for each parent (where one parent is B73 and one is [genotype]).

- requires the bedops toolkit (bedops.readthedocs.io) and the path to this toolkit needs to be updated in the script
- expects SNP file in tsv format from Cactus (or any other 0-based! SNP file in the format CHR POS)
- expects a hal file where one genome is called B73 and the other has the same name as specified in "genotype", e.g. Mo17
- returns two files named B73.[genotype].ID.hallifted.SNPs.bed  [genotype].B73.ID.hallifted.SNPs.bed


# SNP_Gene_files_for_RNA_seq_git.sh

Example command:

`./SNP_Gene_files_for_RNA_seq_git.sh B73.gff parent2.gff SNP_file chromfile Gene_translate_file genotype DIR`

This script generates two files containing all SNPs located in B73 genes that have a corresponding gene model in the paternal genome, in B73 and paternal coordinates, respectively. Each file also includes the corresponding coordinate in the other parent as an entry as well as the name of the gene it is overlapping with (Both B73 and paternal name) and strand information

- expects, for each parental genome, a gff file with exon annotation. Chromosomes must be named B37-[chromosomename] and [genotype]-[chromosomename]
- expects translated SNP files from halliftover_and_depuplicate_of_SNPs_for_RNA_new_lines.sh, chromosomename being the original names used in the hal and bed files in the previous steps
- expects a genome file size for a diploid genome in the form CHROM LENGTH for each Chromosome, where the two parental chromosomes are named B73-[chromosomename] and [genotype]-[chromosomename]
- expects a gene name translation file in the format [paternal gene name] [B73 gene name], if this is not available, use SNP_Gene_files_for_RNA_seq_without_translate_file_git.sh
- generates two files B73.[genotype].SNPs.genes.B73_strand.bed and [genotype].B73.SNPs.genes.B73_strand.bed in format [Chr B73] [Start B73] [Stop B73] [ID] [B73 geneID] [B73 strand]  and [Chr Pat] [Start Pat] [Stop Pat] [ID] [B73 geneID] [Pat strand] [Pat geneID]
- the ID denotes the coordinates in the respective other genome (Chr, position) plus the reference (B73) and alternative allele (Paternal)



# SNP_Gene_files_for_RNA_seq_without_translate_file_git.sh

Example command:

`./SNP_Gene_files_for_RNA_seq_without_translate_file_git.sh B73.gff B73.gff parent2.gff SNP_file chromfile genotype DIR`

This script is employed when no translation file for gene names is available. Such a translation file is generated within this script based on SNP positions overlapping with the genes and only genes on the same strand in both genotypes are retained. Otherwise, the same input as for SNP_Gene_files_for_RNA_seq_git.sh is required.


# Count_SNP_reads_per_gene_for_RNA_seq_git.sh

Example command:

`./Count_SNP_reads_per_gene_for_RNA_seq_new10_22.sh bam1 bam2 bam3 SNPs_genes_B73 SNPs_genes_NAM genotype DIR out`

This script generates read counts for each maternal/paternal gene pair based on counts on SNPs, i.e. only reads on SNPs are counted, double-counting  is prevented. Thus, the result returns the number of reads that overlap with SNP positions (the positions where allele-specific counting is possible since there are sequence differences and the position is 1:1 mappable between the two genomes) in each gene. Only counts from SNPs that carry at least one read in each parent are counted to avoid artefacts due gene length variation.

- expects bam files of RNA-seq reads mapped to a diploid genome, three replicates are assumed
- expects the output files of SNP_Gene_files_for_RNA_seq_git.sh
- generates [genotype].all_genes_with_reads.counts_[out].txt in format [B73 geneID] [B73 counts Rep.1] [B73 counts Rep.2] [B73 counts Rep.3] [Pat counts Rep.1] [Pat counts Rep.2] [Pat counts Rep.3]
- counts can be used for testing of allele-specific expression in e.g. DEseq2 using the parent as treatment

# normalise_allele_counts_git.sh

Example command:

`./normalise_allele_counts.sh infile outfile`

This script normalises allele-specific RNAseq count data and returns counts per million values for each replicate. Note that transcript abundance values counted in this way depend on the number of SNPs present within a gene, thus they only reflect the allelic distribution, not the amount of transcript relative to other genes. The same genes in different conditions, however, can be compared. 

-expects an infile of the format GENE Counts_Rep1_parent1  Counts_Rep2_parent1  Counts_Rep3_parent1  Counts_Rep1_parent2 Counts_Rep2_parent2 Counts_Rep3_parent2 (assuming 3 reps with columns 2 and 5, 3 and 6 and 4 and 7 to be normalised together)


