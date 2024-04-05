#!/bin/bash
#run by ./halliftover_and_depuplicate_of_SNPs_for_RNA_git.sh genotype 
EXPECTED_ARGS=1
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./halliftover_and_depuplicate_of_SNPs_for_RNA_git.sh genotype WorkDIR SNPs.tsv hal"
  exit $E_BADARGS
fi


g=$1
dir=$2
tsv=$3
hal=$4

export PATH="/netscratch/dep_psl/grp_frommer/Thomas/bin/bedops/bin:$PATH"



 cd ${dir}
#tsv file contains 0-based positions make bed file

 gawk -v OFS='\t' '{if(NR>1){print $1, $2, $2+1, $1"."$2+1"."$4"."$5}}' ${tsv} > ${g}.SNPs.for_liftover.bed

echo "halLiftover"
 
halLiftover  --hdf5InMemory ${hal} ${g}.SNPs.for_liftover.bed ${g} ${g}.B73.hallifted.SNPs.bed


echo "get duplicates"


gawk -v OFS='\t' '{print $4}' ${g}.B73.hallifted.SNPs.bed | sort | uniq -d > ${g}.B73.hallifted.SNPs.dups_B73coord.txt
gawk -v OFS='\t' '{print $1$2}' ${g}.B73.hallifted.SNPs.bed | sort | uniq -d > ${g}.B73.hallifted.SNPs.dups_${g}coord.txt

echo "remove duplicates"

echo "check if duplicates in" ${g}  "exist (B73 usually exist)"

if [ -s ${g}.B73.hallifted.SNPs.dups_${g}coord.txt ] 
then

echo " duplicates in " ${g}  "  exist"


gawk -v OFS='\t' 'NR==FNR{a[$1]=$1; next} !($1$2 in a){print $0}' ${g}.B73.hallifted.SNPs.dups_${g}coord.txt  ${g}.B73.hallifted.SNPs.bed > ${g}.B73.hallifted.SNPs.clean1.bed
gawk -v OFS='\t' 'NR==FNR{a[$1]=$1; next} !($4 in a){print $0}'  ${g}.B73.hallifted.SNPs.dups_B73coord.txt ${g}.B73.hallifted.SNPs.clean1.bed > ${g}.B73.hallifted.SNPs.clean2.bed


echo "renaming chromosomes"

gawk -v OFS='\t' -v ALTn=$g '{print ALTn"-"$1,$2,$3,$4}' ${g}.B73.hallifted.SNPs.clean2.bed | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > ${g}.B73.ID.hallifted.SNPs.bed

gawk -v OFS='\t'  '{{split($4,var,"."); print "B73-"var[1],var[2]-1,var[2],$1"."$3"."var[4]"."var[3]}}' ${g}.B73.hallifted.SNPs.clean2.bed  | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt  > B73.${g}.ID.hallifted.SNPs.bed

else

echo " duplicates in " ${g}  " do not exist"

gawk -v OFS='\t' 'NR==FNR{a[$1]=$1; next} !($4 in a){print $0}'  ${g}.B73.hallifted.SNPs.dups_B73coord.txt ${g}.B73.hallifted.SNPs.bed > ${g}.B73.hallifted.SNPs.clean2.bed


echo "renaming chromosomes"

gawk -v OFS='\t' -v ALTn=$g '{print ALTn"-"$1,$2,$3,$4}' ${g}.B73.hallifted.SNPs.clean2.bed | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > ${g}.B73.ID.hallifted.SNPs.bed

gawk -v OFS='\t'  '{{split($4,var,"."); print "B73-"var[1],var[2]-1,var[2],$1"."$3"."var[4]"."var[3]}}' ${g}.B73.hallifted.SNPs.clean2.bed  | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt  > B73.${g}.ID.hallifted.SNPs.bed

fi

