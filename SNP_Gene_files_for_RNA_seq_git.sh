#!/bin/bash
#run by ./SNP_Gene_files_for_RNA_seq_git.sh B73.gff NAM.gff SNP_file chromfile Gene_translate_file genotype DIR
EXPECTED_ARGS=7
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` B73.gff parent2.gff SNP_file chromfile Gene_translate_file genotype DIR"
  exit $E_BADARGS
fi
Bgff=$1
Ngff=$2
SNPs=$3
dat=$4
Trans=$5
g=$6
DIR=$7
# 
 if [ -f "$Bgff" ] 
 then
echo "File $Bgff  exists"
  else
echo "File $Bgff  not found"
fi

 if [ -f "$Ngff" ] 
 then
echo "File $Ngff  exists"
  else
echo "File $Ngff  not found"
fi

  if [ -f "$SNPs" ] 
 then
echo "File $SNPs  exists"
  else
echo "File $SNPs  not found"
fi
 
   if [ -f "$dat" ] 
 then
echo "File $dat  exists"
  else
echo "File $dat  not found"
fi
 
    if [ -f "$Trans" ] 
 then
echo "File $Trans  exists"
  else
echo "File $Trans  not found"
fi

    if [ -d "$DIR" ] 
 then
echo "Folder $DIR  exists"

cd $DIR

echo "generating exon bedfiles B73"

gawk -v OFS='\t' '{if($3=="exon" && $7=="+"){print $1,$4-1,$5,substr($9,8,15)}}' $Bgff  | sortBed > B73.exons.sorted.forward.bed
gawk -v OFS='\t' '{if($3=="exon" && $7=="-"){print $1,$4-1,$5,substr($9,8,15)}}' $Bgff  | sortBed > B73.exons.sorted.reverse.bed

echo "adding B73 gene information to SNPs"

intersectBed -a $SNPs -b B73.exons.sorted.forward.bed -wa -wb | gawk -v OFS='\t' '{print $1,$2,$3,$4,$8}'|sort | uniq |  gawk -v OFS='\t' 'BEGIN{prev=0;warn=0;pc=0}{if(NR==1){prev=$4;pc=$0}else{if($4==prev){warn=1}else{if(warn==0){print pc};warn=0;prev=$4;pc=$0}}}END{if(warn==0){print pc}}' | gawk -v OFS='\t' -v g=${g} '{print $1,$2,$3,g"-"$4,$5,"+"}' | sortBed -g $dat > B73.${g}.SNPs.genes.forward.bed

intersectBed -a $SNPs -b B73.exons.sorted.reverse.bed -wa -wb | gawk -v OFS='\t' '{print $1,$2,$3,$4,$8}'|sort | uniq |  gawk -v OFS='\t' 'BEGIN{prev=0;warn=0;pc=0}{if(NR==1){prev=$4;pc=$0}else{if($4==prev){warn=1}else{if(warn==0){print pc};warn=0;prev=$4;pc=$0}}}END{if(warn==0){print pc}}' | gawk -v OFS='\t' -v g=${g} '{print $1,$2,$3,g"-"$4,$5,"-"}' | sortBed -g $dat > B73.${g}.SNPs.genes.reverse.bed

cat B73.${g}.SNPs.genes.forward.bed B73.${g}.SNPs.genes.reverse.bed > B73.${g}.SNPs.genes.B73_strand.1.bed

echo "writing out " ${g}  " gene strand and official pairs"

gawk -v OFS='\t' '{if($3=="gene"){print substr($9,4,15),$7}}' $Ngff > ${g}.gene.strand.txt

gawk -v OFS='\t' 'NR==FNR {a[$1]=$1;text[$1]=$0; next} ($1 in a){print text[$1],$2}' $Trans ${g}.gene.strand.txt > ${g}.B73.genes.official_pairs.strand.csv

echo "add " ${g}  " information to B73 SNPs"

gawk -v OFS='\t' '{split($4,var,"\."); print var[1],var[2]-1,var[2],$1"."$3"."var[3]"."var[4],$5}' B73.${g}.SNPs.genes.B73_strand.1.bed > ${g}.B73.SNPs.genes.B73.bed

gawk -v OFS='\t' 'NR==FNR {a[$2]=$2;text[$2]=$1;text2[$2]=$3; next} ($5 in a){print $0,text2[$5],text[$5]}' ${g}.B73.genes.official_pairs.strand.csv ${g}.B73.SNPs.genes.B73.bed > ${g}.B73.SNPs.genes.${g}_strand.bed

gawk -v OFS='\t' 'NR==FNR {a[$2]=$2;text[$2]=$1;text2[$2]=$3; next} ($5 in a){print $0}' ${g}.B73.genes.official_pairs.strand.csv B73.${g}.SNPs.genes.B73_strand.1.bed > B73.${g}.SNPs.genes.B73_strand.bed


  else
echo "Folder $DIR  not found"
fi



