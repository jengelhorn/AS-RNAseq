#!/bin/bash
#run by ./Count_SNP_reads_per_gene_for_RNA_seq.sh bam1 bam2 bam3 SNPs_genes_B73 SNPs_genes_NAM genotype DIR out
EXPECTED_ARGS=8
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: `basename $0` bam1 bam2 bam3 SNPs_genes_B73 SNPs_genes_NAM genotype DIR out"
  exit $E_BADARGS
fi
bam1=$1
bam2=$2
bam3=$3
SNPsB=$4
SNPsN=$5
g=$6
DIR=$7
out=$8
bai1=${bam1}.bai
bai2=${bam2}.bai
bai3=${bam3}.bai
# 
 if [ -f "$bam1" ] 
 then
echo "File $bam1  exists"
  else
echo "File $bam1  not found"
exit $E_BADARGS
fi

 if [ -f "$bai1" ] 
 then
echo "File $bam1 index  exists"
  else
echo "File $bam1 index does not exists, creating"
samtools index $bam1
fi

  if [ -f "$bam2" ] 
 then
echo "File $bam2  exists"
  else
echo "File $bam2  not found"
exit $E_BADARGS
fi

 if [ -f "$bai2" ] 
 then
echo "File $bam2 index  exists"
  else
echo "File $bam2 index does not exists, creating"
samtools index $bam2
fi

 if [ -f "$bam3" ] 
 then
echo "File $bam3  exists"
  else
echo "File $bam3  not found"
exit $E_BADARGS
fi

 if [ -f "$bai3" ] 
 then
echo "File $bam3 index  exists"
  else
echo "File $bam3 index does not exists, creating"
samtools index $bam3
fi

   if [ -f "$SNPsB" ] 
 then
echo "File $SNPsB  exists"
  else
echo "File $SNPsB  not found, run SNP_Gene_files_for_RNA_seq.sh to generate"
exit $E_BADARGS
fi
 
   if [ -f "$SNPsN" ] 
 then
echo "File $SNPsN  exists"
  else
echo "File $SNPsN  not found, run SNP_Gene_files_for_RNA_seq.sh to generate"
exit $E_BADARGS
fi


    if [ -d "$DIR" ] 
 then
echo "Folder $DIR  exists"

cd $DIR

echo "generating strand-corrected bed files"

bamToBed -i $bam1 | gawk -v OFS='\t' '{if(substr($4,length($4),1)=="2"){print $0};if(substr($4,length($4),1)=="1"){if($6=="+"){print $1,$2,$3,$4,$5,"-"};if($6=="-"){print $1,$2,$3,$4,$5,"+"}}}' > Rep1_${out}.bed
bamToBed -i $bam2 | gawk -v OFS='\t' '{if(substr($4,length($4),1)=="2"){print $0};if(substr($4,length($4),1)=="1"){if($6=="+"){print $1,$2,$3,$4,$5,"-"};if($6=="-"){print $1,$2,$3,$4,$5,"+"}}}' > Rep2_${out}.bed
bamToBed -i $bam3 | gawk -v OFS='\t' '{if(substr($4,length($4),1)=="2"){print $0};if(substr($4,length($4),1)=="1"){if($6=="+"){print $1,$2,$3,$4,$5,"-"};if($6=="-"){print $1,$2,$3,$4,$5,"+"}}}' > Rep3_${out}.bed


echo "intersecting read beds with SNPs"

intersectBed -a Rep1_${out}.bed -b $SNPsB -wa -wb -s > B73.${g}.reads.Rep1_${out}.bed

intersectBed -a Rep1_${out}.bed -b $SNPsN -wa -wb -s > ${g}.B73.reads.Rep1_${out}.bed

intersectBed -a Rep2_${out}.bed -b $SNPsB -wa -wb -s > B73.${g}.reads.Rep2_${out}.bed

intersectBed -a Rep2_${out}.bed -b $SNPsN -wa -wb -s > ${g}.B73.reads.Rep2_${out}.bed

intersectBed -a Rep3_${out}.bed -b $SNPsB -wa -wb -s > B73.${g}.reads.Rep3_${out}.bed

intersectBed -a Rep3_${out}.bed -b $SNPsN -wa -wb -s > ${g}.B73.reads.Rep3_${out}.bed

echo "Get SNPs that carry reads in both alleles"

cat B73.${g}.reads.Rep1_${out}.bed B73.${g}.reads.Rep2_${out}.bed B73.${g}.reads.Rep3_${out}.bed | gawk -v OFS='\t' '{print $7":"$9":"$11}'|sort|uniq > B73.${g}.SNPs_with_reads_${out}.txt

cat ${g}.B73.reads.Rep1_${out}.bed ${g}.B73.reads.Rep2_${out}.bed ${g}.B73.reads.Rep3_${out}.bed | gawk -v OFS='\t' '{split($10,var,"\.");print var[1]":"var[2]":"$11}'|sort|uniq > ${g}.B73.SNPs_with_reads_${out}.txt

gawk -v OFS='\t' 'NR==FNR {a[$1]=$1; next} ($1 in a){print $0}' ${g}.B73.SNPs_with_reads_${out}.txt B73.${g}.SNPs_with_reads_${out}.txt  > B73.${g}.SNPs_with_reads_both_${out}.txt

echo "cont reads per gene on SNPs with reads in both alleles"

for r in Rep1 Rep2 Rep3; do gawk -v OFS='\t' 'NR==FNR {a[$1]=$1; next} ($7":"$9":"$11 in a){print $4,$11}' B73.${g}.SNPs_with_reads_both_${out}.txt B73.${g}.reads.${r}_${out}.bed |sort|uniq| gawk -v OFS='\t' '{print $2}'|sort| uniq -c  > B73.${g}.readcounts.genes.${r}_${out}.bed; done

for r in Rep1 Rep2 Rep3; do gawk -v OFS='\t' 'NR==FNR {a[$1]=$1; next} {split($10,var,"\.")} (var[1]":"var[2]":"$11 in a){print $4,$11}' B73.${g}.SNPs_with_reads_both_${out}.txt ${g}.B73.reads.${r}_${out}.bed |sort|uniq| gawk -v OFS='\t' '{print $2}'|sort| uniq -c  > ${g}.B73.readcounts.genes.${r}_${out}.bed ; done

echo "paste lists together, add 0 for no counts"

cat B73.${g}.readcounts.genes.Rep1_${out}.bed B73.${g}.readcounts.genes.Rep2_${out}.bed B73.${g}.readcounts.genes.Rep3_${out}.bed ${g}.B73.readcounts.genes.Rep1_${out}.bed ${g}.B73.readcounts.genes.Rep2_${out}.bed ${g}.B73.readcounts.genes.Rep3_${out}.bed | gawk -v OFS='\t' '{print $2}'  |sort|uniq > ${g}.all_genes_with_reads_${out}.txt

for r in Rep1 Rep2 Rep3; do gawk -v OFS='\t' 'NR==FNR {a[$2]=$2;text[$2]=$1; next} ($1 in a){print $1,text[$1]};!($1 in a){print $1,"0"} ' B73.${g}.readcounts.genes.${r}_${out}.bed ${g}.all_genes_with_reads_${out}.txt > B73.${g}.all_genes_with_reads.${r}_counts_${out}.txt; done

for r in Rep1 Rep2 Rep3; do gawk -v OFS='\t' 'NR==FNR {a[$2]=$2;text[$2]=$1; next} ($1 in a){print $1,text[$1]};!($1 in a){print $1,"0"} ' ${g}.B73.readcounts.genes.${r}_${out}.bed ${g}.all_genes_with_reads_${out}.txt > ${g}.B73.all_genes_with_reads.${r}_counts_${out}.txt; done


paste B73.${g}.all_genes_with_reads.Rep1_counts_${out}.txt B73.${g}.all_genes_with_reads.Rep2_counts_${out}.txt B73.${g}.all_genes_with_reads.Rep3_counts_${out}.txt ${g}.B73.all_genes_with_reads.Rep1_counts_${out}.txt ${g}.B73.all_genes_with_reads.Rep2_counts_${out}.txt ${g}.B73.all_genes_with_reads.Rep3_counts_${out}.txt | gawk -v OFS='\t' '{print $1,$2,$4,$6,$8,$10,$12}'  > ${g}.all_genes_with_reads.counts_${out}.txt

  else
echo "Folder $DIR  not found"
exit $E_BADARGS
fi



