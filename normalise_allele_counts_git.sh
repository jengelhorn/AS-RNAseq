#!/bin/bash
#run by ./normalise_allele_counts.sh infile outfile(Assuming 3 reps with columns 2 and 5, 3 and 6 and 4 and 7 to be normalised together)
EXPECTED_ARGS=2
E_BADARGS=2
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./normalise_allele_counts.sh infile outfile"
  exit $E_BADARGS
fi


in=$1
out=$2

cd ${dir}


echo "calculate norm factor" 

N1=`awk -v OFS="\t" 'BEGIN{n1=0} {n1=n1+$2+$5} END { print n1 }' ${in} `;
N2=`awk -v OFS="\t" 'BEGIN{n2=0} {n2=n2+$3+$6} END { print n2 }' ${in} `;
N3=`awk -v OFS="\t" 'BEGIN{n3=0} {n3=n3+$4+$7} END { print n3 }' ${in} `;


echo "normalise"

awk -v OFS="\t" -v N1=${N1}  -v N2=${N2} -v N3=${N3} '{print $1,($2/N1)*1000000,($3/N2)*1000000,($4/N3)*1000000,($5/N1)*1000000,($6/N2)*1000000,($7/N3)*1000000}' ${in} > ${out}
