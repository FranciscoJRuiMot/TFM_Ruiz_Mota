#! /bin/bash

# take sequence and phenotype as input:
seq=$1
phen=$2

# convert sequence into the ancestor organism:
echo $seq | sed 's/./&\n/g' | head -n -1 > tmp.txt
awk 'BEGIN{FS=" "} NR==FNR{a[$2]=$1} NR>FNR{$1=a[$1];print}' OFS='\n' instructions.txt tmp.txt > $phen".org"

# remove temporary file:
rm tmp.txt
