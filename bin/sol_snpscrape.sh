#NOTES:
#Replace "…>49…" with the appropriate read coverage filter threshold
#Replace "….pf" with the appropriate file extension for annotated snp-calling (atlas-snp, &c.) outputs

#SOLiD
mkdir $2
cd $1/somatic
cutter="cut -f1,2,3,8,11,13,14,19,28,38"
for a in `ls *.$3 | cut -d. -f1 | sort -u`; do echo $a; out=$a.somatic.txt; awk -F"\t" '$8>49 && $13>49' $a*.$3 | $cutter >> $out; done
mv *.txt $2
cd $1/germline
cutter="cut -f1,2,3,8,11,13,14,19,28,38"
for a in `ls *.$4 | cut -d. -f1 | sort -u`; do echo $a; out=$a.germline.txt; awk -F"\t" '$8>49 && $13>49' $a*.$4 | $cutter >> $out; done
mv *.txt $2