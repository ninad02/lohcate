#NOTES:
#Replace "…>49…" with the appropriate read coverage filter threshold
#Replace "….pf" with the appropriate file extension for annotated snp-calling (atlas-snp, &c.) outputs

#Illumina // atlas-snp
mkdir $2
cd $1
cutter="cut -f1,2,5,50-55,66-68"
for a in `ls *somatic*.$3 | cut -d. -f1 | sort -u`; do echo $a; out=$a.somatic.txt; awk -F"\t" '$50>49 && $51>49' $a*somatic*.$3 | $cutter >> $out; done
for a in `ls *germline*.$4 | cut -d. -f1 | sort -u`; do echo $a; out=$a.germline.txt; awk -F"\t" '$50>49 && $51>49' $a*germline*.$4 | $cutter >> $out; done
mv *.txt $2