#NOTES:
#Replace "…>49…" with the appropriate read coverage filter threshold
#Replace "….pf" with the appropriate file extension for annotated snp-calling (atlas-snp, &c.) outputs

#Illumina // atlas-snp
#mkdir $2
#cd $1
#cutter="cut -f1,2,5,50-55,66-68"
#for a in `ls *somatic*.pf | cut -d. -f1 | sort -u`; do echo $a; out=$a.somatic.txt; awk -F"\t" '$50>49 && $51>49' $a*somatic*.pf | $cutter >> $out; done
#for a in `ls *germline*.pf | cut -d. -f1 | sort -u`; do echo $a; out=$a.germline.txt; awk -F"\t" '$50>49 && $51>49' $a*germline*.pf | $cutter >> $out; done
#mv *.txt $2

#SOLiD
mkdir $2
cd $1/somatic
cutter="cut -f1,2,3,8,11,13,14,19,28,38"
for a in `ls *.filtered | cut -d. -f1 | sort -u`; do echo $a; out=$a.somatic.txt; awk -F"\t" '$8>49 && $13>49' $a*.filtered | $cutter >> $out; done
mv *.txt $2
cd $1/germline
cutter="cut -f1,2,3,8,11,13,14,19,28,38"
for a in `ls *.filtered | cut -d. -f1 | sort -u`; do echo $a; out=$a.germline.txt; awk -F"\t" '$8>49 && $13>49' $a*.filtered | $cutter >> $out; done
mv *.txt $2