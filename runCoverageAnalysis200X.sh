
#!/bin/bash 

awk 'BEGIN {print "Exon" "\t" "D_0" "\t" "D_5" "\t" "D_10" "\t" "D_20" "\t" "D_50" "\t" "D_100" "\t" "D_150" "\t" "D_200"}' > truseqDepth200X.txt;
awk 'BEGIN {print "Exon" "\t" "D_0" "\t" "D_5" "\t" "D_10" "\t" "D_20" "\t" "D_50" "\t" "D_100" "\t" "D_150" "\t" "D_200"}' > truseqBreadth200X.txt;
awk 'BEGIN {print "chromosome" "\t" "No_Coverage" "\t" "at1X" "\t" "at5X" "\t" "at10X" "\t" "at20X" "\t" "at50X" "\t" "at100X" "\t" "at150X" "\t" "at200X"}' > truseqCoverage200X.txt;


#===================================TruSeq coverage
for i in `find chr* -name "truseqCoverage.bed" -printf "%h\n"`; do

	echo Calculating Truseq Depth and Breadth of $i
	CURRENT=`pwd`
	BASENAME=`basename $CURRENT`
	python  /home/yqin/Pipeline/Orion_QC/PIPELINE/DNA-seq-0.1.3/runOnTruSeqCoverage.py $i $CURRENT
	
	total=`egrep "^$i\W" $i/truseqCoverage.bed | cut -f9 | uniq | awk '{sum1 += $1} END {print sum1}'`
   t0=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7==0 {sum1 += $8} END {print sum1/total}'`
   t1=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=1 {sum1 += $8} END {print sum1/total}'`
   t5=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=5 {sum1 += $8} END {print sum1/total}'`
   t10=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=10 {sum1 += $8} END {print sum1/total}'`
   t20=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=20 {sum1 += $8} END {print sum1/total}'`
   t50=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=50 {sum1 += $8} END {print sum1/total}'`
   t100=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=100 {sum1 += $8} END {print sum1/total}'`
   t150=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=150 {sum1 += $8} END {print sum1/total}'`
   t200=`egrep "^$i\W" $i/truseqCoverage.bed | awk -v total=$total '$7>=200 {sum1 += $8} END {print sum1/total}'`

   awk -v chr=$i -v t0=$t0 -v t1=$t1 -v t5=$t5 -v t10=$t10 -v t20=$t20 -v t50=$t50 -v t100=$t100 -v t150=$t150 -v t200=$t200 'BEGIN {print chr "\t" t0 "\t" t1 "\t" t5 "\t" t10 "\t" t20 "\t" t50 "\t" t100 "\t" t150 "\t" t200}' >> truseqCoverage200X.txt

done


#===================================Exon coverage
awk 'BEGIN {print "Exon" "\t" "D_0" "\t" "D_5" "\t" "D_10" "\t" "D_20" "\t" "D_50" "\t" "D_100" "\t" "D_150" "\t" "D_200"}' > exondepth200X.txt;
awk 'BEGIN {print "Exon" "\t" "D_0" "\t" "D_5" "\t" "D_10" "\t" "D_20" "\t" "D_50" "\t" "D_100" "\t" "D_150" "\t" "D_200"}' > exonbreadth200X.txt;
awk 'BEGIN {print "chromosome" "\t" "No_Coverage" "\t" "at1X" "\t" "at5X" "\t" "at10X" "\t" "at20X" "\t" "at50X" "\t" "at100X" "\t" "at150X" "\t" "at200X"}' > exonCoverage200X.txt;

for i in `find chr* -name "exoncoverage.bed" -printf "%h\n"`; do

	echo Calculating exoncoverage Depth and Breadth of $i
	CURRENT=`pwd`
	BASENAME=`basename $CURRENT`
	python  /home/yqin/Pipeline/Orion_QC/PIPELINE/DNA-seq-0.1.3/runOnExonCoverage.py $i $CURRENT

	total=`egrep "^$i\W" $i/exoncoverage.bed | cut -f9 | uniq | awk '{sum1 += $1} END {print sum1}'`
  t0=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7==0 {sum1 += $8} END {print sum1/total}'`
  t1=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=1 {sum1 += $8} END {print sum1/total}'`
  t5=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=5 {sum1 += $8} END {print sum1/total}'`
  t10=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=10 {sum1 += $8} END {print sum1/total}'`
  t20=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=20 {sum1 += $8} END {print sum1/total}'`
  t50=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=50 {sum1 += $8} END {print sum1/total}'`
  t100=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=100 {sum1 += $8} END {print sum1/total}'`
  t150=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=150 {sum1 += $8} END {print sum1/total}'`
  t200=`egrep "^$i\W" $i/exoncoverage.bed | awk -v total=$total '$7>=200 {sum1 += $8} END {print sum1/total}'`


  awk -v chr=$i -v t0=$t0 -v t1=$t1 -v t5=$t5 -v t10=$t10 -v t20=$t20 -v t50=$t50 -v t100=$t100 -v t150=$t150 -v t200=$t200 'BEGIN {print chr "\t" t0 "\t" t1 "\t" t5 "\t" t10 "\t" t20 "\t" t50 "\t" t100 "\t" t150 "\t" t200}' >> exonCoverage200X.txt

done

#===================================RefSeq coverage
awk 'BEGIN {print "Exon" "\t" "D_0" "\t" "D_5" "\t" "D_10" "\t" "D_20" "\t" "D_50" "\t" "D_100" "\t" "D_150" "\t" "D_200"}' > refseqDepth200X.txt;
awk 'BEGIN {print "Exon" "\t" "D_0" "\t" "D_5" "\t" "D_10" "\t" "D_20" "\t" "D_50" "\t" "D_100" "\t" "D_150" "\t" "D_200"}' > refseqBreadth200X.txt;
awk 'BEGIN {print "chromosome" "\t" "No_Coverage" "\t" "at1X" "\t" "at5X" "\t" "at10X" "\t" "at20X" "\t" "at50X" "\t" "at100X" "\t" "at150X" "\t" "at200X"}' > refseqCoverage200X.txt;

for i in `find chr* -name "refseqCoverage.bed" -printf "%h\n"`; do

	echo Calculating refseqCoverage Depth and Breadth of $i
	CURRENT=`pwd`
	BASENAME=`basename $CURRENT`
	python  /home/yqin/Pipeline/Orion_QC/PIPELINE/DNA-seq-0.1.3/runOnRefSeqCoverage.py $i $CURRENT

   total=`egrep "^$i\W" $i/refseqCoverage.bed | cut -f7 | uniq | awk '{sum1 += $1} END {print sum1}'`
   t0=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5==0 {sum1 += $6} END {print sum1/total}'`
   t1=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=1 {sum1 += $6} END {print sum1/total}'`
   t5=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=5 {sum1 += $6} END {print sum1/total}'`
   t10=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=10 {sum1 += $6} END {print sum1/total}'`
   t20=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=20 {sum1 += $6} END {print sum1/total}'`
   t50=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=50 {sum1 += $6} END {print sum1/total}'`
   t100=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=100 {sum1 += $6} END {print sum1/total}'`
   t150=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=150 {sum1 += $6} END {print sum1/total}'`
   t200=`egrep "^$i\W" $i/refseqCoverage.bed | awk -v total=$total '$5>=200 {sum1 += $6} END {print sum1/total}'`

   awk -v chr=$i -v t0=$t0 -v t1=$t1 -v t5=$t5 -v t10=$t10 -v t20=$t20 -v t50=$t50 -v t100=$t100 -v t150=$t150 -v t200=$t200 'BEGIN {print chr "\t" t0 "\t" t1 "\t" t5 "\t" t10 "\t" t20 "\t" t50 "\t" t100 "\t" t150 "\t" t200}' >> refseqCoverage200X.txt
done
#===================================

#SNP depth
#for i in `find chr* -name "snpCoverage.bed" -printf "%h\n"`; do
#	echo Calculating snpCoverage of $i
        ##awk 'NR>=1' $i/snpCoverage.bed >> snpdepth.txt;
#       awk '{if ($1!="all") print $0}' $i/snpCoverage.bed >> snpdepth.txt;

#done
