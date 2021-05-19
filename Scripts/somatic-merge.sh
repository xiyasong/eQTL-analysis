#use this code to get the PASS mutation in every patients

for i in `find . -maxdepth 1 -type d -printf '%f\n'`;
do 
cd $i;
bcftools view -f PASS $i.somatic.filtered.vcf.gz > $i.somatic.filtered.PASS.vcf.gz;
cd -;
done

#actuually the gained result vcf isn't zip well
#change the name then gzip

for i in `find . -maxdepth 1 -type d -printf '%f\n'`;
do 
cd $i;
mv $i.somatic.filtered.PASS.vcf.gz $i.somatic.filtered.PASS.vcf;
gzip $i.somatic.filtered.PASS.vcf;
cd -;
done



#then build the index
#bgzip rather than gzip
for i in `find . -maxdepth 1 -type d -printf '%f\n'`;
do 
cd $i;
gunzip $i.somatic.filtered.PASS.vcf.gz;
bgzip $i.somatic.filtered.PASS.vcf;
tabix -p vcf $i.somatic.filtered.PASS.vcf.gz
cd -;
done

#process every somatic mutations file into only have tumor sample
for i in `find . -maxdepth 1 -type d -printf '%f\n'`;
do 
cd $i;
zcat $i.somatic.filtered.PASS.vcf.gz | cut -f1-10 > $i.somatic.filtered.PASS.only.tumor.vcf;
bgzip $i.somatic.filtered.PASS.only.tumor.vcf;
tabix -p vcf $i.somatic.filtered.PASS.only.tumor.vcf.gz
cd -;
done


#now move all the tumor sample file to a new folder

for i in `find . -maxdepth 1 -type d -printf '%f\n'`;
do 
cd $i
mv $i.somatic.filtered.PASS.only.tumor.vcf.gz ../vcf-merge-folder
mv $i.somatic.filtered.PASS.only.tumor.vcf.gz.tbi ../vcf-merge-folder
cd -;
done

#vcf-merge in vcf-merge-folder
vcf-merge $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') > somatic.merge.vcf


for i in `ls -1 *.tumor.vcf.gz`;
do
zcat $i | grep -Ev "^#" | wc -l >> mut_count.txt;
done

