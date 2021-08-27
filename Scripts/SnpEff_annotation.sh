# Tutorial: http://pcingola.github.io/SnpEff/examples/
# Step 1: Download and install SnpEff and reference.

# Move to home directory
cd
# Download and install SnpEff
curl -v -L 'https://datasetsnpeff.blob.core.windows.net/dataset/versions/snpEff_latest_core.zip?sv=2019-10-10&st=2020-09-01T00%3A00%3A00Z&se=2050-09-01T00%3A00%3A00Z&si=prod&sr=c&sig=isafOa9tGnYBAvsXFUMDGMTbsG2z%2FShaihzp7JE5dHw%3D' > snpEff_latest_core.zip
unzip snpEff_latest_core.zip

#Tutorial have said: "Once SnpEff is installed, we will enter the following commands to download the pre-built human database (GRCh37.75) that will be used to annotate our data." (Here, the GRCh38.99 is the most recent version for human in SnpEff,so instead we download this version for next step)
cd snpEff
java -jar snpEff.jar download -v GRCh38.99

#(To check all available database : java -jar snpEff.jar databases)

# Step 2: To annotate the input vcf file and create a web page summarizing the annotation results in a html page: "YOUR_VCF_FILE.html" (option -stats):

java -Xmx8g -jar snpEff.jar -v -stats YOUR_VCF_FILE.html GRCh38.99 YOUR_VCF_FILE > YOUR_VCF_FILE.annotated.vcf

#Example in my case:

java -Xmx8g -jar /Users/mstermo/snpEff/snpEff.jar -v -stats WES_ccRCC-chr1-22-X_snp.recalibrated.99.8.QD.QUAL.Inbreed.PPs.GQ.DP.split.PASS.mono.onlyAC.missing0.25.MAF.annotated.remove.vcf.gz.html GRCh38.99 WES_ccRCC-chr1-22-X_snp.recalibrated.99.8.QD.QUAL.Inbreed.PPs.GQ.DP.split.PASS.mono.onlyAC.missing0.25.MAF.annotated.remove.vcf.gz > WES_ccRCC-chr1-22-X_snp.recalibrated.99.8.QD.QUAL.Inbreed.PPs.GQ.DP.split.PASS.mono.onlyAC.missing0.25.MAF.annotated.remove.vcf.Snpeff.gz

#SnpEff produces three output files :

#the HTML file containing summary statistics about the variants and their annotations;
#an annotated VCF file; and
#a text file summarizing the number of variant types per gene.
