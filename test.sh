# Put the test file and model somewhere, might be dropbox I guess?

# Download testfile

#nextflow run main.nf --inputdir testin/ --matcher "*_R{1,2}_001.fastq.gz" --outputdir testout --manifest testin/manifest.tsv -profile lipid,conda -resume
#NXF_VER=20.01.0.5264 nextflow run main.nf -profile lipid,conda,test -resume
nextflow run main.nf -profile lipid,conda,test -resume
