# Put the test file and model somewhere, might be dropbox I guess?

# Download testfile

nextflow run main.nf --inputdir testin/ --outputdir testout --manifest testin/manifest.tsv -profile lipid,conda -resume
