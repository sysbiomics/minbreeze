# Put the test file and model somewhere, might be dropbox I guess?

# Download testfile

nextflow run main.nf --inputdir testin/ --outputdir testout -profile lipid,conda -resume
