# Download testfile

NXF_VER=22.04.0 nextflow run main.nf -profile test,sysbiome,singularity \
    -resume #\
    #-with-report nextflow.html \
    #-with-timeline timeline.html \
    #-with-trace logs/trace.txt \
    #-with-dag dag.png
