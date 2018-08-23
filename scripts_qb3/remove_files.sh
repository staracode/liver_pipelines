rm $DIR_NAME/fastq/*fastqc.zip
rm $DIR_NAME/filtered/*fastqc.zip
rm -r $DIR_NAME/algn/*STARtmp/
gzip $DIR_NAME/algn/*
