
foldername=$1
cd $foldername
folder=`ls ./ | grep bam`
for bam in $folder; 
do 
	echo $bam
	out_index=`echo ${bam} | sed 's/.bam/.bam.bai/'`
	echo $out_index
	samtools index -b $bam $out_index
done; 


