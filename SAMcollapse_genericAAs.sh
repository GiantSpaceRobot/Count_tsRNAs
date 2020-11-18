#!/bin/bash

### Remove tempDir if it exists
if [ -d "tempDir" ]; then
	echo "Deleting tempDir"
fi

### Define variables
cores=$2
chunks=500
if (( $cores > $chunks )); then
	# make sure cores is not set higher than no. of files after split
	cores=$chunks
fi
fileLen=$(< "$1" wc -l)
division1=$((fileLen/chunks))
division=$((division1 + 1))
myFile="tempFile"
mkdir -p tempDir

### Remove header
grep ^@ $1 > tempDir/myHeader.txt &
grep -v ^@ $1 > tempDir/mySAM.sam &
wait

### 
fileToCollapse=tempDir/mySAM.sam

### Split file
split -l $division $fileToCollapse tempDir/splitFile_

### Gather first and last read from every split file and add to separate file. Remove these reads from the split files.
for i in tempDir/splitFile_*; do
	base=$(basename $i)
	first=$(awk 'NR==1' $i | awk '{print $1}') 
	echo $first >> tempDir/${myFile}_HeadsAndTails.txt
	last=$(awk 'END{print}' $i | awk '{print $1}') 
	echo $last >> tempDir/${myFile}_HeadsAndTails.txt
	#grep -P -E -v "$first\t|$last\t" $i > tempDir/edit_${base}  #remove all instances of first and last read from file  
	#grep -P -v "$first\t" $i > tempDir/grepFile_${base}.txt
	#grep -P -v "$last\t" tempDir/grepFile_${base}.txt > tempDir/edit_${base}
done

sort tempDir/${myFile}_HeadsAndTails.txt | uniq > tempDir/${myFile}_HeadsAndTails_uniq.txt #remove duplicates
sed -i 's/$/\t/' tempDir/${myFile}_HeadsAndTails_uniq.txt # Add tab to end of every line to match pattern exactly
grep -f tempDir/${myFile}_HeadsAndTails_uniq.txt $fileToCollapse > tempDir/edit_heads-and-tails #grep all patterns from the heads/tails file

for i in tempDir/splitFile_*; do
	base=$(basename $i)
	grep -v -f  tempDir/${myFile}_HeadsAndTails_uniq.txt $i > tempDir/edit_${base}
done


### Run SAMcollapse2.py. This loop will only run $cores processes at once
COUNTER=0
for i in tempDir/edit_*; 
do
	base=$(basename $i)
	python SAMcollapse-genericAAs.py $i ${fileToCollapse}_${base} & 
	numjobs=($(jobs | wc -l))
	echo Running job number ${COUNTER} of ${chunks}... 
	COUNTER=$[$COUNTER + 1]
	#echo Initialising job number $numjobs
	#echo There are ${#background[@]} jobs now
    while (( $numjobs == $cores )); do
    	echo There are $numjobs jobs now. Waiting for jobs to finish...
		numjobs=($(jobs | wc -l))
		sleep 2 #Enter next loop iteration
	done
done
wait

### Concatenate results
cat tempDir/myHeader.txt ${fileToCollapse}*_edit_* > Collapsed.sam
grep -v ^@ Collapsed.sam | grep chr > $3
### Sort SAM to BAM
#samtools sort Collapsed.sam > Collapsed.bam

rm -rf tempDir/

