#!/bin/bash

function SAMcollapse () {
	echo "Collapsing SAM file..."
	### How many chunks to create:
	readNumber=$(grep -v ^@ $1 | awk '{print $1}' | uniq | wc -l)
	if (( $readNumber > 100000 )); then
		chunksRaw=$(( readNumber / 10000 )) # 10,000 lines per SAM chunk
		chunks=$( echo $chunksRaw | awk '{print int($1+0.5)}' )
		echo "Over 100,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	elif (( readNumber>=10000 && readNumber<=100000 )); then
		chunksRaw=$(( readNumber / 5000 )) # 5000 lines per SAM chunk
		chunks=$( echo $chunksRaw | awk '{print int($1+0.5)}' )
		echo "Between 10,000 and 100,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	elif (( readNumber>=1000 && readNumber<=10000 )); then
		chunks=10 #
		echo "Between 1,000 and 10,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	else
		chunks=1
		echo "Less than 1,000 unique reads in SAM file. Splitting SAM into $chunks files..."
	fi
	### Define variables
	threads_available_for_chunks=30
	if (( $threads_available_for_chunks > $chunks )); then
		# make sure the number of threads is not set higher than no. of files (chunks) after split
		threads_available_for_chunks=$chunks
	fi
	echo "
		Threads: $threads_available_for_chunks
		# Split SAM files: $chunks"
	fileLen=$(< "$1" wc -l)
	division1=$((fileLen/chunks))
	division=$((division1 + 1))
	myFile="tempFile"
	mkdir -p tempDir
	mkdir -p tempDir/tRNAsAlmostMapped

	### Remove header
	grep ^@ $1 > tempDir/myHeader.txt &
	grep -v ^@ $1 > tempDir/mySAM.sam &
	wait
	### 
	fileToCollapse=tempDir/mySAM.sam
	### Split file
	echo "Splitting SAM..."
	split -l $division $fileToCollapse tempDir/splitFile_
	### Gather first and last read from every split file and add to separate file. Remove these reads from the split files.
	echo "Gathering the names of the first and last reads from every SAM chunk..."
	for i in tempDir/splitFile_*; do
		base=$(basename $i)
		first=$(awk 'NR==1' $i | awk '{print $1}') 
		echo $first >> tempDir/${myFile}_HeadsAndTails.txt
		last=$(awk 'END{print}' $i | awk '{print $1}') 
		echo $last >> tempDir/${myFile}_HeadsAndTails.txt
	done
	echo "Gathering unique set of read names from first/last read names..."
	sort tempDir/${myFile}_HeadsAndTails.txt | uniq \
		> tempDir/${myFile}_HeadsAndTails_uniq.txt #remove duplicates
	sed -i 's/$/\t/' tempDir/${myFile}_HeadsAndTails_uniq.txt # Add tab to end of every line to match pattern exactly
	grep -f tempDir/${myFile}_HeadsAndTails_uniq.txt $fileToCollapse \
		> tempDir/edit_heads-and-tails #grep all patterns from the heads/tails file
	echo "Extracting alignments for first/last reads from all files..."
	for i in tempDir/splitFile_*; do
		base=$(basename $i)
		grep -v -f  tempDir/${myFile}_HeadsAndTails_uniq.txt $i \
			> tempDir/edit_${base}
	done
	### Run SAMcollapse.py. This loop will only run $threads_available_for_chunks processes at once
	COUNTER=1
	if (( $chunks == 1)); then # If there are so few reads in input that only one SAM chunk is present
		chunksDiv=1
	else # This is the most common situation: Very large SAM has been split into many chunks, divide the chunks by 10 to allow printing of progress
		chunksDiv=$((chunks/10))
	fi
	echo "Collapsing every chunk of SAM..."
	for i in tempDir/edit_*; 
	do
		base=$(basename $i)
		if (( $chunks > 1)); then
			python2 SAMcollapse.py \
				$i \
				${fileToCollapse}_${base} \
				>> tempDir/tRNAsAlmostMapped/collapsed-reads.txt & 
			if (( $COUNTER % $chunksDiv == 0 )); then # If SAM is split into 100 chunks, print progress when the counter reaches 10, 20, 30... 
				echo "Started job $COUNTER of $chunks" # i.e. print progress every 10%
			fi
			numjobs=($(jobs | wc -l))
			COUNTER=$[$COUNTER + 1]
			while (( $numjobs == $threads_available_for_chunks )); do
				numjobs=($(jobs | wc -l))
				sleep 2 #Enter next loop iteration
			done
		else
			python2 SAMcollapse.py \
				$i \
				${fileToCollapse}_${base} \
				>> tempDir/tRNAsAlmostMapped/collapsed-reads.txt & 
			echo "Started job $COUNTER of $chunks" 
		fi
	done
	wait
	readsCollapsedSpecies=$(awk '{split($0,a," "); sum += a[1]} END {print sum}' tempDir/tRNAsAlmostMapped/collapsed-reads.txt)
	readsCollapsedGroup=$(awk '{split($0,a," "); sum += a[2]} END {print sum}' tempDir/tRNAsAlmostMapped/collapsed-reads.txt)
	echo -e "SAM collapse results:\n\t$readsCollapsedSpecies reads collapsed at the tRNA species level (e.g. 2 gene copies of ProCCG)\n\t$readsCollapsedGroup reads collapsed at the tRNA group level (e.g. ProCCG and ProAAG)"

	### Concatenate results
	echo "Gathering reads that were mapped to similar tRNAs..."
	echo -e "tRNA.group\tread.start\tread.end.approx\tread.name" \
		> tempDir/tRNAsAlmostMapped/tRNAs-almost-mapped.txt
	cat tempDir/*tRNAs-almost-mapped* | sort \
		>> $filename/STAR-output/tRNAs-almost-mapped.txt
	#mkdir tempDir/tRNAsAlmostMapped
	mv tempDir/*tRNAs-almost-mapped* tempDir/tRNAsAlmostMapped/
	echo "Concatenating SAM header with collapsed files..."
	cat tempDir/myHeader.txt ${fileToCollapse}*_edit_* \
		> $filename/STAR-output/Collapsed.sam
	rm -rf tempDir/ # Remove temp directory
	echo "Finished collapsing SAM file"
}


for i in $1/*; do
	my_file=$(basename $i)
	substr=""
	case $my_file in *.fq.gz) filename=$(echo "${my_file/.fq.gz/$substr}");; esac # If fq.gz suffix in filename, remove
	echo $filename
	mkdir -p $filename
	mkdir -p $filename/FastQC
	mkdir -p $filename/Trimmed_data
	mkdir -p $filename/STAR-output
	mkdir -p $filename/Results
	if [[ $i == *".gz"* ]]; then
		suffix="fq.gz"
		STARparam="--readFilesCommand zcat"
	else
		suffix="fq"
		STARparam=""
	fi
	echo "Trimming..."
	trim_galore \
		--stringency 10 \
		-o $filename/Trimmed_data/ \
		--fastqc_args "--outdir $filename/FastQC/" \
		$i
	echo "Running STAR..."
	STAR \
	--runThreadN 15 \
	--genomeDir /home/paul/Documents/Pipelines/tsRNAsearch/DBs/species_index/mouse-ncRNAs/ \
	--readFilesIn $filename/Trimmed_data/${filename}_trimmed.fq.gz \
	--outFileNamePrefix $filename/STAR-output/ \
	--outSAMattributes AS nM HI NH \
	--outFilterMultimapScoreRange 0 \
	$STARparam
	### Generate standard collapsed SAM file and generate depth files for plotting
	echo "SAM collapse..."
	SAMcollapse $filename/STAR-output/Aligned.out.sam #Collapse reads aligned to the same tRNA species 
	echo "Processing files..."
	mv $filename/STAR-output/Collapsed.sam $filename/STAR-output/aligned_tRNAdb.sam # match hisat2 naming convention
	grep -v "tRNA.group" $filename/STAR-output/tRNAs-almost-mapped.txt \
		| awk '{print $1}' \
		| uniq -c \
		| awk '{print $2"\t"$1}' \
		> $filename/STAR-output/tRNAs-almost-mapped.count
	grep ^@ $filename/STAR-output/aligned_tRNAdb.sam > $filename/STAR-output/SamHeader.sam &
	grep ENS $filename/STAR-output/aligned_tRNAdb.sam > $filename/STAR-output/ncRNAs.sam &
	grep -v ENS $filename/STAR-output/aligned_tRNAdb.sam > $filename/STAR-output/tsRNAs_aligned.sam &
	wait
	cat $filename/STAR-output/SamHeader.sam $filename/STAR-output/ncRNAs.sam \
		> $filename/STAR-output/ncRNAs_aligned.sam
	wait
	echo "Converting to BAM and indexing..."
	samtools view -bS $filename/STAR-output/tsRNAs_aligned.sam \
		| samtools sort \
		> $filename/STAR-output/tRNA-alignment_accepted_hits.bam
	samtools index $filename/STAR-output/tRNA-alignment_accepted_hits.bam &
	#mv $filename/STAR-output/unmapped.fastq $filename/STAR-output/$myFile
	### Move ncRNA BAM to directory
	samtools view -bS $filename/STAR-output/ncRNAs.sam \
		| samtools sort \
		> $filename/STAR-output/ncRNA-alignment_accepted_hits.bam
	samtools index $filename/STAR-output/ncRNA-alignment_accepted_hits.bam &
	#rm \
	#	$filename/STAR-output/SamHeader.sam \
	#	$filename/STAR-output/ncRNAs.sam \
	#	$filename/STAR-output/aligned_tRNAdb.sam &|
	### Generate counts for ncRNAs
	echo "Counting using featureCount..."
	featureCounts \
		-T 15 \
		-a /home/paul/Documents/Pipelines/tsRNAsearch/DBs/mouse_ncRNAs_relative_cdhit.gtf \
		-o $filename/Results/ncRNA-alignment.fcount \
		$filename/STAR-output/ncRNA-alignment_accepted_hits.bam
	grep -v featureCounts $filename/Results/ncRNA-alignment.fcount \
		| grep -v ^Geneid \
		| awk -v OFS='\t' '{print $1, $7}' \
		> $filename/Results/ncRNA-alignment.count
	### Generate counts for tRNAs
	featureCounts \
		-T 15 \
		-a /home/paul/Documents/Pipelines/tsRNAsearch/DBs/mouse_tRNAs_relative_cdhit.gtf \
		-o $filename/Results/tRNA-alignment.fcount \
		$filename/STAR-output/tRNA-alignment_accepted_hits.bam
	grep -v featureCounts $filename/Results/tRNA-alignment.fcount \
		| grep -v ^Geneid \
		| awk -v OFS='\t' '{print $1, $7}' \
		> $filename/Results/tRNA-alignment.count
	### Get number of mapped reads
	#echo "Get total number of reads mapped"
	echo "Getting number of mapped reads..."
	cat $filename/Results/*.count | grep -v ^__ | sort -k1,1 \
		> $filename/Results/$filename.all-features.count
	sed -i '1s/^/Features\t'"$filename"'\n/' $filename/Results/$filename.all-features.count # Add column headers
	direct_mapped=$(awk '{sum+=$2} END{print sum;}' $filename/Results/$filename.all-features.count)
	my_var=$(wc -l $filename/STAR-output/tRNAs-almost-mapped.txt | cut -f1 -d' ') # Count number of multimapping reads salvaged
	multimapped_and_reclaimed=$(( my_var - 1 )) # Reduce by 1 to remove header
	mapped=$(( direct_mapped + multimapped_and_reclaimed ))
	### Generate depth files
	echo "Generating depth file..."
	samtools depth \
		-d 100000000 \
		-aa $filename/STAR-output/tRNA-alignment_accepted_hits.bam \
		> $filename/Results/tRNA-alignment_accepted_hits.depth   # A lot faster than bedtools genomecov
	cp $filename/Results/tRNA-alignment_accepted_hits.depth $filename/Results/tRNA-alignment_accepted_hits_raw.depth
	### Normalise by reads per million (RPM)
	python2 /home/paul/Documents/Pipelines/tsRNAsearch/bin/Depth-to-Depth_RPM.py \
		$filename/Results/tRNA-alignment_accepted_hits_raw.depth \
		$mapped \
		$filename/Results/tRNA-alignment_accepted_hits_RPM.depth
	
	### Generate counts of all possible tRNA fragments
	python3 SAM-to-tsRNA-count-v2.py /home/paul/Documents/Pipelines/tsRNAsearch/DBs/mouse_tRNAs_relative_cdhit.gtf $filename/STAR-output/tsRNAs_aligned.sam $filename/Results/${filename}_every-coordinate.tsv &

	### Clean up dir (remove big files)
	#rm $filename/STAR-output/*sam 
done

wait

mkdir -p EveryCoordinate
find -iname '*every-coordinate.tsv' -exec cp {} EveryCoordinate/ \;

exit 1
 
### Process tsRNA counts to RPM text files and plots
mkdir -p Results
cd Results
Rscript /home/paul/Dropbox/Work/scripts/R/tsRNA-Conditional-Differences_25-9-20.R ../MZ_ANG0_1.tsv ../MZ_ANG0_2.tsv ../MZ_ANG0_3.tsv ../MZ_ANG500_1.tsv ../MZ_ANG500_2.tsv ../MZ_ANG500_3.tsv MZ_ANG0-vs-ANG500

