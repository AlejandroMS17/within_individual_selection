outputf="/data/within_individual_selection/"

mkdir $outputf

threads=25 # number of threads to use
kmer_size=31 # kmer size to use for minhash sketches
m=5 # Minimum copies of each k-mer required to pass noise filter for reads. You can estimate this by running genomescope first...
s=10000000 #Sketch size. Each sketch will have at most this many non-redundant min-hashes. Bigger is better but slower.
sg=5000000000 # G+Gcek from the manual: G = 450M (genome size), c = 30x, e = 1%, k = ~32, so sg ~ 6000000000
read_length=150 # length of reads in your input fastq files
kmer_max=500 # ignore kmers represented more the kmer_max times (avoids issues with e.g. cpDNA)

cd $outputf

for tree in RB5 RB7 RB98; do
	for branch in A B C D E F G H; do
		# this uses info in a tsv file to find all fastq files from one branch. Should be 12 files in most of our datasets
		# 12 files are forward and reverse reads, for each of 3 replicate samples, sequenced independently on two lanes (2*3*2=12)
		fastqs=$(awk -v tree=$tree -v branch=$branch -F"\t" '$1 == tree && $2 == branch { print "/data/raw_data/"tree"/"$3"/"$4 }' sequencing_file_details.tsv)
		prefix=$tree"_"$branch
		# cat the fastqs and pipe to mash (stdin input is a - in -r)
		# NB - cat seems so slow thta mash only runs on one CPU
		#cat $fastqs | mash sketch -k 21 -r - -p $threads -m $m -k $kmer_size -s $s -o $prefix
		# cat first, mash second
		cat $fastqs > $prefix".fastq.gz"
	done
done

# first we do genomescoepe to help choose parameters for the mash analysis and QC the data

#for readf in $(find *.fastq.gz); do
#	# count kmers in jellyfish, then make a histogram. NB: use zcat because jellyfish needs raw fastq, not .gz
#	zcat $readf | jellyfish count -C -m $kmer_size -s $sg -t $threads -o $outputf$readf'.jf' /dev/fd/0
#	jellyfish histo -t $threads $outputf$readf'.jf' > $readf'.histo'
#	# run genomescope and put the output in a new directory
#	id=`echo "$readf" | cut -d'.' -f1`
#	genomescope.R $readf'.histo' $kmer_size $read_length $outputf$id $kmer_max verbose
#done


# now we mash, genome size estimated at 450MB
# do it with different error cutoffs and kmer sizes (m and k)
gzs=$(find *.fastq.gz)
for m in 3 4 5; do
	for kmer_size in 30 31 32; do
		id="m"$m"_k"$kmer_size
		mash sketch -p $threads -m $m -k $kmer_size -s $s -g 450000000 -o $id $gzs 
		mash info $id".msh"
		mash dist $id".msh" $id".msh" > $id".tab"
		Rscript --vanilla heatmap.r $outputf$id".tab" $outputf$id"_heatmap.pdf"
	done
done

rm *.fastq.gz
rm *.msh
