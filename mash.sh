outputf="/data/within_individual_selection/mash/"

mkdir $outputf

threads=25 # number of threads to use
kmer_size=21 # kmer size to use for minhash sketches
m=5 # Minimum copies of each k-mer required to pass noise filter for reads. You can estimate this by running genomescope first...
s=1000000 #Sketch size. Each sketch will have at most this many non-redundant min-hashes. Bigger is better but slower.
sg='100M' 
read_length=150 # length of reads in your input fastq files
kmer_max=5000 # ignore kmers represented more the kmer_max times (avoids issues with e.g. cpDNA)

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

for readf in $(find *.fastq.gz); do

	# count kmers in jellyfish, then make a histogram. NB: use zcat because jellyfish needs raw fastq, not .gz
	zcat $readf | jellyfish count -C -m $kmer_size -s $sg -t $threads -o $outputf$readf'.jf' /dev/fd/0
	jellyfish histo -t $threads $outputf$readf'.jf' > $readf'.histo'

	# run genomescope and put the output in a new directory
	id=`echo "$readf" | cut -d'.' -f1`
	genomescope.R $readf'.histo' $kmer_size $read_length $outputf$id $kmer_max verbose

done


# now we mash
gzs=$(find *.fastq.gz)
mash sketch -p $threads -m $m -k $kmer_size -s $s $gzs

# use mash paste to join the sketches together
fs=$(find *.msh)
mash paste raw $fs   

# verify it worked
mash info raw.msh

# calculate all pairwise distances
mash dist raw.msh raw.msh > distances.tab

Rscript --vanilla heatmap.r $outputf"distances.tab" $outputf"heatmap.pdf"