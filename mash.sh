outputf="/data/within_individual_selection/mash/"

mkdir $outputf

threads=5 # number of threads to use
kmer_size=32 # kmer size to use for minhash sketches
m=5 # Minimum copies of each k-mer required to pass noise filter for reads. You can estimate this by running genomescope first...
s=1000000 #Sketch size. Each sketch will have at most this many non-redundant min-hashes. Bigger is better but slower.

cd $outputf

for tree in RB5 RB7 RB98; do

	for branch in A B C D E F G H; do

		# this uses info in a tsv file to find all fastq files from one branch. Should be 12 files in most of our datasets
		# 12 files are forward and reverse reads, for each of 3 replicate samples, sequenced independently on two lanes (2*3*2=12)
		fastqs=$(awk -v tree=$tree -v branch=$branch -F"\t" '$1 == tree && $2 == branch { print "/data/raw_data"tree"/"$3"/"$4 }' sequencing_file_details.tsv)

		prefix=$tree"_"$branch

		# cat the fastqs and pipe to mash (stdin input is a - in -r)
		cat $fastqs | mash sketch -k 21 -r - -p $threads -m $m -k $kmer_size -s $s -o $prefix

	done


done

# use mash paste to join the sketches together
fs=$(find *.msh)
mash paste raw $fs   

# verify it worked
mash info raw.msh

# calculate all pairwise distances
mash dist raw.msh raw.msh > distances.tab

Rscript --vanilla heatmap.r $outputf"distances.tab" $outputf"heatmap.pdf"