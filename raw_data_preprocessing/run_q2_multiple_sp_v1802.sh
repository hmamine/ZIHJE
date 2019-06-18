#!/bin/bash

# script for pre-processing 16S (V5-V7) amplicons sequencing using Qiime 2 v2.0

# originally written by M. Amine Hassani - ahassani@bot.uni-kiel.de - Nov2018

#SBATCH --job-name=qiime2_multiple
#SBATCH --mail-user=ahassani@bot.uni-kiel.de
#SBATCH --mail-type=ALL
#SBATCH --output=q2_multi.out
#SBATCH --error=q2_multi.err
#SBATCH --nodes=2
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60000
#SBATCH --time=240:00:00
#SBATCH --partition=angus

# loading the configuration file 
config_file=$1
source $config_file

# loading modules
module load qiime2-2018.2.0
module load intelmpi16.0.0 
export OMP_NUM_THREADS="$nbr_thr"

# clearing files
rm -fr "$path_to_data/output_*" "$path_to_data/mapping_files"

# making new directories
mkdir -p "$path_to_data/mapping_files"

# looping over the files

for l in "${list[@]}" 

do
	mkdir -p "$path_to_data/output_$l"   	 
	cp "$path_to_raw/sp$l/sample_metadata.tsv" "$path_to_data/mapping_files/sample_$l.tsv"
	mv "$path_to_raw/sp$l/sample_metadata.tsv" "$path_to_data/output_$l"	
	
	echo -e "\n 1-importing sequencing data from sp$l ..."
	
	qiime tools import --type EMPPairedEndSequences \
   	 --input-path "$path_to_raw/sp$l" \
   	 --output-path "$path_to_data/output_$l/emp-paired-end-sequences_$l.qza"
	
	echo "emp-paired-end-sequences_$l.qza saved in $path_to_data/output_$l" 
	
	mv "$path_to_data/output_$l/sample_metadata.tsv" "$path_to_raw/sp$l"

	echo -e "\n 2-demultiplexing sequences form sp$l ..."

	qiime demux emp-paired \
	  --m-barcodes-file "$path_to_raw/sp$l/sample_metadata.tsv" \
	  --m-barcodes-column BarcodeSequence \
	  --i-seqs "$path_to_data/output_$l/emp-paired-end-sequences_$l.qza" \
	  --o-per-sample-sequences "$path_to_data/output_$l/demux_$l" \
	  --p-rev-comp-mapping-barcodes
	
	echo -e "\n 3-sequence QC and generating feature table and representative sequences for sp$l ..."

	srun qiime dada2 denoise-paired \
	  --i-demultiplexed-seqs "$path_to_data/output_$l/demux_$l.qza" \
	  --p-trunc-len-f "$len" --p-trunc-len-r "$len" \
	  --p-n-threads "$nbr_thr" --p-n-reads-learn 9999 \
	  --o-table "$path_to_data/output_$l/table_$l" \
	  --o-representative-sequences "$path_to_data/output_$l/rep-seqs_$l"

	inargA+=( --i-tables "$path_to_data/output_$l/table_$l.qza" )
	inargB+=( --i-data "$path_to_data/output_$l/rep-seqs_$l.qza" )

done

mkdir -p "$path_to_data/output_merge"

echo -e "\n merging feature tables ..."

	qiime feature-table merge \
	"${inargA[@]}" \
	--o-merged-table "$path_to_data/output_merge/table_merge.qza"

echo -e "\n merging representative sequences ..."

	qiime feature-table merge-seqs \
	"${inargB[@]}" \
	--o-merged-data "$path_to_data/output_merge/rep-seqs_merge.qza"

echo -e "\n merging sample data files ..."

	( head -2 $path_to_data/mapping_files/sample_005.tsv && \
	tail -n +3 -q $path_to_data/mapping_files/sample_*.tsv ) \
	> $path_to_data/output_merge/sample_data_merge.tsv

	echo -e "sample_data_merge.tsv saved in $path_to_data/output_merge"

echo -e "\n multi-alignment of representative sequences ..."

	srun qiime alignment mafft \
	 --i-sequences "$path_to_data/output_merge/rep-seqs_merge.qza" \
	 --o-alignment "$path_to_data/output_merge/aligned-rep-seqs.qza" \
	 --p-n-threads "$nbr_thr" 

echo -e "\n filtering alignment ..."

	qiime alignment mask \
	  --i-alignment "$path_to_data/output_merge/aligned-rep-seqs.qza" \
	  --o-masked-alignment "$path_to_data/output_merge/masked-aligned-rep-seqs.qza" 

echo -e "\n generating phylogenetic tree ..."

	qiime phylogeny fasttree \
	  --i-alignment "$path_to_data/output_merge/masked-aligned-rep-seqs.qza" \
	  --o-tree "$path_to_data/output_merge/unrooted-tree.qza" 

echo -e "\n rooting the tree at midpoint ..."

	qiime phylogeny midpoint-root \
	  --i-tree "$path_to_data/output_merge/unrooted-tree.qza" \
	  --o-rooted-tree "$PWD/output_merge/rooted-tree.qza"

echo -e "\n taxonomic classification of representative sequences ..."

	mpirun -np 4 qiime feature-classifier classify-sklearn  --i-classifier "$class_path" \
	 --i-reads "$path_to_data/output_merge/rep-seqs_merge.qza" \
	 --o-classification "$path_to_data/output_merge/taxonomy.qza" --verbose \
	 --p-reads-per-batch 0  --p-n-jobs 4 --p-pre-dispatch 3*n_jobs

echo -e "\n exporting .qza files..."

qiime tools export "$path_to_data/output_merge/taxonomy.qza" --output-dir export
qiime tools export "$path_to_data/output_merge/rooted-tree.qza" --output-dir export
qiime tools export "$path_to_data/output_merge/table_merge.qza" --output-dir export

sed -i '1s/^/\#/' $path_to_data/export/taxonomy.tsv

biom add-metadata -i "$path_to_data/export/feature-table.biom" --output-fp "$path_to_data/export/table-with-taxonomy.biom" --observation-metadata-fp "$path_to_data/export/taxonomy.tsv" --sc-separated taxonomy

biom convert -i "$path_to_data/export/table-with-taxonomy.biom" -o "$path_to_data/export/table-with-taxonomy.txt" --to-tsv 

echo -e "\n archiving export ..."

tar czvf archive_export.tar.gz export/

echo -e "\n*** done ***\n"
