#!/bin/bash

base=/cluster/projects/pughlab/projects/CHARM/LFS/griffin
input=/cluster/projects/pughlab/projects/CHARM/LFS/griffin/sites
outdir=$base/output/sites
shdir=$base/sh_scripts/sites
griffin=/cluster/projects/pughlab/bin/Griffin/v0.1.0

mkdir -p $outdir
mkdir -p $shdir

cd $input
ls *.txt > $shdir/sites
cd $shdir

sed 's/....$//' sites > site
mv site sites

for site in $(cat sites); do
echo -e "#!/bin/bash
source activate base
conda activate griffin" > $shdir/${site}.sh

echo -e "$griffin/scripts/griffin_filter_sites.py \
-i $input/${site}.txt \
--name $site \
-o $outdir \
-m $griffin/Ref/k50.Umap.MultiTrackMappability.hg38.bw \
-c Chrom \
-s start \
-e end \
--strand_column NA \
--chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
--window_values -5000 5000 \
--targeted_panel FALSE \
--targeted_window_columns NA NA \
--threshold 0.95" >> $shdir/${site}.sh

done

cd $shdir
ls *.sh > files
for file in $(cat files); do
sbatch --mem 4G -t 1:00:00 $file
done
