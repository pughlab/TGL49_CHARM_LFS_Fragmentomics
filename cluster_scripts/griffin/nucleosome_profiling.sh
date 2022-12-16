#!/bin/bash

analysis=LFS
CPU=1
mem=4G

griffin=/cluster/projects/pughlab/bin/Griffin/v0.1.0
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/griffin
sites=$griffin/site_configs/${analysis}_sites.yaml
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
input=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
counts=$basedir/output/GC_correction/repeat_masker.mapable.k50.Umap.hg38/GC_bias
outdir=${basedir}_all/output/nucleosome_profiling/$analysis
shdir=${basedir}_all/sh_scripts/nucleosome_profiling/$analysis

mkdir -p $outdir
mkdir -p $shdir

cd $input
ls *bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
mv bam bams

for bam in $(cat bams);do

name=$(echo $bam | sed 's/\_WG.*/_WG/')
echo $bam
echo $name
echo -e "#!/bin/bash\n
source activate base\n
conda activate griffin\n" > $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_calc_coverage.py \
--sample_name $name \
--bam $input/${bam}.bam \
--GC_bias $counts/${name}.GC_bias.txt \
--background_normalization None \
--sites_yaml $sites \
--reference_genome $ref \
--results_dir $outdir \
--chrom_column Chrom \
--chroms chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 \
--norm_window -5000 5000 \
--plot_window -1000 1000 \
--fragment_length 165 \
--step 15 \
--size_range 35 220 \
--map_quality 20 \
--strand_column Strand \
--individual False \
--smoothing True \
--num_sites none \
--sort_by none \
--ascending none \
--cpu $CPU" >> $shdir/${name}.sh

done

cd $shdir
ls *.sh > files
for file in $(cat files);do
sbatch -c $CPU --mem $mem -t 24:00:00 $file
done
