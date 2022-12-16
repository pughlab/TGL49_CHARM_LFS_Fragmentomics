#!/bin/bash

griffin=/cluster/projects/pughlab/bin/Griffin/v0.1.0
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/griffin
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
input=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
outdir=$basedir/output/GC_correction
shdir=$basedir/sh_scripts/GC_correction

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

echo -e "$griffin/scripts/griffin_GC_counts.py \
--bam_file $input/${bam}.bam \
--bam_file_name $name \
--mapable_regions $griffin/Ref/repeat_masker.mapable.k50.Umap.hg38.bedGraph \
--ref_seq $ref \
--chrom_sizes $griffin/Ref/hg38.standard.chrom.sizes \
--out_dir $outdir \
--map_q 20 \
--size_range 15 500 \
--CPU 8\n" >> $shdir/${name}.sh

echo -e "$griffin/scripts/griffin_GC_bias.py \
--bam_file_name $name \
--mapable_name repeat_masker.mapable.k50.Umap.hg38 \
--genome_GC_frequency $griffin/Ref/genome_GC_frequency \
--out_dir $outdir/ \
--size_range 15 500" >> $shdir/${name}.sh

done

cd $shdir
ls *.sh > files
for file in $(cat files);do
sbatch -p all -c 8 --mem 8G -t 24:00:00 $file
done
