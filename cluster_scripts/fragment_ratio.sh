INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
PDIR=/cluster/projects/pughlab/projects/CHARM/LFS/fragment_ratio
shdir=$PDIR/sh_scripts
outdir=$PDIR/output
DELFI=/cluster/projects/pughlab/bin/fragmentomics/v2/ratio

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *.bam > $shdir/bams

cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}.sh
#!/bin/bash
#
#$ -cwd

module load R/4.0.0
Rscript $DELFI/runFrag.R\
 --id $bam\
 --bamdir $INPUTDIR\
 --filters $DELFI/extdata/filters.hg38.rda\
 --gaps $DELFI/extdata/gaps.hg38.rda\
 --VNTRs $DELFI/extdata/VNTRs.hg38.rda\
 --tiles $DELFI/extdata/hg38_tiles.bed\
 --healthy $DELFI/extdata/healthy.median.hg38.rda\
 --outdir $outdir\
 --libdir $DELFI

EOF

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 -p himem --mem 40G $file

done
