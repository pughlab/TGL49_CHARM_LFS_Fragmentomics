INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
PDIR=/cluster/projects/pughlab/projects/CHARM/LFS/telomerehunter
shdir=$PDIR/sh_scripts
outdir=$PDIR/output

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

module load python2/2.7.15
module load R/3.3.0

telomerehunter\
 -ibt $INPUTDIR/${bam}.bam\
 -o $outdir\
 -p $bam\
 -b /cluster/projects/pughlab/bin/TelomereHunter/hg38_chr_bands.txt\
 -rt 6\
 -con

EOF

done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 24:00:00 --mem 8G $file
done
