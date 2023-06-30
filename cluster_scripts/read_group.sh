INPUTDIR=/cluster/projects/pughlab/hereditary/external_data/TGL49_CHARM/LFS/LFS_WG/bams
basedir=/cluster/projects/pughlab/hereditary/projects/CHARM/LFS/read_group
shdir=$basedir/sh_scripts
outdir=$basedir/output

mkdir -p $outdir
mkdir -p $shdir

cd $INPUTDIR
ls *.bam > $shdir/bams
cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

echo -e "#!/bin/bash\n
module load samtools/1.10\n" > $shdir/script.sh

for bam in $(cat bams); do
echo -e "samtools view -H $INPUTDIR/${bam}.bam |grep '^@RG' | grep -o 'PU:[^[:space:]]*' > $outdir/${bam}.txt" >> $shdir/script.sh

done 

cd $shdir

sh script.sh
