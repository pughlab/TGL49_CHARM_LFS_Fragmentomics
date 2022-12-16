INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
basedir=/cluster/projects/pughlab/projects/CHARM/LFS/peak_reads
ref=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
picard_dir=/cluster/tools/software/picard/2.10.9

site=DHS

bed=/cluster/projects/pughlab/projects/CHARM/LFS/peak_reads/${site}_peaks_hg38.bed
outdir=$basedir/output/$site
shdir=$basedir/sh_scripts/$site

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *.bam > $shdir/bams
cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
### Load in the modules
echo -e "#!/bin/bash\n
module load samtools/1.10
module load picard/2.10.9
module load bedtools/2.27.1\n" > $shdir/${bam}.sh

### Make bam file of only site regions
echo -e "samtools view -hb -f 0x2 -L $bed $INPUTDIR/${bam}.bam > $outdir/${bam}_site.bam
samtools index $outdir/${bam}_site.bam\n" >> $shdir/${bam}.sh

### Remove duplicates
echo -e "java -jar $picard_dir/picard.jar MarkDuplicates \
I=$outdir/${bam}_site.bam \
O=$outdir/${bam}_deduped.bam \
M=$outdir/${bam}_metrics.txt \
TMP_DIR=$outdir \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true

samtools sort $outdir/${bam}_deduped.bam -o $outdir/${bam}_deduped_sorted.bam
samtools index $outdir/${bam}_deduped_sorted.bam\n" >> $shdir/${bam}.sh

### Get Insert sizes
echo -e "java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
I=$outdir/${bam}_deduped_sorted.bam \
O=$outdir/${bam}_picard.txt \
H=$outdir/${bam}.pdf \
M=0 \
W=600\n" >> $shdir/${bam}.sh

### Remove intermediate files
echo -e "if [[ -f "$outdir/${bam}_picard.txt" ]]
then
echo -e \"Output completed sucessfully\"
rm $outdir/${bam}*.bam*
rm $outdir/${bam}_metrics.txt
else
echo -e \"Errors in run\"
fi" >> $shdir/${bam}.sh
done 

cd $shdir

ls *.sh > files
for file in $(cat files);do
sbatch -c 1 -t 10:00:00 -p all --mem 24G $file
done
