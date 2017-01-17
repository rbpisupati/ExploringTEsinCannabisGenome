# Run the script from assembly folder of RepeatExplorer
# First input the contigs generated from the RE
# Second the raw reads fastq files or trimmed files
#       Froward read then reverse
# Fourth input as the SRA id

# Requirements:
#  large_file_fasta_min_length_histogram.pl
#  changeQualScore_to_illumina_remove_INDELs_pileup.pl
#  subSample_100cov_pileup_based_on_allelFractions.pl
#  join_avg_cover_pi.pl
#  extract_Reads_atMaxPosition_Contig.sh
#  bwa , samtools (also old for pileup) installed
#  popoolation at ~/Apps/popoolation_1.2.2


#Example run for finola 
# bash ~/cannabis/rahul/potato_skin_Divergence_plots_for_contigs_RE.sh contigs ~/cannabis/rahul/finola/sativa_finola_SRR351929_1.fastq ~/cannabis/rahul/finola/sativa_finola_SRR351929_2.fastq SRR351929

# Filtering contigs with min length 500
perl ~/cannabis/rahul/large_file_fasta_min_length_histogram.pl $1 500 > hist.contigs

# Generating pileup file for contigs and mapping to raw reads
bwa index ${1}.longest.min500.fa
bwa mem -t 8 ${1}.longest.min500.fa $2 $3 > ${4}.sam
samtools view -b -o ${4}.bam -S ${4}.sam
samtools sort ${4}.bam ${4}.sorted
samtools index ${4}.sorted.bam
~/Apps/samtools-0.1.15/samtools pileup ${4}.sorted.bam > ${4}.pileup

# Changing quality scores and removing INDELs in the pileup before giving that to pileup
perl ~/cannabis/rahul/changeQualScore_to_illumina_remove_INDELs_pileup.pl ${4}.pileup ${4}.modified.pileup
perl ~/cannabis/rahul/subSample_100cov_pileup_based_on_allelFractions.pl ${4}.modified.pileup ${4}.subcov100.pileup

# Running Popoolation on the subsampled pileup file
# Calculating divergence and D values
perl ~/Apps/popoolation_1.2.2/Variance-sliding.pl --measure pi --input ${4}.subcov100.pileup --output ${4}.subcov100.pi --pool-size 1000 --min-count 2 --min-covered-fraction 0
sed -i 's/Contig/\tContig/' ${4}.subcov100.pi

perl ~/Apps/popoolation_1.2.2/Variance-sliding.pl --measure D --input ${4}.subcov100.pileup --output ${4}.subcov100.D --pool-size 1000 --min-count 2 --min-covered-fraction 0
sed -i 's/Contig/\tContig/' ${4}.subcov100.D

# To get average coverage in pileup
awk ' NF > 0 {totcov[$1] = totcov[$1] + $4; count[$1] = count[$1] + 1; } END { for (contig in totcov) print contig "\t" totcov[contig]/count[contig]; }' ${4}.modified.pileup > ${4}.avg.cov
sed -i 's/Contig/\tContig/'  ${4}.avg.cov
perl join_avg_cover_pi.pl  ${4}.subcov100.pi ${4}.avg.cov ${4}.subcov100.mod.pi
perl join_avg_cover_pi.pl  ${4}.subcov100.D ${4}.avg.cov ${4}.subcov100.mod.D

mv ${4}.subcov100.mod.pi ${4}.subcov100.pi
mv ${4}.subcov100.mod.D ${4}.subcov100.D

# Extract the reads at a position with hightest coverage from the sam file

mkdir contig_positions
cd contig_positions
bash ~/cannabis/rahul/extract_Reads_atMaxPosition_Contig.sh ../${4}.sam
cd ../

# Run Aaron's tool to give the phylogeny between the elements foreach contig


# Removing the sam and bam files
rm  ${4}.sam  ${4}.bam
