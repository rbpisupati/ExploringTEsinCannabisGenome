#!/bin/bash
#PBS -S /bin/bash
#PBS -N divPlot_REcontig
#PBS -P nordborg_common
#PBS -l select=1:ncpus=8:mem=80gb
#PBS -l walltime=48:00:00
#PBS -o logs.oexploreRE
#PBS -e logs.oexploreRE

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
# bash ~/SCRIPT_PATH/potato_skin_Divergence_plots_for_contigs_RE.sh contigs ~/finola/sativa_finola_SRR351929_1.fastq ~/finola/sativa_finola_SRR351929_2.fastq SRR351929

cd $PBS_O_WORKDIR
module load BWA/0.7.13-intel-2016a
module load SAMtools/1.3.1-intel-2016a
module load picard
tempFol="/lustre/scratch/users/rahul.pisupati/tempFiles/"
nt=8
SCRIPT_PATH="/home/GMI/rahul.pisupati/mygit/ExploringTEsinCannabisGenome/scripts/"  ## gives the path

REcontigs="Repeat_contigs_RE_filtered_min500bp.fa"
sraid=`ls *fastq | cut -f1 -d "_" | uniq -d | head -n 1`

# Generating pileup file for contigs and mapping to raw reads
#bwa index $REcontigs
#bwa mem -t $nt $REcontigs ${sraid}_1.fastq ${sraid}_2.fastq > ${sraid}.sam
#samtools view -b -o ${sraid}.aligned.bam -S ${sraid}.sam
#samtools sort --threads $nt -o ${sraid}.sorted.bam ${sraid}.aligned.bam
#java -Djava.io.tmpdir=$tempFol -jar ${EBROOTPICARD}/picard.jar MarkDuplicates INPUT=${sraid}.sorted.bam OUTPUT=${sraid}.dedup.bam METRICS_FILE=${sraid}.metrics
#java -Djava.io.tmpdir=$tempFol -jar ${EBROOTPICARD}/picard.jar AddOrReplaceReadGroups INPUT=${sraid}.dedup.bam OUTPUT=${sraid}.modified.bam ID=$sraid LB=$sraid PL=illumina PU=none SM=$sraid
#samtools index ${sraid}.modified.bam

#if [ -f ${sraid}.modified.bam ]
#then
#        rm ${sraid}.sam ${sraid}.aligned.bam ${sraid}.sorted.bam ${sraid}.dedup.bam ${sraid}.metrics
#fi

# Changing quality scores and removing INDELs in the pileup before giving that to pileup
#samtools mpileup ${sraid}.modified.bam > ${sraid}.pileup
module load BioPerl/1.6.924-foss-2015a-Perl-5.20.3

perl ${SCRIPT_PATH}/01.01.changeQualScore_to_illumina_remove_INDELs_pileup.pl ${sraid}.pileup ${sraid}.modified.pileup
perl ${SCRIPT_PATH}/01.02.subSample_100cov_pileup_based_on_allelFractions.pl ${sraid}.modified.pileup ${sraid}.subcov.modified.pileup
## get the average coverage
awk ' NF > 0 {totcov[$1] = totcov[$1] + $4; count[$1] = count[$1] + 1; } END { for (contig in totcov) print contig "\t" totcov[contig]/count[contig]; }' ${sraid}.pileup | sed 's/Contig/\tContig/' > ${sraid}.avg.cov.txt

if [ -f ${sraid}.subcov.modified.pileup ]
then
        rm ${sraid}.pileup ${sraid}.modified.pileup
fi

# Running Popoolation on the subsampled pileup file
# Calculating divergence and D values
perl ~/Softwares/popoolation_1.2.2/Variance-sliding.pl --measure pi --input ${sraid}.subcov.modified.pileup --output ${sraid}.subcov100.pi --pool-size 1000 --min-count 2 --min-covered-fraction 0
sed -i 's/Contig/\tContig/' ${sraid}.subcov100.pi
#perl ~/Apps/popoolation_1.2.2/Variance-sliding.pl --measure D --input ${4}.subcov100.pileup --output ${4}.subcov100.D --pool-size 1000 --min-count 2 --min-covered-fraction 0
#sed -i 's/Contig/\tContig/' ${4}.subcov100.D

#avgDepth=`cat ${sraid}_avgDepth.txt`
perl ${SCRIPT_PATH}/01.03.join_avg_cover_pi.pl  ${sraid}.subcov100.pi ${sraid}.avg.cov.txt ${avgDepth} ${sraid}.subcov100.mod.pi
#perl ${SCRIPT_PATH}/01.03.join_avg_cover_pi.pl  ${sraid}.subcov100.D ${sraid}.avg.cov.txt ${avgDepth} ${sraid}.subcov100.mod.D

if [ -f ${sraid}.subcov100.mod.pi ]
then
        rm ${sraid}.subcov100.pi
fi

# Extract the reads at a position with hightest coverage from the sam file
#mkdir contig_positions
#cd contig_positions
#bash ~/SCRIPT_PATH/extract_Reads_atMaxPosition_Contig.sh ../${4}.sam
#cd ../
# Run Aaron's tool to give the phylogeny between the elements foreach contig
# Removing the sam and bam files
#rm  ${4}.sam  ${4}.bam
