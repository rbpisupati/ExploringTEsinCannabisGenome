# Given the file having divergence calculations generated for each of the contig
#from the pipeline potato_skin_Divergence_plots_for_contigs_RE.sh
#
# Input:	.pi file
#
# Run:		bash potatoSkin_editDivergence_Clusterwise.sh __.pi summary.csv contigs_longest_500.fa

# copied the table from summary.html and make tab delimited summary.csv file

# Output:	final.txt file (table with all the calculations)
# 		TE family wise contig ids
#		final.txt.TEfam.txt (final table with top hit from summary file and its %content in genome)
# use perl one liner to generate fasta files for each contig
#	perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file fasta.file

SCRIPT_PATH="/home/GMI/rahul.pisupati/mygit/ExploringTEsinCannabisGenome/scripts/"

## generates TEfamilies_CLwise_final_table.csv
Rscript ${SCRIPT_PATH}/02.01.Export_AveCovDiv_SlopeHalfLife_Clusterwise.R $1 $2


#  Give contigs fasta file to extract specific contigs to each family
perl ${SCRIPT_PATH}/02.02.classifyContigsTEfamilies.pl TEfamilies_CLwise_final_table.csv $3
