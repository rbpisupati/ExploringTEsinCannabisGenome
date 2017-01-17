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


allCLs=( $(cut -f1 $1| sort | uniq) )
clsLen=`cut -f1 $1| sort | uniq | wc -l`

#echo ${allCLs[3]}
#echo $clsLen

for(( i=0; i<$clsLen; i++ ));do
	clName=${allCLs[$i]}
	grep ^$clName $1 > temp.txt
	Rscript ~/Pipeline_makeDivergenceTable_CLwise_TEfam_fromSummaryFile/Export_AveCovDiv_SlopeHalfLife_Clusterwise.R temp.txt
	rPrint=`cat temp.fin`
	echo $clName $rPrint >> final.txt	
done

# adding the TE families to the Divergence table
sed -i 's/ /\t/g' final.txt
perl ~/Pipeline_makeDivergenceTable_CLwise_TEfam_fromSummaryFile/add_TEfamily_clusterDivergence_table.pl final.txt $2
sed -i 's/ /\t/g' final.txt.TEfam.txt

#  Give contigs fasta file to extract specific contigs to each family
perl ~/Pipeline_makeDivergenceTable_CLwise_TEfam_fromSummaryFile/classifyContigsTEfamilies.pl final.txt.TEfam.txt $3

