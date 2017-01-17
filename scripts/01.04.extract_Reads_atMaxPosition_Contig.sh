#  Input to the script is the SAM file generatedafter aligning reads to contigs (RE)
# bash max_position_contigs.sh _____.sam
# Output is bunch of files (each for a contig) containing the reads at a position of max depth 

grep -v '@' $1 |cut -f3,4|sort -k1,1n|uniq -c > sorted

contigs=(`awk '{print $2}' sorted | sort | uniq`)
length=${#contigs[@]}

for((i=0;i<$length;i++));do
        reqContig=`grep ${contigs[$i]} sorted |sort -k 1,1nr | head -n 1 |awk '{print $2}'`
        reqPos=`grep ${contigs[$i]} sorted |sort -k 1,1nr | head -n 1 |awk '{print $3}'`
        awk -v con=$reqContig -v pos=$reqPos '$5 == "60" && $3 == con && $4 == pos {print ">" $1 "\n" $10}' $1 > ${reqContig}.fasta

done

rm sorted

