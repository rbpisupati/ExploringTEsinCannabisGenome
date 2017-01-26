#!/usr/bin/perl
#Takes a fasta file and produces a filtered fasta file containing only the sequences longer than a minimum length, 100
$FASTA = shift @ARGV;
@fasta = ();
$count = 0;
$count2 = 0;
$max = 0;
$min = shift @ARGV;
$lengthtot=0;
@lengths = ();

open FASTA or die "No file $FASTA";
open OUTFILE, ">$FASTA.longest.min$min.fa";
{
        local $/ = '>';
        while (<FASTA>) {


                @temp = split /\n/, $_;
                $name = shift @temp;
                $seq2 = join '', @temp;
                $length = length $seq2;
                $lengthtot= $lengthtot+$length;
                push @lengths, $length;
                if ($length > $max) {
                        $max = $length;
                        print "$max\n";
                }
                if ($length > $min){
                        $count ++;
                $seq2 =~ s/>//g;
                print OUTFILE '>', "$name\n$seq2\n";
        }
                else {
                        $count2 ++;
                }

        }
}
close FASTA;

close OUTFILE;
$seqnum = $count + $count2;
$average = $lengthtot / ($seqnum);

print "There are $count sequences at least $min nucleotides in length and $count2 less than $min nucleotides in length. Total length is $lengthtot .
The longest sequence is $max nucleotides in length.\nAverage sequence length of all $seqnum sequences is $average.\n";


$n50 = 0;
$n90 = 0;
@lengths = sort { $b <=> $a }  @lengths;
$total = $count + $count2;
$avg = $lengthtot / $total;
$n50 = 0;
$y = 0;
for ($x=0; $x < ($lengthtot / 2); $x = ($x + $lengths[$y])){
        $n50 = $lengths[$y];
        $y=$y++;
}

for ($x=0; $x < ($lengthtot * 0.9); $x = ($x + $lengths[$y])){
        $n90 = $lengths[$y];
        $y=$y++;
}


print "\nAverage length $avg\nTotal sequence length $totseq\n\n";

@hist=();
for ($i = 0; $i < $max; $i++){
        $hist[$i]=0;
}
foreach (@lengths){
        $hist[$_]++;
}
for ($i = 0; $i < $max; $i++){
        print "$i\t$hist[$i]\n";
}

