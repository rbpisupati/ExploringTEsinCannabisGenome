#!/usr/bin/perl
use warnings;
use strict;

my $in_file = $ARGV[0];
chomp $in_file;
my $out_file = $ARGV[1];
open (pileup_in, "$in_file") or die $!;
open (OUTFILE,">$out_file") ;

foreach my $each_pos (<pileup_in>)
{
  chomp $each_pos;
  my @pile_line = split("\t",$each_pos) ;
  my @bases =();my @quality =();
  # convert Sanger-style phred quality scores to more typical illumina quality scores
  my @ASCII = unpack ("C*", $pile_line[5]);
  foreach my $val (@ASCII) {
        $val = $val +31;
  }
  my $mod_qual = pack ("C*", @ASCII);
  my $g=0;
  for (my $i=0;$i<length($pile_line[4]);$i++){
    my $base = substr($pile_line[4], $i, 1);
    if ($base =~ /[atgcnNATGC]/) {
        push @bases, $base;
        my $qual = substr($mod_qual, $g, 1);
        $g ++;
        push @quality, $qual;
    }
    elsif($base eq '$' or $base eq '*' or $base eq '<' or $base eq '>'){
        next;
    }
    elsif($base eq '^'){
        $i = $i+1;
        next;
    }
    elsif($base eq '-' or $base eq '+'){
        my @nums = ();
        my $num = 0;
        while((substr($pile_line[4], $i+1, 1)) =~ /[0-9]/){
            push @nums, (substr($pile_line[4], $i+1, 1));
            $i = $i+1;
        }
        $num = join '', @nums;
        $i = $i+$num    ;
#        print "number is $num\n";
        next;
    }
  
    else { print "not accounted for: base equals $base\n";}
  }  
  my $seq = join '', @bases;
  my $qual_seq = join '',@quality;
  my $cov = $g ;
  print OUTFILE "$pile_line[0]\t$pile_line[1]\t$pile_line[2]\t$cov\t$seq\t$qual_seq\n";
#    my $seq_sub = substr($seq,0,50);
#    my $qual_sub = substr ($qual_seq,0,50);
#    my $cov = length($seq_sub);
#    print OUTFILE "$pile_line[0]\t$pile_line[1]\t$pile_line[2]\t$cov\t$seq_sub\t$qual_sub\n";
}
