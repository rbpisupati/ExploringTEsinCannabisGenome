#!/usr/bin/perl
use warnings;
use strict;

my $in_file=$ARGV[0];chomp $in_file;
open (in_pileup,"$in_file") or die $!;
my $out_file = $ARGV[1];
open (OUTFILE,">$out_file") ;

foreach my $each_line (<in_pileup>)
{
  chomp $each_line;
  my @each_pos = split("\t",$each_line);
  my $depth = length($each_pos[4]);
  my $cou_a;my $cou_t;my $cou_g;my $cou_c;
 if ($each_pos[3] > 100)
 {
  my @ASCII = unpack ("C*", $each_pos[5]);
  my @passed_qual;
  for(my $i =0; $i < $depth;$i++)
  {
    my $val = $ASCII[$i];
    my $base = substr($each_pos[4],$i,1);
#    print "$base\n";
    if ($val -31 >= 20) 
    {
      my $qual = pack("C*", $val);push (@passed_qual, $qual);
#     print "$qual\n";
      if ($base =~ /[Aa]/) { $cou_a ++; }
      elsif ($base =~ /[Tt]/) { $cou_t ++; }
      elsif ($base =~ /[Gg]/) { $cou_g ++; }
      elsif ($base =~ /[Cc]/) { $cou_c ++; }
    }
  }
  my $per_a = int (100*($cou_a/$depth) +0.5); my $per_t = int (100*($cou_t/$depth) +0.5);
  my $per_g = int (100*($cou_g/$depth) +0.5); my $per_c = int (100*($cou_c/$depth) +0.5);
  my @total_str;
  for (1..$per_a) {push (@total_str,'A');}
  for (1..$per_t) {push (@total_str,'T');}
  for (1..$per_g) {push (@total_str,'G');}
  for (1..$per_c) {push (@total_str,'C');}
  
  my $seq = join '', @total_str;
  my $taken = length($seq);
#  my $taken_seq = substr ($seq,0,$taken);
#  my $req_seq = $passed_qual[0..$taken]
  my $qual_seq = join '', @passed_qual ;
  my $req_seq = substr ($qual_seq,0,$taken);
  print OUTFILE "$each_pos[0]\t$each_pos[1]\t$each_pos[2]\t$taken\t$seq\t$req_seq\n";
 }
 else {
 print OUTFILE "$each_line\n";
 }
}
