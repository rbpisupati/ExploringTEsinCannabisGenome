#!/usr/bin/perl
use warnings;
use strict;

my $in_file = $ARGV[0]; chomp $in_file;
my $avg_file= $ARGV[1]; chomp $avg_file;


open (in_pi,"$in_file") or die $!;
open (OUTFILE,">$ARGV[2]") ;

foreach my $each_line(<in_pi>)
{
  chomp $each_line;
  my @cover_1 = split("\t",$each_line);
  my $avg_cov = `grep -w '$cover_1[0]' $avg_file|grep -w '$cover_1[1]'|cut -f3`;chomp $avg_cov;

  print OUTFILE "$each_line\t$avg_cov\n";
}

