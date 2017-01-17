#!/usr/bin/perl
use warnings;
use strict;

my $divFile=$ARGV[0];chomp $divFile;
my $summFile=$ARGV[1]; chomp $summFile;

open (diverCL, "$divFile") or die $!;
open (outFile, ">$divFile.TEfam.txt");

foreach my $eachLine (<diverCL>)
{
  chomp $eachLine;
  my @eachCL = split("\t",$eachLine);
  my $CLno = $eachCL[0];
  
  my $repFam = `grep '$CLno	' $summFile |cut -f7`; chomp $repFam;
  my $clPer  = `grep '$CLno	' $summFile |cut -f5`; chomp $clPer;
  
  print outFile "$eachLine\t$repFam\t$clPer\n";

}
