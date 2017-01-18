#!/usr/bin/perl
use warnings;
use strict;

my $finTE = $ARGV[0]; chomp $finTE;
my $fastain = $ARGV[1]; chomp $fastain;

my @TEfamilies = `cut -f8 $finTE |cut -f1 -d' '| sort | uniq`;chomp @TEfamilies;

foreach my $tefam (@TEfamilies)
{
my @reqContigs;
if ($tefam)
{
  my @CLids = `grep '$tefam     ' $finTE | cut -f1 `;chomp @CLids;
  foreach my $clid (@CLids)
  {
    my @contigids = `grep '$clid' $fastain |sed 's/>//'`;chomp @contigids;
    push @reqContigs, @contigids;
  }
#  $tefam =~ s/\.//g;
  open (outFAM, ">$tefam.txt") or die $!;
  foreach my $contig (@reqContigs)
  { print outFAM "$contig\n";}
  close outFAM;
}
}

