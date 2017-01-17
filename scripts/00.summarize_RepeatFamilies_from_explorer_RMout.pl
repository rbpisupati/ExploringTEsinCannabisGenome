#!/usr/bin/perl
#use warnings;
use strict;

my $in_rm = $ARGV[0]; chomp $in_rm;
my $repcon = $ARGV[1];chomp $repcon;
my $out_file= $ARGV[2]; chomp $out_file;
open (rm_output,"$in_rm") or die $!;
open (OUTFILE, ">$out_file");

my @all_calc;my @tot_base;my $total_reads;
my @all_bp = `grep '^content' $in_rm|cut -f2`;
my $total;
foreach(@all_bp)
{
	$total = $total + eval $_;
}
my $header = `head -1 $in_rm`;chomp $header;
my @header_tab = split ("\t",$header);
foreach (<rm_output>)
{
	if (/^content/)
	{
		my @each_con = split("\t",$_);
		$total_reads += eval $each_con[2];
		my $top=0;my $max=0;
		for(my $i=3;$i<@each_con;$i++)
		{
			if($each_con[$i]>$max){$max= eval $each_con[$i];$top= $i;}
		}
		my $add = (eval $each_con[1]/$total)*$repcon;
		$all_calc[$top] = $add + eval $all_calc[$top];
	}
}
my $classfied;
for(my $i=3; $i< @all_calc; $i++) 
{
	if($all_calc[$i]) { $classfied = $classfied + $all_calc[$i]; }
}
my $unknown = $repcon - $classfied ;
my $flag=0;my $tot_unkn;
for(my $i=3;$i<@all_calc;$i++)
{
	if($header_tab[$i] =~ /Unknown/i)
	{	
		$tot_unkn = $all_calc[$i] + $unknown;
		print OUTFILE "Unknown\t$tot_unkn\n";
		$flag =1;
	}
	elsif ($all_calc[$i])
	{
		print OUTFILE "$header_tab[$i]\t$all_calc[$i]\n";
	}
}
if($flag == 0 and $unknown >0 )
{
	print OUTFILE "Unknown\t$unknown\n";
}
print OUTFILE "Total content\t$repcon\n";
