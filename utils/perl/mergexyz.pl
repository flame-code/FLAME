#!/usr/bin/env perl
#;-*- Perl -*-

use strict;
#use FindBin qw($Bin);
#use lib "$Bin";
#use Vasp;
#use Math::Trig;
#$fact=180/pi;
my @args;
my $inputfilename;
my $outputfilename;
my $nl;
my $nw;
my $errmess;
#my $errmess1;
#my $errmess2;
my @w;
my $line;
my $ifile;
my $nfile;
my $nat;
my $iat;

#print "\n";
@args=@ARGV;
$errmess="usage: mergexyz.pl input_filename output_filename \n";
$errmess=$errmess."       input can be multiple files as standard Linux.\n";
@args>=2 || die $errmess;
$nfile=scalar(@args);


#print $args[0],$args[1],"\n";
$outputfilename=$args[$nfile-1];

$errmess="ERROR: output file exists. This script does not overwrite files \n";
$errmess=$errmess."       in order to avoid possible mistake while ";
$errmess=$errmess."providing multiple input files.\n";
#$errmess=$errmess1.$errmess2;
if(-e $outputfilename) {die $errmess};
$errmess="ERROR: cannot open file ".$outputfilename."\n";
open OUTFILE, ">$outputfilename" or die $errmess;
for($ifile=0;$ifile<$nfile-1;$ifile++) {
	$inputfilename=$args[$ifile];
	print "Reading file ",$inputfilename," and writing to ",$outputfilename,"\n";
	$errmess="ERROR: cannot open file ".$inputfilename."\n";
	open INFILE, "<$inputfilename" or die $errmess;
	$nl=0;
	$nw=0;
	$iat=0;
	$nat=0;
	$errmess="ERROR: first of the first line of xyz file must be # of atoms.";
	#while($line=<INFILE> && $iat<=$nat) {
	while($line=<INFILE>) {
		$nl++;
		if($nl==1) {
			@w=split(" ",$line);
			@w>=1 or die $errmess;
			$nat=@w[0];
			chop $line;
			$line=$line."  ".$inputfilename."  \n";
		}
		if($nl>2) {$iat++};
		#$nw=scalar(@w);
		print OUTFILE $line;
		#print $nl,"  ",$nat,"\n";
		#print $nl,"  ",$nw,"\n";
		#print $inputfilename,"  ",$nat,"\n";
		last if($iat==$nat);
	}
	$iat==$nat or die "ERROR: number of atoms not as indicated in the first line.\n";
	close(INFILE);
}
close(OUTFILE);
