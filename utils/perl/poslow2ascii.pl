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
my $fnbase;
my $fn;
my $nl;
my $i;
#my $nw;
my $errmess;
my $firstline;
#my $errmess2;
my $dxx;
my $dyy;
my $dzz;
my $dyx;
my $dzx;
my $dzy;
my @w;
my @sat;
my @xat;
my @yat;
my @zat;
my $line;
my $boundcond;
my $complete;
my $nconf;
my $nat;
my $iat;
my $trans;

#print "\n";
@args=@ARGV;
$errmess="usage: poslow2ascii.pl input_filename \n";
@args>=1 || die $errmess;
#$nfile=scalar(@args);


$inputfilename=$args[0];

$errmess="ERROR: cannot open file ".$inputfilename."\n";
open INFILE, "<$inputfilename" or die $errmess;

$fnbase = substr($inputfilename,0,-4);

$nl=0;
$nconf=0;
while($line=<INFILE>) {
	$complete=0;
	$nl++;
	if($nl==1) {
		$firstline=$line;
		@w=split(" ",$line);
		$errmess="ERROR: first of the first line of xyz file must be # of atoms.";
		@w>=1 or die $errmess;
		$nat=@w[0];
		#chop $line;
		#$line=$line."  ".$inputfilename."  \n";
		$iat=0;
	}
	elsif($nl==2) {
		@w=split(" ",$line);
		$errmess="ERROR: second line cannot be parsed.";
		#@w==7 or die $errmess;
		$i=0;
		if(@w==7) {
			$boundcond=@w[$i]; $i++;
		}
		elsif(@w==6) {
		}
		else {
			die $errmess;
		}
		$dxx=@w[$i]; $i++;
		$dyy=@w[$i]; $i++;
		$dzz=@w[$i]; $i++;
		$dyx=@w[$i]; $i++;
		$dzx=@w[$i]; $i++;
		$dzy=@w[$i]; $i++;
	}
	else {
		@w=split(" ",$line);
		@w>=4 or die $errmess;
		@sat[$iat]=@w[0];
		@xat[$iat]=@w[1];
		@yat[$iat]=@w[2];
		@zat[$iat]=@w[3];
		$iat++;
		if($iat==$nat) {$complete=1;}
	}
	if($complete==1) {
		$nl=0;
		$nconf++;
		$fn=sprintf "%5.5d",$nconf;
		$outputfilename=$fnbase.$fn.".ascii";
		$errmess="ERROR: cannot open file ".$outputfilename."\n";
		open OUTFILE, ">$outputfilename" or die $errmess;
		print OUTFILE $firstline;
		printf OUTFILE  "%24.15E%24.15E%24.15E\n",$dxx,$dyx,$dyy;
		printf OUTFILE  "%24.15E%24.15E%24.15E\n",$dzx,$dzy,$dzz;
		for($iat=0;$iat<$nat;$iat++) {
			$trans=0.0;
			printf OUTFILE "%24.15E%24.15E%24.15E%5s\n",@xat[$iat]+$trans,@yat[$iat]+$trans,@zat[$iat]+$trans,@sat[$iat];
		}
		close(OUTFILE);
	}
}
close(INFILE);
