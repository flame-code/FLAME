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
#my $outputfilename;
#my $fn;
#my $nl;
my $errmess;
#my $firstline;
my @w;
#my @sat;
my @xat;
my @yat;
my @zat;
my $dnorm;
my $line;
#my $complete;
#my $nconf;
my $nat;
my $iat;
#my $xo;
#my $yo;
#my $zo;
#my $dx;
#my $dy;
#my $dz;
#my $ix;
#my $iy;
#my $iz;
#my $nx;
#my $ny;
#my $nz;
#my $nw;
#my $comment1;
#my $comment2;
#my $i;
#my @v;
#my $avgxyplanes;
#my $bohr=0.529177249;

#print "\n";
@args=@ARGV;
$errmess="usage: normalize.pl input_filename \n";
@args>=1 || die $errmess;
#$nfile=scalar(@args);


#print $args[0],$args[1],"\n";
$inputfilename=$args[0];

#$errmess="ERROR: output file exists. This script does not overwrite files \n";
#$errmess=$errmess."       in order to avoid possible mistake while ";
#$errmess=$errmess."providing multiple input files.\n";
#if(-e $outputfilename) {die $errmess};

$errmess="ERROR: cannot open file ".$inputfilename."\n";
open INFILE, "<$inputfilename" or die $errmess;


$nat=0;
#$ix=0;
#$iy=0;
#$iz=0;
$dnorm=0.0;
while($line=<INFILE>) {
	#$complete=0;
	#$nl++;
	@w=split(" ",$line);
	$xat[$nat]=$w[0];
	$yat[$nat]=$w[1];
	$zat[$nat]=$w[2];
	$dnorm+=$w[0]**2+$w[1]**2+$w[2]**2;
	$nat++;
}
close(INFILE);
$dnorm=sqrt($dnorm);
#print "dnorm ",$dnorm,"\n";
for($iat=0;$iat<$nat;$iat++) {
	$xat[$iat]=$xat[$iat]/$dnorm;
	$yat[$iat]=$yat[$iat]/$dnorm;
	$zat[$iat]=$zat[$iat]/$dnorm;
	printf "%25.15E  %25.15E  %25.15E  \n",$xat[$iat],$yat[$iat],$zat[$iat];
}


#if($ix!=$nx-1 || $ix!=$nx-1 || $ix!=$nx-1) {
#	print "ERROR: incomplete input file:\n";
#	print "nx,ny,nz ",$nx," ",$ny," ",$nz,"\n";
#	print "ix,iy,iz ",$ix," ",$iy," ",$iz,"\n";
#}
#for($iz=0;$iz<$nz;$iz++) {
#	$avgxyplanes=0.0;
#	for($iy=0;$iy<$ny;$iy++) {
#		for($ix=0;$ix<$nx;$ix++) {
#			$avgxyplanes+=$v[$ix,$iy,$iz];
#			#print $ix," ",$iy," ",$iz," ",$v[$ix,$iy,$iz],"\n";
#		}
#	}
#	printf "%20.10f  %20.5E  \n",$zo+$iz*$dz,$avgxyplanes/($nx*$ny);
#}
