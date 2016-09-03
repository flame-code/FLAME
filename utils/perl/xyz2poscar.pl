#!/usr/bin/env perl
#;-*- Perl -*-

use strict;
#use FindBin qw($Bin);
#use lib "$Bin";
#use Vasp;
#use Math::Trig;
#$fact=180/pi;
my @args;
my $fninp;
my $nl;
my $errmess;
my $firstline;
my $boundcond;
my @w;
my @sat;
my @xat;
my @yat;
my @zat;
my @xms;
my @yms;
my @zms;
my $dxx;
my $dyy;
my $dzz;
my $dyx;
my $dzx;
my $dzy;
my $x;
my $y;
my $z;
my $line;
my $nat;
my $iat;
my $i;
my @natarr;
my @typat;
my $ntypat;
my $new;

#print "\n";
@args=@ARGV;
$errmess="usage: xyz2poscar.pl input_filename\n";
@args>=1 || die $errmess;
#$nfile=scalar(@args);


#print $args[0],$args[1],"\n";
$fninp=$args[0];

#print $args[0];
$errmess="ERROR: cannot open file ".$fninp."\n";
open INFILE, "<$fninp" or die $errmess;

$errmess="ERROR: something wrong with file contents.\n";
$nl=0;
while($line=<INFILE>) {
	$nl++;
	if($nl==1) {
        $iat=0;
        $firstline=$line;
        @w=split(" ",$line);
        @w>=1 or die $errmess;
        $nat=@w[0];
    }
    elsif($nl==2) {
        @w=split(" ",$line);
        @w>=7 or die $errmess;
        $boundcond=@w[0];
        $dxx=@w[1]; $dyy=@w[2]; $dzz=@w[3];
        $dyx=@w[4]; $dzx=@w[5]; $dzy=@w[6];
    }
    else {
        @w=split(" ",$line);
        @w>=5 or die $errmess;
        @sat[$iat]=@w[0];
        @xat[$iat]=@w[1];
        @yat[$iat]=@w[2];
        @zat[$iat]=@w[3];
		@xms[$iat]=substr(@w[4],0,1); #split("",@w[4]);
		@yms[$iat]=substr(@w[4],1,1);
		@zms[$iat]=substr(@w[4],2,1);
        $iat++;
    }
}
close(INFILE);
if($iat!=$nat) {printf("ERROR: nat differes from number of atoms listed in the file: %6d %6d\n",$nat,$iat);}

$ntypat=0;
for($iat=0;$iat<$nat;$iat++) {
	$new=1;
	for($i=0;$i<$ntypat;$i++) {
		if(@sat[$iat] eq @typat[$i]) {$new=0;last;}
	}
	#printf("%4d  %4d  %2s  %2s\n",$iat,$new,@sat[$iat],@typat[$i]);
	if($new==1) {
		@typat[$ntypat]=@sat[$iat];
		#printf("%5s\n",@typat[$ntypat]);
		$ntypat++;
		@natarr[$ntypat-1]=1;
	}
	else {
		@natarr[$ntypat-1]++;
	}
}
for($i=0;$i<$ntypat;$i++) {printf("%5s",@typat[$i]);}
printf("\n");
printf("  1.0  \n");
printf("  %17.10f  %17.10f  %17.10f\n",$dxx,0.00,0.00);
printf("  %17.10f  %17.10f  %17.10f\n",$dyx,$dyy,0.00);
printf("  %17.10f  %17.10f  %17.10f\n",$dzx,$dzy,$dzz);
for($i=0;$i<$ntypat;$i++) {printf("%5s",@typat[$i]);}
printf("\n");
for($i=0;$i<$ntypat;$i++) {printf("%5d",@natarr[$i]);}
printf("\n");
printf("selective dynamics\n");
printf("cartesian\n");
#printf("Direct\n");
for($iat=0;$iat<$nat;$iat++) {
	$x=@xat[$iat];
	$y=@yat[$iat];
	$z=@zat[$iat];
    #if($x<0.0) {$x+=$dxx};
	#if($y<0.0) {$y+=$dyy};
	#if($z<0.0) {$z+=$dzz};
	#if($x>$dxx) {$x-=$dxx};
	#if($y>$dyy) {$y-=$dyy};
	#if($z>$dzz) {$z-=$dzz};
	printf("  %24.15E  %24.15E  %24.15E  %1s  %1s  %1s  \n",$x,$y,$z,@xms[$iat],@yms[$iat],@zms[$iat]);
	#$x=@xat[$iat]/$dxx;
	#$y=@yat[$iat]/$dyy;
	#$z=@zat[$iat]/$dzz;
    #printf("  %24.15f  %24.15f  %24.15f  %1s  %1s  %1s  \n",$x,$y,$z,@xms[$iat],@yms[$iat],@zms[$iat]);
}
