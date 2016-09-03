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
my ($dxx, $dyy, $dzz);
my $frac;
#my $x, $y, $z;
my @w;
my (@sat,@xat, @yat, @zat);
my (@satp,@xatp, @yatp, @zatp);
my $line;
my $boundcond;
my $complete;
my $nconf;
my ($iat, $nat, $natp);

#print "\n";
@args=@ARGV;
$errmess="usage: xyz-expand.pl input_filename frac\n";
@args>=2 || die $errmess;
#$nfile=scalar(@args);


$inputfilename=$args[0];
$frac=$args[1];

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
		$i=0;
		if(@w>=4) {
			$boundcond=@w[$i]; $i++;
		}
		else {
			die $errmess;
		}
		$dxx=@w[$i]; $i++;
		$dyy=@w[$i]; $i++;
		$dzz=@w[$i]; $i++;
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
        $natp=0;
        for($iat=0;$iat<$nat;$iat++) {
            my $l1=check_in_cell($dxx,@xat[$iat]);
            my $l2=check_in_cell($dyy,@yat[$iat]);
            my $l3=check_in_cell($dzz,@zat[$iat]);
            if(!$l1 || !$l2 || !$l3) {printf "ERROR: atom not in cell %5d%2d%2d%2d\n",$iat,$l1,$l2,$l3;die;}
            @satp[$natp]=@sat[$iat];
            @xatp[$natp]=@xat[$iat];
            @yatp[$natp]=@yat[$iat];
            @zatp[$natp]=@zat[$iat];
            $natp++;
        }
        for(my $iz=0;$iz<=1;$iz++) {
        for(my $iy=0;$iy<=1;$iy++) {
        for(my $ix=0;$ix<=1;$ix++) {
            if($ix==0 && $iy==0 && $iz==0) {next;}
            for($iat=0;$iat<$nat;$iat++) {
                if(@sat[$iat] eq 'H') {next;}
                my $x=@xat[$iat]+$ix*$dxx;
                my $y=@yat[$iat]+$iy*$dyy;
                my $z=@zat[$iat]+$iz*$dzz;
                if(!($x>(1.0+$frac)*$dxx) && !($y>(1.0+$frac)*$dyy) && !($z>(1.0+$frac)*$dzz)) {
                    @satp[$natp]=@sat[$iat];
                    @xatp[$natp]=$x;
                    @yatp[$natp]=$y;
                    @zatp[$natp]=$z;
                    $natp++;
                }
            }
        }
        }
        }
        #print the expanded system on screen
        printf "%8d  angstroemd0\n",$natp;
        printf "%s%24.15E%24.15E%24.15E\n",$boundcond,(1.0+$frac)*$dxx,(1.0+$frac)*$dyy,(1.0+$frac)*$dzz;
        for($iat=0;$iat<$natp;$iat++) {
            printf "%5s%24.15E%24.15E%24.15E\n",@satp[$iat],@xatp[$iat],@yatp[$iat],@zatp[$iat];
        }

	}
}
close(INFILE);


sub check_in_cell
{
    my $l=0;
    if(!($_[1]<0.0) && $_[1]<$_[0]) {$l=1}
    return $l;
}
