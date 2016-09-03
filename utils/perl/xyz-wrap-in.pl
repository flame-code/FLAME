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
my @bemoved;
my $line;
my $boundcond;
my $complete;
my $nconf;
my $nat;
my $iat;

#print "\n";
@args=@ARGV;
$errmess="usage: xyz-wrap-in.pl input_filename \n";
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
		@w>=5 or die $errmess;
		@sat[$iat]=@w[0];
		@xat[$iat]=@w[1];
		@yat[$iat]=@w[2];
		@zat[$iat]=@w[3];
		@bemoved[$iat]=@w[4];
		$iat++;
		if($iat==$nat) {$complete=1;}
	}
	if($complete==1) {
		$nl=0;
		$nconf++;
		$fn=sprintf "%5.5d",$nconf;
		$outputfilename=$fnbase.$fn.".ascii";
		$errmess="ERROR: cannot open file ".$outputfilename."\n";
        #open OUTFILE, ">$outputfilename" or die $errmess;
        #print OUTFILE $firstline;
        #print $firstline;
        printf "%8d  angstroemd0\n",$nat;
        #printf OUTFILE  "%24.15E%24.15E%24.15E\n",$dxx,$dyx,$dyy;
        #printf OUTFILE  "%24.15E%24.15E%24.15E\n",$dzx,$dzy,$dzz;
        printf "%s%24.15E%24.15E%24.15E\n",$boundcond,$dxx,$dyy,$dzz;
        for($iat=0;$iat<$nat;$iat++) {
            @xat[$iat]=wrap_in($dxx,@xat[$iat]);
            @yat[$iat]=wrap_in($dyy,@yat[$iat]);
            @zat[$iat]=wrap_in($dzz,@zat[$iat]);
            printf "%5s%24.15E%24.15E%24.15E%5s\n",@sat[$iat],@xat[$iat],@yat[$iat],@zat[$iat],@bemoved[$iat];
        }
        #close(OUTFILE);
        #printf "%5.1f\n",wrap_in(10.0,0.0);
	}
}
close(INFILE);


sub wrap_in
{
    my $res;
    $res=$_[1];
    while($res<0.0) {
        $res+=$_[0];
        #printf "A: %5.1f%5.1f\n",$res,$_[0];
    }
    while(!($res<$_[0])) {
        $res-=$_[0];
        #printf "B: %5.1f%5.1f\n",$res,$_[0];
    }
    return $res;
}
