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
my $dxx;
my $dyy;
my $dzz;
my $dyx;
my $dzx;
my $dzy;
my $scale;
my $line;
my $nat;
my $iat;

#print "\n";
@args=@ARGV;
$errmess="usage: scale.pl input_filename scale_value\n";
@args>=2 || die $errmess;
#$nfile=scalar(@args);


#print $args[0],$args[1],"\n";
$fninp=$args[0];
$scale=$args[1];

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
        $dxx=@w[1]*$scale; $dyy=@w[2]*$scale; $dzz=@w[3]*$scale;
        $dyx=@w[4]*$scale; $dzx=@w[5]*$scale; $dzy=@w[6]*$scale;
    }
    else {
        @w=split(" ",$line);
        @w>=4 or die $errmess;
        @sat[$iat]=@w[0];
        @xat[$iat]=@w[1]*$scale;
        @yat[$iat]=@w[2]*$scale;
        @zat[$iat]=@w[3]*$scale;
        $iat++;
    }
}
close(INFILE);
if($iat!=$nat) {printf("ERROR: nat differes from number of atoms listed in the file: %6d %6d\n",$nat,$iat);}
printf("%s",$firstline);
printf("%s  %17.10f  %17.10f  %17.10f  %17.10f  %17.10f  %17.10f\n",$boundcond,$dxx,$dyy,$dzz,$dyx,$dzx,$dzy);
for($iat=0;$iat<$nat;$iat++) {
    printf("%5s  %24.15E  %24.15E  %24.15E  \n",@sat[$iat],@xat[$iat],@yat[$iat],@zat[$iat]);
}
