#!/usr/bin/env perl

use FindBin qw($Bin);

$filename="o1";
$nargs = $#ARGV + 1;
if($nargs>0) {
    if($ARGV[0] eq "-last") {$last_energ=1;}
    else {$filename=$ARGV[0];}
    if($nargs==2) {
        $filename=$ARGV[0];
        if($ARGV[1] eq "-last") {$last_energ=1;}
    }
}
#printf "%d  %d  %s\n\n",$last_energ,$nargs,$ARGV[0];

unless(-f $filename) {
    printf "\nERROR: file %s does not exit.\n",$filename;
    if($nargs<1) {printf "Notice: no arguement was provided so default filename o1 is considered.\n";}
    printf "\n";
    exit 0;
}

$forces = `grep -a 'Maximum force component' $filename` ;
$energy = `grep -a 'Total energy uncorrected' $filename` ;

@forces = split /\n/ , $forces ;
@energy = split /\n/ , $energy ; 

$num = @energy;
#open OUT , ">fe.dat" ;

if($last_energ) {$iprint=$num-2;}
else {$iprint=-1;}
for($i=0; $i<$num; $i++){
    $line = $forces[$i] ;
    chomp($line) ;
    $line=~s/^\s+//g;
    @line=split /\s+/,$line;
    $f = $line[4],"\n" ;

    $line = $energy[$i] ; 
    chomp($line) ;
    $line=~s/^\s+//g;
    @line=split /\s+/,$line;
    $e = $line[5],"\n" ;
    if(!$i) {$e0 = $e ;$eold=$e;}
    #printf OUT "%5i %20.8f %20.6f %20.6e %20.6e \n",$i,$f,$e,$e-$eold,$e-$e0;
    if($i>$iprint) {
        printf "%5i %20.8f %20.6f %20.6e %20.6e \n",$i,$f,$e,$e-$eold,$e-$e0;
    }
    $eold=$e;
}

#close OUT ;

