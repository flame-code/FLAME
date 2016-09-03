#!/usr/bin/awk -f
function isnum(x){return(x==x+0)}
BEGIN {
	if(ARGC<2) {print "ERROR: argument(s) are missing";exit;}
	if(ARGC>2) {print "WARNING: too many arguments";}
    if(scale==0.0) {print "usage: scaleposcar.awk -v scale=some_value filename";exit;}
    #scale=ARGV[2];
    #FILENAME=ARGV[1];
	#printf("%10.5f\n",scale);
    xyz_lines=0;
    selective="none";
}
{
    #printf("%s\n",$0);
    nline++;
    if(nline==1) {print $0;}
    else if(nline==2) {vscale=$1;}
    else if(nline==3) {
        printf("%22.16f%22.16f%22.16f\n",scale*$1,scale*$2,scale*$3);
    }
    else if(nline==4) {
        printf("%22.16f%22.16f%22.16f\n",scale*$1,scale*$2,scale*$3);
    }
    else if(nline==5) {
        printf("%22.16f%22.16f%22.16f\n",scale*$1,scale*$2,scale*$3);
    }
    else if(nline==6) {
        if(isnum($1)) {version=4;}
        else {version=5;}
        print $0;
        #print isnum("hello"),isnum("-42"),isnum($1)
    }
    else if(nline>6 && nline<10 && xyz_lines==0) {
        print $0
        ch=substr($1,1,1);
        if(ch=="s" || ch=="S") {selective="Selective dynamics";}
        if(ch=="d" || ch=="D") {coordinates="Direct";xyz_lines=1;next;}
        else if(ch=="c" || ch=="C") {coordinates="Cartesian";xyz_lines=1;next;}
        #print "xyz_lines ",xyz_lines;
    }
    if(xyz_lines==1) {
        if(coordinates=="Cartesian") {
            printf("%22.16f%22.16f%22.16f",scale*$1,scale*$2,scale*$3);
            if(selective=="Selective dynamics") {printf("%2s%2s%2s\n",$4,$5,$6);}
            else {printf("\n");}
        }
        else {
            print $0
        }
    }
}
