#!/usr/bin/awk -f
function isnum(x){return(x==x+0)}
BEGIN {
	if(ARGC<2) {print "ERROR: argument(s) are missing";exit;}
	if(ARGC>2) {print "WARNING: too many arguments";}
    if(cvxx==0.0 || cvyy==0.0 || cvzz==0.0) {
        print "usage: xyz2xyz.awk -v cvxx=value -v cvyy=value -v cvzz=value filename";
        exit;
    }
}
{
    nline++;
    if(nline==1) {print $0;}
    else if(nline==2) {printf("free%22.16f%22.16f%22.16f  0.0  0.0  0.0\n",cvxx,cvyy,cvzz);}
    else {printf("%s  TTT \n",$0);}
}
