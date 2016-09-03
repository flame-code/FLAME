#!/usr/bin/python
#print ( 'set_rcov module imported' )
rcov=0.0
def setrcov (sat):
    if sat=='LJ':
        rcov=2.0**(1.0/6.0)/2.0
    elif sat=='H' : 
        rcov=0.31
    elif sat=='He':
        rcov=0.28
    elif sat=='Li':
        rcov=1.28
    elif sat=='Be': 
        rcov=0.96
    elif sat=='B' :
        rcov=0.84
    elif sat=='C' :
        rcov=0.76
    elif sat=='N' :
        rcov=0.71
    elif sat=='O' :
        rcov=0.66
    elif sat=='F' :
        rcov=0.57
    elif sat=='Ne':
        rcov=0.58
    elif sat=='Na':
        rcov=1.66
    elif sat=='Mg':
        rcov=1.41
    elif sat=='Al':
        rcov=1.21
    elif sat=='Si':
        rcov=1.11
    elif sat=='P' :
        rcov=1.07
    elif sat=='S' :
        rcov=1.05
    elif sat=='Cl':
        rcov=1.02
    elif sat=='Ar':
        rcov=1.06
    elif sat=='K' :
        rcov=2.03
    elif sat=='Ca':
        rcov=1.76
    elif sat=='Sc':
        rcov=1.70
    elif sat=='Ti':
        rcov=1.60
    elif sat=='V' :
        rcov=1.53
    elif sat=='Cr':
        rcov=1.39
    elif sat=='Mn':
        rcov=1.39
    elif sat=='Fe':
        rcov=1.32
    elif sat=='Co':
        rcov=1.26
    elif sat=='Ni':
        rcov=1.24
    elif sat=='Cu':
        rcov=1.32
    elif sat=='Zn':
        rcov=1.22
    elif sat=='Ga':
        rcov=1.22
    elif sat=='Ge':
        rcov=1.20
    elif sat=='As':
        rcov=1.19
    elif sat=='Se':
        rcov=1.20
    elif sat=='Br':
        rcov=1.20
    elif sat=='Kr':
        rcov=1.16
    elif sat=='Rb':
        rcov=2.20
    elif sat=='Sr':
        rcov=1.95
    elif sat=='Y' :
        rcov=1.90
    elif sat=='Zr':
        rcov=1.75
    elif sat=='Nb':
        rcov=1.64
    elif sat=='Mo':
        rcov=1.54
    elif sat=='Tc':
        rcov=1.47
    elif sat=='Ru':
        rcov=1.46
    elif sat=='Rh':
        rcov=1.42
    elif sat=='Pd':
        rcov=1.39
    elif sat=='Ag':
        rcov=1.45
    elif sat=='Cd':
        rcov=1.44
    elif sat=='In':
        rcov=1.42
    elif sat=='Sn':
        rcov=1.39
    elif sat=='Sb':
        rcov=1.39
    elif sat=='Te':
        rcov=1.38
    elif sat=='I' :
        rcov=1.39
    elif sat=='Xe':
        rcov=1.40
    elif sat=='Cs':
        rcov=2.44
    elif sat=='Ba':
        rcov=2.15
    elif sat=='La':
        rcov=2.07
    elif sat=='Ce':
        rcov=2.04
    elif sat=='Pr':
        rcov=2.03
    elif sat=='Nd':
        rcov=2.01
    elif sat=='Pm':
        rcov=1.99
    elif sat=='Sm':
        rcov=1.98
    elif sat=='Eu':
        rcov=1.98
    elif sat=='Gd':
        rcov=1.96
    elif sat=='Tb':
        rcov=1.94
    elif sat=='Dy':
        rcov=1.92
    elif sat=='Ho':
        rcov=1.92
    elif sat=='Er':
        rcov=1.89
    elif sat=='Tm':
        rcov=1.90
    elif sat=='Yb':
        rcov=1.87
    elif sat=='Lu':
        rcov=1.87
    elif sat=='Hf':
        rcov=1.75
    elif sat=='Ta':
        rcov=1.70
    elif sat=='W' :
        rcov=1.62
    elif sat=='Re':
        rcov=1.51
    elif sat=='Os':
        rcov=1.44
    elif sat=='Ir':
        rcov=1.41
    elif sat=='Pt':
        rcov=1.36
    elif sat=='Au':
        rcov=1.36
    elif sat=='Hg':
        rcov=1.32
    elif sat=='Tl':
        rcov=1.45
    elif sat=='Pb':
        rcov=1.46
    elif sat=='Bi':
        rcov=1.48
    elif sat=='Po':
        rcov=1.40
    elif sat=='At':
        rcov=1.50
    elif sat=='Rn':
        rcov=1.50
    elif sat=='Fr':
        rcov=2.60
    elif sat=='Ra':
        rcov=2.21
    elif sat=='Ac':
        rcov=2.15
    elif sat=='Th':
        rcov=2.06
    elif sat=='Pa':
        rcov=2.00
    elif sat=='U' :
        rcov=1.96
    elif sat=='Np':
        rcov=1.90
    elif sat=='Pu':
        rcov=1.87
    elif sat=='Am':
        rcov=1.80
    elif sat=='Cm':
        rcov=1.69
    else :
        print 'ERROR: no covalent radius stored for atomtype=',sat
    return rcov

if __name__ == '__main__':
    print ( 'test: radius of covalent bond of "H" is: ', setrcov('H')  )
else: 
    print('name is not main')
