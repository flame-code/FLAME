#*****************************************************************************************
def get_atomic_symbol(iatom):
    elemens=['H','He','Li','Be','B','C','N','O','F','Ne','Na',
          'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti',
          'V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
          'Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru',
          'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs',
          'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
          'Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir',
          'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
          'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es',
          'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt']
    if iatom<1 or iatom>109:
        print "ERROR: unknown atomic number in get_atomic_symbol %5d" % iatom
        exit()
    else:
        return elemens[iatom-1]
#*****************************************************************************************
def get_atomic_number(sat):
    if sat=='H':
        iatom=1
    elif sat=='He':
        iatom=2
    elif sat=='C':
        iatom=6
    elif sat=='N':
        iatom=7
    elif sat=='O':
        iatom=8
    elif sat=='Na':
        iatom=11
    elif sat=='Cl':
        iatom=17
    elif sat=='Cu':
        iatom=29
    else:
        print "ERROR: unknown atomic symbol in get_atomic_number %5s" % sat
    return iatom
#*****************************************************************************************
