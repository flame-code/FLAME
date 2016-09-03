#*****************************************************************************************
def get_atomic_symbol(iatom):
    if iatom==1:
        sat='H'
    elif iatom==2:
        sat='He'
    elif iatom==6:
        sat='C'
    elif iatom==7:
        sat='N'
    elif iatom==8:
        sat='O'
    elif iatom==11:
        sat='Na'
    elif iatom==17:
        sat='Cl'
    elif iatom==29:
        sat='Cu'
    else:
        print "ERROR: unknown atomic number in get_atomic_symbol %5d" % iatom
    return sat
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
