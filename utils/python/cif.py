import sys
import atoms
import math as mt
import numpy as np
from atoms import *
#*****************************************************************************************
def cif_write(atoms,filename):
    a=np.sqrt(atoms.cellvec[0][0]**2+atoms.cellvec[0][1]**2+atoms.cellvec[0][2]**2)
    b=np.sqrt(atoms.cellvec[1][0]**2+atoms.cellvec[1][1]**2+atoms.cellvec[1][2]**2)
    c=np.sqrt(atoms.cellvec[2][0]**2+atoms.cellvec[2][1]**2+atoms.cellvec[2][2]**2)
    alpha=mt.acos(np.dot(np.array(atoms.cellvec[1][:]),np.array(atoms.cellvec[2][:]))/(b*c))*180./mt.pi
    beta =mt.acos(np.dot(np.array(atoms.cellvec[2][:]),np.array(atoms.cellvec[0][:]))/(c*a))*180./mt.pi
    gamma=mt.acos(np.dot(np.array(atoms.cellvec[0][:]),np.array(atoms.cellvec[1][:]))/(a*b))*180./mt.pi
    if filename=="screen":
        print "%s %10.15f" % ("_cell_length_a",a)
        print "%s %10.15f" % ("_cell_length_b",b)
        print "%s %10.15f" % ("_cell_length_c",c)
        print ""
        print "%s %10.15f" % ("_cell_angle_alpha", alpha)
        print "%s %10.15f" % ("_cell_angle_beta ", beta )
        print "%s %10.15f" % ("_cell_angle_gamma", gamma)
        print ""
        print "_symmetry_space_group_name_H-M    'P 1'"
        print "_symmetry_int_tables_number        1"
        print ""
        print "loop_"
        print "  _symmetry_equiv_pos_as_xyz"
        print "  'x, y, z'"
        print ""
        print "loop_"
        print "  _atom_site_label"
        print "  _atom_site_type_symbol"
        print "  _atom_site_occupancy"
        print "  _atom_site_fract_x"
        print "  _atom_site_fract_y"
        print "  _atom_site_fract_z"
        latinv=np.matrix('1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0')
        if(atoms.coordinates=="Cartesian"):
            latmat=np.asmatrix(atoms.cellvec)
            latinv=np.linalg.inv(latmat)
        elno = dict()
        for i in range(atoms.nat):
            if atoms.sat[i] in elno:
                elno[atoms.sat[i]]+=1
            else:
                elno[atoms.sat[i]]=1
            ind=elno[atoms.sat[i]]
            xred=np.dot(latinv.T,np.array(atoms.rat[i][:]))
            x=xred.item(0)
            y=xred.item(1)
            z=xred.item(2)
            elind = str(atoms.sat[i])+str(ind)
            print "%5s %s %0.2f %10.15f %10.15f %10.15f" % (elind,atoms.sat[i],1.,x,y,z)
    else:
        f= open(filename,"w")
        f.write( "%s %10.15f\n" % ("_cell_length_a",a))
        f.write( "%s %10.15f\n" % ("_cell_length_b",b))
        f.write( "%s %10.15f\n" % ("_cell_length_c",c))
        f.write( "\n")
        f.write( "%s %10.15f\n" % ("_cell_angle_alpha", alpha))
        f.write( "%s %10.15f\n" % ("_cell_angle_beta ", beta ))
        f.write( "%s %10.15f\n" % ("_cell_angle_gamma", gamma))
        f.write( "\n")
        f.write( "_symmetry_space_group_name_H-M    'P 1'\n")
        f.write( "_symmetry_int_tables_number          1 \n")
        f.write( "\n")
        f.write( "loop_\n")
        f.write( "  _symmetry_equiv_pos_as_xyz\n")
        f.write( "  'x, y, z'\n")
        f.write( "\n")
        f.write( "loop_\n")
        f.write( "  _atom_site_label\n")
        f.write( "  _atom_site_type_symbol\n")
        f.write( "  _atom_site_occupancy\n")
        f.write( "  _atom_site_fract_x\n")
        f.write( "  _atom_site_fract_y\n")
        f.write( "  _atom_site_fract_z\n")
        latinv=np.matrix('1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0')
        if(atoms.coordinates=="Cartesian"):
            latmat=np.asmatrix(atoms.cellvec)
            latinv=np.linalg.inv(latmat)
        elno = dict()
        for i in range(atoms.nat):
            if atoms.sat[i] in elno:
                elno[atoms.sat[i]]+=1
            else:
                elno[atoms.sat[i]]=1
            ind=elno[atoms.sat[i]]
            xred=np.dot(latinv.T,np.array(atoms.rat[i][:]))
            x=xred.item(0)
            y=xred.item(1)
            z=xred.item(2)
            elind = str(atoms.sat[i])+str(ind)
            f.write( "%5s %s %0.2f %10.15f %10.15f %10.15f\n" % (elind,atoms.sat[i],1.,x,y,z))
        f.close()
