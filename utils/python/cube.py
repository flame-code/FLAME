from atoms import *
from atomic_data import *
from poisson import *
#*****************************************************************************************
def cube_advance_index(poisson,igpx,igpy,igpz):
    igpz+=1
    if igpz==poisson.ngpz:
        igpz=0
        igpy+=1
        if igpy==poisson.ngpy:
            igpy=0
            igpx+=1
    return igpx,igpy,igpz
#*****************************************************************************************
def cube_read(filename):
    f=open(filename,"r")
    atoms=Atoms()
    poisson=Poisson()
    atoms.boundcond="bulk"
    iline=0
    atoms.nat=0
    iat=0
    for line in f.readlines():
        iline+=1
        if iline==1:
            comment1=line
        elif iline==2:
            comment2=line
        elif iline==3:
            atoms.nat=int(line.split()[0])
            poisson.shiftx=float(line.split()[1])
            poisson.shifty=float(line.split()[2])
            poisson.shiftz=float(line.split()[3])
        elif iline==4:
            poisson.ngpx=int(line.split()[0])
            poisson.hx=float(line.split()[1])
        elif iline==5:
            poisson.ngpy=int(line.split()[0])
            poisson.hy=float(line.split()[2])
        elif iline==6:
            poisson.ngpz=int(line.split()[0])
            poisson.hz=float(line.split()[3])
            poisson.allocate()
        elif iline<=atoms.nat+6:
            iatom=int(line.split()[0])
            atoms.sat.append(get_atomic_symbol(iatom))
            atoms.rat.append([])
            atoms.qat.append(float(line.split()[1]))
            atoms.rat[-1].append(float(line.split()[2]))
            atoms.rat[-1].append(float(line.split()[3]))
            atoms.rat[-1].append(float(line.split()[4]))
            if iline==atoms.nat+6:
                igpx=0
                igpy=0
                igpz=0
        else:
            num=len(line.split())
            if num>0:
                poisson.rho[igpz][igpy][igpx]=float(line.split()[0])
                igpx,igpy,igpz=cube_advance_index(poisson,igpx,igpy,igpz)
            if num>1:
                poisson.rho[igpz][igpy][igpx]=float(line.split()[1])
                igpx,igpy,igpz=cube_advance_index(poisson,igpx,igpy,igpz)
            if num>2:
                poisson.rho[igpz][igpy][igpx]=float(line.split()[2])
                igpx,igpy,igpz=cube_advance_index(poisson,igpx,igpy,igpz)
            if num>3:
                poisson.rho[igpz][igpy][igpx]=float(line.split()[3])
                igpx,igpy,igpz=cube_advance_index(poisson,igpx,igpy,igpz)
            if num>4:
                poisson.rho[igpz][igpy][igpx]=float(line.split()[4])
                igpx,igpy,igpz=cube_advance_index(poisson,igpx,igpy,igpz)
            if num>5:
                poisson.rho[igpz][igpy][igpx]=float(line.split()[5])
                igpx,igpy,igpz=cube_advance_index(poisson,igpx,igpy,igpz)
    return atoms,poisson
#*****************************************************************************************
def cube_write(filename,atoms,ngpx,ngpy,ngpz,arr,hx,hy,hz,frmt1,frmt2):
    f=open(filename,"w")
    f.write("\n")
    f.write("\n")
    f.write(frmt1 % (int(atoms.nat),0.0,0.0,0.0))
    f.write(frmt1 % (ngpx,hx,0.0,0.0))
    f.write(frmt1 % (ngpy,0.0,hy,0.0))
    f.write(frmt1 % (ngpz,0.0,0.0,hz))
    for iat in range(int(atoms.nat)):
        iatom=get_atomic_number(atoms.sat[iat])
        x=atoms.rat[iat][0]
        y=atoms.rat[iat][1]
        z=atoms.rat[iat][2]
        f.write("%5d%13.5f%13.6f%13.6f%13.6f\n" % (iatom,atoms.qat[iat],x,y,z))
    item=0
    for igpx in range(ngpx):
        for igpy in range(ngpy):
            for igpz in range(ngpz):
                f.write(frmt2 % arr[igpz][igpy][igpx])
                item+=1
                if(item==6 or (igpx==ngpx-1 and igpy==ngpy-1 and igpz==ngpz-1)):
                    item=0
                    f.write("\n")
    f.close()
#*****************************************************************************************
