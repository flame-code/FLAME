class Atoms:
    def __init__(self):
        self.sat=[]
        self.rat=[]
        self.fat=[]
        self.bemoved=[]
        self.qat=[]
        self.nat=-1
        self.nmolecule=-1
        self.cellvec=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
        self.boundcond="unknown"
        self.epot=1.E100
        self.qtot=0.0
        self.coordinates="Cartesian"
        self.pattern=-1
