class Atoms:
    def __init__(self):
        self.sat=[]
        self.rat=[]
        self.posred=[]
        self.fat=[]
        self.vat=[]
        self.bemoved=[]
        self.qat=[]
        self.nat=-1
        self.nmolecule=-1
        self.cellvec=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
        self.qpm=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
        self.dpm=[0.0,0.0,0.0]
        self.boundcond="unknown"
        self.epot=1.E100
        self.qtot=0.0
        self.coordinates="Cartesian"
        self.pattern=-1
        self.units_length_io='atomic'
        self.cell_present=False
        self.epot_present=False
        self.fat_present=False
        self.qtot_present=False
        self.bemoved_present=False
        self.vat_present=False
        self.dpm_present=False
        self.qpm_present=False
