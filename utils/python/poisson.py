import numpy as np

class Poisson:
    def __init__(self): #,ngpx,ngpy,ngpz):
        self.ngpx=0
        self.ngpy=0
        self.ngpz=0
        self.hx=0.0
        self.hy=0.0
        self.hz=0.0
        self.rho=None
    def allocate(self):
        self.rho=np.empty([self.ngpz,self.ngpy,self.ngpx])
        #=np.empty() #([ngpz,ngpy,ngpx])
        #self.pot=np.empty([ngpz,ngpy,ngpx])
