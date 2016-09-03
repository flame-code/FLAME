import random
import math 
print ( 'pattern module imported')
def gen_pattern(x_cm, y_cm, z_cm, rcov_max,n_m,pat_num):
    dx=random.uniform(3,5)*rcov_max
    dy=random.uniform(3,5)*rcov_max
    dz=random.uniform(3,5)*rcov_max
    for t in range(1,n_m):
        if pat_num==1:
            x_cm+=t*dx
            #print x_cm, y_cm, z_cm
            return x_cm, y_cm, z_cm
        elif pat_num==2:
            r_circle=random.uniform(4,6)*rcov_max
            x_cm=t*dx*math.cos(pi/2.)
            y_cm=dy*math.sin(pi/2.)
            return x_cm, y_cm, z_cm
        else: break
