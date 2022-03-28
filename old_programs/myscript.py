#source /home/afornara/py/2022_03_04/miniconda/bin/activate 

import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp
import json
import pandas as pd
from cpymad.madx import Madx
#from matplotlib import pyplot as plt
import NAFFlib
from math import modf
#from scipy.constants import physical_constants

#input
i_mo = 350
p0c = 7000e9
normal_emitt_x = 3e-6
normal_emitt_y = 2.6e-6
#proton_mass = physical_constants['proton mass energy equivalent in MeV'][0]/1000 

#ctx_gpu = xo.ContextCupy()
ctx_cpu = xo.ContextCpu()

if False:
    mad = Madx()
    mad.option(echo=False)
    mad.call('andrea.madx')
    mad.use(sequence="lhcb1")
    line = xt.Line.from_madx_sequence(mad.sequence['lhcb1'],
                                  deferred_expressions=True
                                  )
    with open('line.json', 'w') as fid:
        json.dump(line.to_dict(), fid, cls=xo.JEncoder)

else:
    with open('line.json', 'r') as fid:
        loaded_dct = json.load(fid)
    line = xt.Line.from_dict(loaded_dct)

particle_0 = xp.Particles(_context=ctx_cpu, p0c=p0c)

line.vars['i_mo'] = i_mo
line.vars['cmrs.b1_sq']=0.5*line.vars['cmrs.b1_sq']
assert np.isclose(line.vars['i_mo']._get_value(),350)

tracker = xt.Tracker(_context=ctx_cpu, line=line)
tw = tracker.twiss(particle_0)

betx_at_ip3 = tw['betx'][0]
bety_at_ip3 = tw['bety'][0]
#assert np.isclose(tw['betx'][0],120.29057)

sigma_x = np.sqrt(betx_at_ip3*normal_emitt_x/(particle_0.gamma0*particle_0.beta0))
sigma_y = np.sqrt(bety_at_ip3*normal_emitt_y/(particle_0.gamma0*particle_0.beta0))
#assert np.isclose(sigma_x,219.93e-6)
#assert np.isclose(sigma_x,219.93e-6)

p0 = tw['particle_on_co']

N_particles=13

particles = xp.build_particles(_context=ctx_cpu, 
                               particle_ref=p0, 
                               x=[mysigma*sigma_x[0] for mysigma 
                                 in np.linspace(0.2,2,N_particles)])

#particles = xp.build_particles(_context=ctx_cpu, particle_ref=p0, x=[p0.x-1.8*sigma_x,p0.x-1.5*sigma_x,p0.x-1.2*sigma_x,p0.x-0.9*sigma_x,p0.x-0.6*sigma_x,p0.x-0.3*sigma_x,p0.x+1e-5,p0.x+0.3*sigma_x,p0.x+0.6*sigma_x,p0.x+0.9*sigma_x,p0.x+1.2*sigma_x,p0.x+1.5*sigma_x,p0.x+1.8*sigma_x])


N=1000
#I'm tracking 11 particles in this case
n_turns=1
my_result = {}
for ii in ['x','px','y','py','zeta','delta']:
    my_result[ii] = []

for ii in range(N):
    tracker.track(particles, num_turns=n_turns,turn_by_turn_monitor=False)
    for jj in ['x','px','y','py','zeta','delta']:
        my_result[jj].append(getattr(particles,jj).copy())
        
for jj in ['x','px','y','py','zeta','delta']:
    my_result[jj] = np.array(my_result[jj])
    #first index is turn second index is particle index, my_result['x'][:,1]
    #are all turns for the second particle


        



qx = []
x0 = []
for ii in range(N_particles):
    qx.append(NAFFlib.get_tune(my_result['x'][:,ii]))
    x0.append(my_result['x'][0,ii])
    

#Do we obtain the same result with a new tracker per particle?


#particles = xp.build_particles(_context=ctx_cpu, particle_ref=p0, x=[p0.x-1.8*sigma_x,p0.x-1.5*sigma_x,p0.x-1.2*sigma_x,p0.x-0.9*sigma_x,p0.x-0.6*sigma_x,p0.x-0.3*sigma_x,p0.x+1e-5,p0.x+0.3*sigma_x,p0.x+0.6*sigma_x,p0.x+0.9*sigma_x,p0.x+1.2*sigma_x,p0.x+1.5*sigma_x,p0.x+1.8*sigma_x])

if False:
    tracker2 = xt.Tracker(_context=ctx_cpu, line=line)
    tw2 = tracker2.twiss(particle_0)
    part_2 = tw2['particle_on_co']
    part_2.x = part_2.x+1.5*sigma_x

    x_part2 = np.zeros(N)


    for ii in range(N):
        tracker2.track(part_2, num_turns=n_turns,turn_by_turn_monitor=False)
        x_part2[ii] = part_2.x



stop 
x_0=np.zeros(N)       
x_1=np.zeros(N)
x_2=np.zeros(N)
x_3=np.zeros(N)
x_4=np.zeros(N)
x_5=np.zeros(N)
x_6=np.zeros(N)
x_7=np.zeros(N)
x_8=np.zeros(N)
x_9=np.zeros(N)
x_10=np.zeros(N)
x_11=np.zeros(N)
x_12=np.zeros(N)
for kk in range(N):
    x_0[kk]=x[0][kk]
    x_1[kk]=x[1][kk]
    x_2[kk]=x[2][kk]
    x_3[kk]=x[3][kk]
    x_4[kk]=x[4][kk]
    x_5[kk]=x[5][kk]
    x_6[kk]=x[6][kk]
    x_7[kk]=x[7][kk] 
    x_8[kk]=x[8][kk]
    x_9[kk]=x[9][kk]
    x_10[kk]=x[10][kk]
    x_11[kk]=x[11][kk]
    x_12[kk]=x[12][kk]
        #px[jj][ii] = ctx_cpu.nparray_from_context_array(particles.px)[jj]  






q0 = NAFFlib.get_tune(x_0)
q1 = NAFFlib.get_tune(x_1)
q2 = NAFFlib.get_tune(x_2)
q3 = NAFFlib.get_tune(x_3)
q4 = NAFFlib.get_tune(x_4)
q5 = NAFFlib.get_tune(x_5)
q6 = NAFFlib.get_tune(x_6)
q7 = NAFFlib.get_tune(x_7)
q8 = NAFFlib.get_tune(x_8)
q9 = NAFFlib.get_tune(x_9)
q10 = NAFFlib.get_tune(x_10)
q11 = NAFFlib.get_tune(x_11)
q12 = NAFFlib.get_tune(x_12)

print('q0 = ',q0,'q1 =',q1,'q2 =', q2, ', q3 =', q3, ', q4 =', q4,', q5 =', q5,'q6',q6,'q7 =',q7,'q8 =', q8, ', q9 =', q9, ', q10 =', q10,', q11 =', q11, 'q12 = ',q12)
                                                                
                                                               








