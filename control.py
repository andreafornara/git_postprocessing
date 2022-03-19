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

p0c = 7000e9

ctx_cpu = xo.ContextCpu()

mad = Madx()
mad.option(echo=False)
mad.call('lhcandrea4.madx')
mad.use(sequence="lhcb1")
line = xt.Line.from_madx_sequence(mad.sequence['lhcb1'],
deferred_expressions=True)
with open('hl_line.json', 'w') as fid:
    json.dump(line.to_dict(), fid, cls=xo.JEncoder)
if False:
    with open('hl_line.json', 'r') as fid:
        loaded_dct = json.load(fid)
    line = xt.Line.from_dict(loaded_dct)

particle_0 = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=p0c, x=1e-3, y=1e-3)
print(line.vars['i_mo']._get_value())
tracker_normal = xt.Tracker(_context=ctx_cpu, line=line)
ref = tracker_normal.find_closed_orbit(particle_0)
print(ref.x,ref.y,ref.zeta)
