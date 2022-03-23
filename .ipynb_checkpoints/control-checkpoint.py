import xobjects as xo
import xtrack as xt
import xpart as xp
import json
from cpymad.madx import Madx
import NAFFlib
from math import modf

p0c = 7000e9
i_mo = -350
ctx_cpu = xo.ContextCpu()

if False:
    mad = Madx()
    mad.option(echo=False)
    mad.call('lhcandrea5.madx')
    mad.use(sequence="lhcb1")
    line = xt.Line.from_madx_sequence(mad.sequence['lhcb1'],
    deferred_expressions=True)
    with open('hl_line.json', 'w') as fid:
        json.dump(line.to_dict(), fid, cls=xo.JEncoder)
        
with open('hl_line.json', 'r') as fid:
        loaded_dct = json.load(fid)
line = xt.Line.from_dict(loaded_dct)
    
particle_0 = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, p0c=p0c, x=1e-3, y=1e-3)

print(f"\n\nBefore changing Octupole Current, i_mo = {line.vars['i_mo']._get_value()}\n\n")
tracker_normal = xt.Tracker(_context=ctx_cpu, line=line)
ref = tracker_normal.find_closed_orbit(particle_0)
print(f'\n\nref.x = {ref.x}, ref.y = {ref.y}, ref.zeta = {ref.zeta} \n\n')


line.vars['i_mo'] = i_mo
print(f"\n\nAfter changing Octupole Current, i_mo = {line.vars['i_mo']._get_value()}\n\n")
tracker_normal = xt.Tracker(_context=ctx_cpu, line=line)
ref = tracker_normal.find_closed_orbit(particle_0)
print(f'\n\nref.x = {ref.x}, ref.y = {ref.y}, ref.zeta = {ref.zeta} \n\n')