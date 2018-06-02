from __future__ import division
from psana import *
ds = DataSource('exp=cxig3614:run=231')

ebeamDet = Detector('EBeam')

for nevent,evt in enumerate(ds.events()):
    ebeam = ebeamDet.get(evt)
    if ebeam is None: continue
    print(nevent, ebeam.ebeamPhotonEnergy())
    #break
