from __future__ import division
from psana import *
#ds = DataSource('exp=cxig3614:run=231:smd')
ds = DataSource('exp=cxig3614:run=231')
nevent = 0
for evt in ds.events():
    nevent+=1
    #if nevent==3: break
    print('Processed',nevent,'events.',evt)
