import psana
import sys
from xfel.cxi.cspad_ana import cspad_tbx
import numpy as np
from scitbx.array_family import flex
runno = int(sys.argv[1])
dataset_name = "exp=mfxls4916:run=%d:idx"%(runno)

ds = psana.DataSource(dataset_name)
#address = "MfxEndstation-0|Rayonix-0"
#src = psana.Source('DetInfo(%s)'%address)
# libtbx.python rund.py 23
# detnames exp=mfxls4916:run=23

from matplotlib import pyplot as plt
plt.ion()
icount=1
for run in ds.runs():
  times = run.times()

  for i, t in enumerate(times):
    evt = run.event(t)
    ts = cspad_tbx.evt_timestamp((t.seconds(),t.nanoseconds()/1e6))
    print icount,ts,
    icount+=1
    for key in evt.keys():
      if key.alias() != "":

        print key.alias(),
      if 0 and 'FEE-SPEC0' in str(key):
        print key.src(),
        d = evt.get(key.type(), key.src())
        plt.cla()
        plt.plot(range(len(d.hproj())), d.hproj(), '-')
        plt.draw()
        plt.pause(0.1)
      if "FEE_Spec" in key.alias():
        #from IPython import embed; embed()
        d = evt.get(key.type(), key.src())
        thisdata = d.data16()
        plt.cla()
        plt.imshow(thisdata)
        plt.draw()
        plt.pause(0.01)
    print
