from __future__ import division, print_function
from six.moves import range
import psana
from xfel.cxi.cspad_ana import cspad_tbx
from scitbx.array_family import flex # implicit import
dataset_name = "exp=mfxls4916:run=17:idx"

ds = psana.DataSource(dataset_name)
#address = "MfxEndstation-0|Rayonix-0"
#src = psana.Source('DetInfo(%s)'%address)

from matplotlib import pyplot as plt
plt.ion()

for run in ds.runs():
  times = run.times()

  for i, t in enumerate(times):
    evt = run.event(t)
    ts = cspad_tbx.evt_timestamp((t.seconds(),t.nanoseconds()/1e6))
    print(ts, end=' ')
    for key in evt:
      if key.alias() in ['Rayonix', 'Jungfrau1M']:
        print(key.alias(), end=' ')
      if 'FEE-SPEC0' in str(key):
        print(key.src(), end=' ')
        d = evt.get(key.type(), key.src())
        plt.cla()
        plt.plot(range(len(d.hproj())), d.hproj(), '-')
        plt.draw()
        plt.pause(0.1)
    print()
