from __future__ import division, print_function
from psana import *
import matplotlib.pyplot as plt

ds = DataSource('exp=xpptut15:run=54:smd')
det = Detector('cspad')

for nevent,evt in enumerate(ds.events()):
    if nevent>=2: break
    # includes pedestal subtraction, common-mode correction, bad-pixel
    # suppresion, and uses geometry to position the multiple CSPAD panels
    # into a 2D image
    print('Fetching event number',nevent)
    img = det.image(evt)

    plt.imshow(img,vmin=-2,vmax=2)
    plt.show()
print('Done.')
