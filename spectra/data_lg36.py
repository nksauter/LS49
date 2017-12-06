#source /reg/g/psdm/etc/psconda.sh
from psana import *
ds = DataSource('exp=cxig3614:run=210')
det = Detector('Ds1CsPad')
#det = Detector('Fee_Orca_Spectrometer')
ebeamDet = Detector('EBeam')
for nevent,evt in enumerate(ds.events()):
    # includes pedestal subtraction, common-mode correction, bad-pixel
    # suppresion, and returns an "unassembled" 3D array of cspad panels
    calib_array = det.calib(evt)
    # this is the same as the above, but also uses geometry to
    # create an "assembled" 2D image (including "fake pixels" in gaps)
    img = det.image(evt)
    if img is None:
      print 'None',nevent
      continue
    #break
    ebeam = ebeamDet.get(evt)

    print nevent, ebeam.ebeamPhotonEnergy()
    #if img.mean() < 10: continue
    import matplotlib.pyplot as plt
    #from IPython import embed; embed()
    plt.imshow(img,vmin=-2,vmax=2)
    #plt.imshow(img, interpolation='nearest')
    plt.show()

