from __future__ import division, print_function
import psana
ds = psana.DataSource('exp=cxig3614:run=201')
det = psana.Detector('Ds1CsPad')
#det = psana.Detector('Fee_Orca_Spectrometer')
det = psana.Detector('CsPad140k')
##-----------------------------
# graphics initialization
import matplotlib.pyplot as plt
fig  = plt.figure(figsize=(13,12), dpi=80, facecolor='w', edgecolor='w', frameon=True)
axim = fig.add_axes((0.05,  0.03, 0.87, 0.93))
axcb = fig.add_axes((0.923, 0.03, 0.02, 0.93))
plt.ion() # do not hold control on show() in the event loop
##-----------------------------

for nevent,evt in enumerate(ds.events()):
    #if nevent>10 : break
    calib_array = det.calib(evt)
    img = det.image(evt)
    #from IPython import embed; embed()
    if img is None:
      print(nevent,"None"); continue
    ##-----------------------------
    # graphics
    axim.cla()                                         # clear image axes if necessary...
    ave, rms = img.mean(), img.std()                   # evaluate average intensity and rms over image pixels
    imsh = axim.imshow(img, interpolation='nearest', aspect='auto', origin='upper',\
                       vmin=ave-1*rms, vmax=ave+2*rms) # make object which produces image in axes
    colb = fig.colorbar(imsh, cax=axcb)                # make colorbar object associated with figure
    fig.canvas.set_window_title('Event %d' % nevent)   # set window title
    fig.canvas.draw()                                  # redraw canvas
    plt.show()                                         # show plot (need it for the 1st plot only)

plt.ioff() # hold control on show() at the end
plt.show()
##-----------------------------
