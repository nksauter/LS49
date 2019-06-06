from __future__ import division, print_function
"""Read in a 3000 x 3000 png image
1. Confirm it is greyscale
2. Convert to 1000 x 1000
3. Set pixels to max(3x3 box)
4. Output new png.
"""

from PIL import Image

im = Image.open("fig4_.png")
nw = im.convert(mode="L")
px = nw.load() # original image pixel access, 255=white, 0=black
stride = 3 # corresponds to input 3K x 3K, output 1K x 1K

ot = nw.resize((1000,1000),resample=Image.LANCZOS) # will not use pixel values in final output
pxout = ot.load() # output image pixel access
modulus = 2
for irow in range(modulus+1,1000-modulus-1):
  for icol in range(modulus+1,1000-modulus-1):
    srcset = []
    for ioffr in range(-modulus,modulus+1):
      for ioffc in range(-modulus,modulus+1):
        srcset.append(px[icol*stride+ioffc,irow*stride+ioffr])
    pxout[icol,irow]= min(srcset)
    #pxout[icol,irow]=10

ot.show()
ot.save("fig4_mod%1d.png"%modulus)

#from IPython import embed; embed()
