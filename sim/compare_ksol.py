from __future__ import division

from scitbx.array_family import flex

pdb_lines = open("1m2a.pdb","r").read()

if __name__=="__main__":
  from LS49.sim.util_fmodel import gen_fmodel
  GF = gen_fmodel(resolution=10.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.73424)
  ksol=0.01*flex.double([float(u) for u in xrange(70)])
  sumarray = []
  for u in ksol:
    GF.params2.fmodel.k_sol = u
    amplitudes = GF.get_amplitudes()
    asum = flex.sum(amplitudes.data())
    #for i in xrange(amplitudes.size()):
      #print i,amplitudes.indices()[i],amplitudes.data()[i]
    print "ksol = %f sum is %f"%(u,asum)
    sumarray.append(asum)
  from matplotlib import pyplot as plt
  plt.plot(ksol,sumarray,'r.')

  GF = gen_fmodel(resolution=7.0,pdb_text=pdb_lines,algorithm="fft",wavelength=1.73424)
  ksol=0.01*flex.double([float(u) for u in xrange(70)])
  sumarray = []
  for u in ksol:
    GF.params2.fmodel.k_sol = u
    amplitudes = GF.get_amplitudes()
    asum = flex.sum(amplitudes.data())
    #for i in xrange(amplitudes.size()):
      #print i,amplitudes.indices()[i],amplitudes.data()[i]
    print "ksol = %f sum is %f"%(u,asum)
    sumarray.append(asum*0.7) # arbitrary scaling factor for better looking plot
  plt.plot(ksol,sumarray,'b.')
  plt.show()


  plt.show()
