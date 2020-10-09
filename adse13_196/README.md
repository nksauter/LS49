Quick test to see if simulation results match up to a reference:
libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/compare.py testdir=??? refdir=??? maxdelta=1 maxcount=100
Furthermore, the specific reference used for the ADSE13-196 report was:
refdir=$(libtbx.find_in_repositories ls49_big_data)/adse13_196/report1
Executing within the image directory you can avoid supplying the dir keywords:
libtbx.python $(libtbx.find_in_repositories LS49)/adse13_196/compare.py maxdelta=1 maxcount=100

