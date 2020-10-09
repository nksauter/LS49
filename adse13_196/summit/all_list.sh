#!/bin/bash -f

for m in 353397 353395 353370 353326 353325 353323 353319 353715 353716 353717 353718 353719 353720 353721;
do echo ${m};
cd ${m};
libtbx.python ../cctbx_deployment/opt/cctbx/modules/LS49/adse13_196/weather2.py show_plot=False;
cd -;
echo "...";
done

for m in 354113 354114 354115 354116 354117 354118 354119 354175 354124 354125 354126 354127 354128 354129;
do echo ${m};
cd ${m};
libtbx.python ../cctbx_deployment/opt/cctbx/modules/LS49/adse13_196/weather2.py show_plot=False;
cd -;
echo "...";
done

