#/bin/bash

libtbx.python step5_ang_misset_to_fine_detail.py fine LS49_integ_allrestr > allrestr_to_fine_detail.digest

libtbx.python step5_ang_misset_to_fine_detail.py coarse LS49_integ_allrestr > allrestr_to_coarse_detail.digest

libtbx.python step5_ang_misset_to_fine_detail.py fine LS49_integ_betarestr > betarestr_to_fine_detail.digest

libtbx.python step5_ang_misset_to_fine_detail.py coarse LS49_integ_betarestr > betarestr_to_coarse_detail.digest

libtbx.python step5_ang_misset_to_fine_detail.py fine LS49_integ_step5cori > norestr_to_fine_detail.digest

libtbx.python step5_ang_misset_to_fine_detail.py coarse LS49_integ_step5cori > norestr_to_coarse_detail.digest
