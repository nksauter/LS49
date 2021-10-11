#!/bin/bash -f
export CFSX=${SCRATCH}/cfsx/cdirs/m3562/swissfel
for run in {795..805}
do
    echo $run
    libtbx.python ${MODULES}/cctbx_project/xfel/swissfel/jf16m_cxigeom2nexus.py \
    unassembled_file=${CFSX}/data/Cyt/data/dark/run_000${run}.JF07T32V01.h5 \
    raw_file=${CFSX}/raw/Cyt/run_000${run}.JF07T32V01.h5 \
    geom_file=${CFSX}/16M_bernina_backview_optimized_adu_quads.geom \
    detector_distance=163.9 \
    include_spectra=True \
    pedestal_file=${CFSX}/../der/pedestals/pedestal_20191029_2033.JF07T32V01.res.h5 \
    gain_file=${CFSX}/gainMaps/JF07T32V01/gains.2019-06.h5 \
    output_file=${BERNINA}/run_000${run}.JF07T32V01_temp_readout_master.h5 \
    raw=False \
    beam_file=${CFSX}/raw/Cyt/run_000${run}.BSREAD.h5

    libtbx.python ${MODULES}/cctbx_project/xfel/swissfel/sync_jf16m_geometry.py \
    nexus_master=${BERNINA}/run_000${run}.JF07T32V01_temp_readout_master.h5 \
    new_geometry=${BERNINA}/split_0000.expt \
    output_file=${BERNINA}/run_000${run}.JF07T32V01_readout_master.h5
done
