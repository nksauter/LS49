#!/bin/bash -f
export OLDCFSX=/global/cfs/cdirs/m3562/swissfel
export MIRCFSX=${SCRATCH}/cfsx/cdirs/m3562/swissfel
mkdir -p ${MIRCFSX}
mkdir -p ${MIRCFSX}/data/Cyt/data/dark
mkdir -p ${MIRCFSX}/raw/Cyt
mkdir -p ${MIRCFSX}/gainMaps/JF07T32V01
mkdir -p ${MIRCFSX}/../der/pedestals

for run in {795..805}
do
    echo $run
    # unassembled_file
    cp -pr ${OLDCFSX}/data/Cyt/data/dark/run_000${run}.JF07T32V01.h5 ${MIRCFSX}/data/Cyt/data/dark
    # raw_file
    cp -pr ${OLDCFSX}/raw/Cyt/run_000${run}.JF07T32V01.h5 ${MIRCFSX}/raw/Cyt
    # geom_file
    cp -pr ${OLDCFSX}/16M_bernina_backview_optimized_adu_quads.geom ${MIRCFSX}
    # pedestal_file
    cp -pr ${OLDCFSX}/../der/pedestals/pedestal_20191029_2033.JF07T32V01.res.h5 ${MIRCFSX}/../der/pedestals
    # gain_file
    cp -pr ${OLDCFSX}/gainMaps/JF07T32V01/gains.2019-06.h5 ${MIRCFSX}/gainMaps/JF07T32V01
    # beam_file
    cp -pr ${OLDCFSX}/raw/Cyt/run_000${run}.BSREAD.h5 ${MIRCFSX}/raw/Cyt

done
