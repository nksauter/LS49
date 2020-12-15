
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
import pandas
import os
from dxtbx.model.experiment_list import ExperimentListFactory
from LS49 import special_util

"""
Regarding Cytochrome SwissFEL data Code name: boop

Simulates optimized diffraction models and compares with real measurements

After running stage 1 diffBragg we obtain optimal diffraction models
Model parameters are stored in a pandas dataframe pickle file, they include
  - Ncells_abc (anisotropic)
  - Amatrix
  - scale factor (per image)

Background is simulated separately and then optimized against the measurements
before being added to the simulated spots
"""

# TODO philify
exascale = False
num_devices = 5
ADU_PER_PHOTON = 9.481
oversample = 0
df = pandas.read_pickle("/global/cfs/cdirs/m3562/der/braggnanimous/best_betty_boopz.pkl")
outdir = "/global/cfs/cdirs/m3562/collab_view/best_boops_trial0"
number_of_sim = 15


if COMM.rank == 0:
    while os.path.exists(outdir):
        trial_number = int(outdir.split("_trial")[1])+1
        new_outdir = outdir.split("_trial")[0] + "_trial%d" % trial_number
        print("%s exists, making new directory %s" % (outdir, new_outdir))
        outdir = new_outdir
    os.makedirs(outdir)
outdir = COMM.bcast(outdir)

df['basename'] = [os.path.basename(f) for f in df.exp_name]
df["oversample"] = oversample
top_shots = ["top_%d.expt" % x for x in range(number_of_sim)]

from iotbx.reflection_file_reader import any_reflection_file
merge_file = "/global/cfs/cdirs/m3562/der/cyto_init_merge.mtz"
Fmerge = any_reflection_file(merge_file).as_miller_arrays()[0].as_amplitude_array()

if exascale:
    from simtbx.gpu import gpu_energy_channels
    simulator_F = gpu_energy_channels(deviceId=COMM.rank % num_devices)
else:
    simulator_F = Fmerge

shot_info = []
for i_shot, top_exp in enumerate(top_shots):
    if i_shot % COMM.size != COMM.rank:
        continue

    shot = df.query("basename=='%s'" % top_exp)

    print("Simulating shot %d" % i_shot)
    if exascale:
        if simulator_F.get_nchannels() == 0:
            F_P1 = Fmerge.expand_to_p1()
            simulator_F.structure_factors_to_GPU_direct_cuda(0, F_P1.indices(), F_P1.data())

    expt = ExperimentListFactory.from_json_file(shot.exp_name.values[0], check_format=True)[0]
    print("Simulating background")
    background = special_util.exascale_sim(DETECTOR=expt.detector, BEAM=expt.beam,
                                           background_total_flux=shot.total_flux.values[0],
                                           include_background=True, include_spots=False, device_Id=0)
    print("Loading experimental data")
    imgs = special_util.image_data_from_expt(expt)
    imgs /= ADU_PER_PHOTON

    print("fitting background")
    bkgrnd_fit = special_util.determine_bkgrnd_scale(imgs, background)
    if not bkgrnd_fit.success:
        bg_scale = 0.06  # rough
    else:
        bg_scale = bkgrnd_fit.x[0]  # optimized
    background *= bg_scale

    spots = special_util.spots_from_pandas(shot, simulator_F, oversample=oversample,
                    cuda=True, device_Id=0, time_panels=True,
                    njobs=1,
                    exascale=exascale)

    spotmask = spots < 1
    prefix = os.path.join( outdir, "best_boop_%d" % i_shot)
    spot_mask_name = "%s.mask" % prefix
    special_util.save_numpy_mask_as_flex(spotmask, spot_mask_name)

    img_file = "%s.h5" % prefix
    special_util.save_model_to_image(expt, spots+background, img_file, save_experiment_data=True)

    shot["bkgrnd_scale"] = bg_scale
    shot["img_file"] = img_file
    shot["mask_file"] = spot_mask_name
    shot_info.append(shot)

shot_info = COMM.reduce(shot_info)

if COMM.rank == 0:
    shot_info = pandas.concat(shot_info)
    shot_info.to_pickle(os.path.join(outdir, "shot_pandas.pkl"))
    print("Done with betty boop!")
