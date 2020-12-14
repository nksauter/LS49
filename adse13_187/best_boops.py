
from libtbx.mpi4py import MPI
COMM = MPI.COMM_WORLD
import pandas
import os
from dxtbx.model.experiment_list import ExperimentListFactory
from LS49 import special_util

exascale = False
num_devices = 1

ADU_PER_PHOTON = 9.481
oversample = 8
df = pandas.read_pickle("/global/cfs/cdirs/m3562/der/braggnanimous/best_betty_boopz.pkl")
outdir = "/global/cfs/cdirs/m3562/collab_view/best_boops_trial0"
if COMM.rank == 0:
    while os.path.exists(outdir):
        trial_number = int(outdir.split("_trial")[1])+1
        new_outdir = outdir.split("_trial")[0] + "_trial%d" % trial_number
        print("%s exists, making new directory %s" % (outdir, new_outdir))
        outdir = new_outdir
    os.makedirs(outdir)
outdir = COMM.bcast(outdir)

df['basename'] = [os.path.basename(f) for f in df.exp_name]
top_5 = ["top_%d.expt" % x for x in range(5)]


from iotbx.reflection_file_reader import any_reflection_file
merge_file = "/global/cfs/cdirs/m3562/der/cyto_init_merge.mtz"
Fmerge = any_reflection_file(merge_file).as_miller_arrays()[0].as_amplitude_array()

if exascale:
    from simtbx.gpu import gpu_energy_channels
    simulator_F = gpu_energy_channels(deviceId=COMM.rank % num_devices)
else:
    simulator_F = Fmerge

model_info = []
for i_shot, top_exp in enumerate(top_5):
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
    expt = ExperimentListFactory.from_json_file(shot.exp_name.values[0], check_format=True)[0]
    imgs = special_util.image_data_from_expt(expt)
    imgs /= ADU_PER_PHOTON

    print("fitting background")
    bkgrnd_fit = special_util.determine_bkgrnd_scale(imgs, background)
    if not bkgrnd_fit.success:
        background *= 0.06  # rough
    else:
        background *= bkgrnd_fit.x[0]  # scale the background

    spots = special_util.spots_from_pandas(shot, simulator_F, oversample=oversample,
                    cuda=True, device_Id=0, time_panels=True,
                    njobs=1,
                    exascale=exascale)

    spotmask = spots < 1
    prefix = os.path.join( outdir, "best_boop_%d" % i_shot)
    special_util.save_numpy_mask_as_flex(spotmask, "%s.mask" % prefix)

    img_file = "%s.h5" % prefix
    special_util.save_model_to_image(expt, spots+background, img_file, save_experiment_data=True)
    model_info.append(shot)
    shot["oversample"] = oversample
    shot["img_file"] = img_file

model_info = pandas.concat(model_info)
model_info.to_piklle("%s.pandas.pkl" % prefix)

print("Done with betty boop!")
