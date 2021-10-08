srun -n1 -c10 libtbx.python \
$MODULES/LS49/adse13_187/adse13_221/mcmc_runner.py \
trusted_mask=$BERNINA/trusted_Py3.mask \
cryst=/global/cfs/cdirs/m3562/der/braggnanimous/top8_newlam2/expers/rank0/stg1_top_0_0.expt \
refl=$BERNINA/split2b/split_0309.refl \
expt=$BERNINA/split_c/split_0648.expt \
output.label=mcmc3 output.index=0 model.mosaic_spread.value=0.01 \
model.Nabc.value=50,50,15 model.Nabc.hyperparameter=0.8 \
model.rot.refine=True \
model.cell.covariance="$BERNINA/../covariance_cytochrome_form.pickle" simplex.cycles=200 mcmc.cycles=1000

# ... a baseline test to exercise mcmc.

#cryst=/global/cscratch1/sd/nksauter/adse13_187/bernina/split2b/split_0309.expt \
#cryst=/global/cfs/cdirs/m3562/der/braggnanimous/top8_newlam2/expers/rank0/stg1_top_0_0.expt \
#model.Nabc.value=130,30,10
