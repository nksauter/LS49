for p in 0.0064 0.0128 0.0256 0.0384 0.0512; do \ # loop over eta values in degrees
for q in 24 48 96; do \ # loop over NcellsAB assumes tetragonal or hexagonal
for r in 8 16 24; do \ # loop over NcellsC
sbatch $MODULES/LS49/LY99/spread1/trial8_2923740_series.sh ${p} ${q} ${r}; sleep 10; done; done; done
