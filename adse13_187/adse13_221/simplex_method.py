from __future__ import division
from scitbx.array_family import flex
from scitbx.simplex import simplex_opt
from LS49.adse13_187.cyto_batch import multipanel_sim
import copy
from time import time

class simplex_detail:
  def __init__(selfOO,alt_crystal,Ncells_abc,host_runner,PP,n_cycles,s_cycles):
          selfOO.n_cycles = n_cycles
          selfOO.PP = PP
          selfOO.crnm = host_runner
          # count the full number of parameters
          selfOO.n = 0
          selfOO.iteration = 0
          selfOO.alt_crystal = alt_crystal
          selfOO.Ncells_abc = Ncells_abc
          initial_values = flex.double()
          selfOO.starting_simplex=[initial_values]
          bounding_values = flex.double()
          for key in selfOO.crnm.ref_params:
            selfOO.crnm.ref_params[key].accept()
            for label in selfOO.crnm.ref_params[key].display_labels:
              selfOO.n += 1
              print(label, selfOO.crnm.ref_params[key].chain[label][-1])
              initial_values.append(selfOO.crnm.ref_params[key].chain[label][-1])
            vals = selfOO.crnm.ref_params[key].generate_simplex_interval()
            for il, label in enumerate(selfOO.crnm.ref_params[key].display_labels):
              print(label, vals[il])
              bounding_values.append(vals[il])
          #print("there are %d parameters"%selfOO.n)
          selfOO.crnm.plot_all(selfOO.iteration+1,of=n_cycles)

          for ii in range(selfOO.n):
            vertex = copy.deepcopy(initial_values)
            vertex[ii] = bounding_values[ii]
            selfOO.starting_simplex.append(vertex)
          print(selfOO.starting_simplex)
          selfOO.optimizer = simplex_opt( dimension=selfOO.n,
                                        matrix  = selfOO.starting_simplex,
                                        evaluator = selfOO,
                                        monitor_cycle=20,
                                        max_iter=s_cycles-1,
                                        tolerance=1e-7)
          selfOO.x = selfOO.optimizer.get_solution()
          #print(list(selfOO.x))
          #input()

  def target(selfOO, vector):
          selfOO.iteration+=1
          ptr=0
          for key in selfOO.crnm.ref_params:
            n_values = len(selfOO.crnm.ref_params[key].display_labels)
            these_values = vector[ptr:ptr+n_values]
            selfOO.crnm.ref_params[key].set_proposal_from_simplex(these_values)
            ptr+=n_values

          BEG=time()
          if "cell" in selfOO.crnm.ref_params:
            selfOO.alt_crystal = selfOO.crnm.ref_params["cell"].get_current_crystal_model(selfOO.alt_crystal)
          if "rot" in selfOO.crnm.ref_params:
            selfOO.alt_crystal = selfOO.crnm.ref_params["rot"].get_current_crystal_model(selfOO.alt_crystal)
          if "ncells" in selfOO.crnm.ref_params:
            selfOO.Ncells_abc = selfOO.crnm.ref_params["ncells"].get_current_model()

          whitelist_only, TIME_BG, TIME_BRAGG, selfOO.crnm.exascale_mos_blocks = multipanel_sim(
        CRYSTAL=selfOO.alt_crystal, DETECTOR=selfOO.PP["detector"], BEAM=selfOO.PP["beam"],
        Famp = selfOO.crnm.gpu_channels_singleton,
        energies=selfOO.PP["energies"], fluxes=selfOO.PP["weights"],
        cuda=True,
        oversample=selfOO.PP["oversample"], Ncells_abc=selfOO.Ncells_abc,
        mos_dom=selfOO.PP["mosaic_spread_samples"], mos_spread=selfOO.crnm.parameters["etaa"].proposal,
        mos_aniso=(selfOO.crnm.parameters["etaa"].proposal,
                   selfOO.crnm.parameters["etab"].proposal,
                   selfOO.crnm.parameters["etac"].proposal),
        beamsize_mm=selfOO.PP["beamsize_mm"],
        profile=selfOO.PP["shapetype"],
        show_params=False,
        time_panels=False, verbose=selfOO.PP["verbose"],
        spot_scale_override=selfOO.PP["spot_scale"],
        include_background=False,
        mask_file=selfOO.PP["mask_array"], skip_numpy=True,
        relevant_whitelist_order=selfOO.crnm.relevant_whitelist_order
      )
          Rmsd,sigZ,LLG = selfOO.PP["Z"](kernel_model=whitelist_only, plot=False)
          #print ("Old NLL ",selfOO.crnm.llg_chain[-1], "NEW LLG",LLG, "diff",selfOO.crnm.llg_chain[-1] - LLG)
          for key in selfOO.crnm.ref_params:
              selfOO.crnm.ref_params[key].accept()
          selfOO.crnm.accept.append(1)
          selfOO.crnm.rmsd_chain.append(Rmsd); selfOO.crnm.sigz_chain.append(sigZ); selfOO.crnm.llg_chain.append(LLG)
          selfOO.crnm.plot_all(selfOO.iteration+1,of=selfOO.n_cycles)
          TIME_EXA = time()-BEG
          #print("\t\tExascale: time for Bragg sim: %.4fs; total: %.4fs" % (TIME_BRAGG, TIME_EXA))
          return LLG

