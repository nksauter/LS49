from __future__ import division
from collections import OrderedDict

application_slots = OrderedDict(
  [
    #Output 1.  Lunus pixel-assimilated image
    # 3D numpy array
    ("lunus_filtered_data", None),

    #Output 2.  Modified copy of the Lunus image, with 1st-order Taylor shoeboxes
    ("plane_shoeboxes", None),

    #Output 3. analyze proposal and add background
    ("bragg_plus_background", None),

    #Output 4. renormalize the proposal
    ("renormalize_bragg_plus_background", None),

    #Output 5. Mockup simulation laid on top of 1st-Taylor background
    ("sim_mock", None),

    #Output 6. Figure the Z-plot
    ("Z_plot", None),

    #Output 7. Experimental res-data
    # list of 2D numpy arrays
    ("exp_data", None),
  ]
)
