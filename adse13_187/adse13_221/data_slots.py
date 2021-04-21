from __future__ import division
from collections import OrderedDict

application_slots = OrderedDict(
  [
    #Output 1.  Lunus pixel-assimilated image
    # 3D numpy array
    ("lunus_filtered_data", None),
    #Output 7. Experimental res-data
    # list of 2D numpy arrays
    ("exp_data", None),
  ]
)
