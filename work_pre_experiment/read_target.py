from __future__ import print_function
from scitbx.array_family import flex
import pickle

V = open("data.pickle","rb")
while 1:
  image = pickle.load(V)
  print (image["image"])
  highcc = flex.double(image["cc"]) > 0.70
  if highcc.count(True)<4: continue
  for i in range(len(image["cc"])):
    if image["cc"][i]<0.7: continue
    print ("CC>70%%: %20s %5.2f"%(image["millers"][i],
           image["cc"][i]))

