from __future__ import division
"""
Find correlation per-event for two performance times, either accessing
data from scratch disk (all-flash) or community file system (spinning disk).
"""
run1={}
run2={}
with open("/pscratch/sd/n/nksauter/adse13_249/work/112721/perftimes","r") as F:
    for line in F:
        tokens=line.strip().split()
        task=int(tokens[0])
        time=float(tokens[1][:-1])
        run1[task]=time
print(run1)
with open("/pscratch/sd/n/nksauter/adse13_249/work/112746/perftimes","r") as F:
    for line in F:
        tokens=line.strip().split()
        task=int(tokens[0])
        time=float(tokens[1][:-1])
        run2[task]=time

task1 = [run1[T] for T in range(3000)]
task2 = [run2[T] for T in range(3000)]

from matplotlib import pyplot as plt
plt.plot(task1, task2, "r.")
plt.show()
