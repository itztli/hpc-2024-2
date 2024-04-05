#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

with open("spectrum-mpi-1.dat") as f:
    data = f.readlines()

sampling=1024
signal = []
f = []
i = 0
for line in data:
    if 1 <= i <= sampling/2:
        signal.append(float(line))
    i = i+1

N=sampling # esto solo ocurre con datos por segundo.
delta = 1.0/N
f_c = 1.0 / (2.0*delta) 

for i in np.arange(0,int(N/2)):
    n = (-N/2) + i    
    f.append(f_c + (n/(N*delta)))

print(f)
plt.plot(f,signal)
#plt.xlim([0,511])
plt.show()
