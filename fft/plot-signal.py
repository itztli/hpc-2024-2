#!/usr/bin/env python3

import matplotlib.pyplot as plt

with open("signal.dat") as f:
    data = f.readlines()

signal = []
for line in data:
    signal.append(float(line))

plt.plot(signal)
plt.plot(signal,'or')


plt.show()
