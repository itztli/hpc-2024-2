#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import wave
import sys


spf = wave.open("mono.wav", "r")

# Extract Raw Audio from Wav File
signal = spf.readframes(-1)
#print(signal)
signal = np.frombuffer(signal,dtype=np.int16)


# If Stereo
if spf.getnchannels() == 2:
    print("Just mono files")
    sys.exit(0)

#for data in signal:
#    print(data)

plt.figure(1)
plt.title("Signal Wave...")
plt.plot(signal)
plt.show()
