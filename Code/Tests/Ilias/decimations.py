import numpy as np
from itertools import product
import pickle

def PSD(seq):
    return np.rint(np.abs(np.fft.fft(seq))**2)[1:]

def Compresion(seq, d):
    result = []
    T = len(seq) // d
    for j in range(d):
        result.append(sum(seq[j + i*d] for i in range(T)))
    return result

with open("99-pairs-3", 'r') as f:
    for l in f.readlines():
        seq = list(map(int, l.split()))
        s1, s2 = seq[:len(seq)//2], seq[len(seq)//2:]
        print(PSD(s1) + PSD(s2))

            
