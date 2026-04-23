from collections import defaultdict
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


values = set()
with open("7137_solution_strings.txt", 'r') as f:
    for l in f.readlines():
        seq = int(l.strip())
        values.add(seq)

lp = defaultdict(list)
with open("b9.w200.res", 'r') as f:
    for l in f.readlines():
        a = l.strip().split()[0]
        if a not in values:
            continue
        seq = list(map(int, l.strip().split())[1:])
        lp[a].append(seq)
with open("a9.w200.res", 'r') as f:
    for l in f.readlines():
        a = l.strip().split()[0]
        if int(a) not in values:
            continue
        seq = list(map(int, l.strip().split()[1:]))
        lp[a].append(seq)

psd_values = dict()
legendre_pairs = []
for k,v in lp.items():
    for seq in v:
        psd_values[tuple(PSD(seq))] = seq
        if tuple(200-PSD(seq)) in psd_values:
            legendre_pairs.append((seq, psd_values[tuple(200-PSD(seq))]))

f = open("../99-pairs-11", 'w')
for a, b in legendre_pairs:
    s1 = [(11+i)//2 for i in a]
    s2 = [(11+i)//2 for i in b]
    f.write(" ".join(map(str, s1+s2)) + "\n")
    
