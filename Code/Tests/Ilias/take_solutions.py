from collections import defaultdict
import numpy as np
from itertools import product
import pickle
from math import gcd

def PSD(seq):
    return np.rint(np.abs(np.fft.fft(seq))**2)[1:]

def Compresion(seq, d):
    result = []
    T = len(seq) // d
    for j in range(d):
        result.append(sum(seq[j + i*d] for i in range(T)))
    return result

legendre_pairs = []
with open("99-pairs-3", 'r') as f:
    for l in f.readlines():
        seq = list(map(int, l.split()))
        s1, s2 = seq[:len(seq)//2], seq[len(seq)//2:]
        s1 = [(3 +i)//2 for i in s1]
        s2 = [(3 +i)//2 for i in s2]
        if np.all(PSD(s1) + PSD(s2) != 50):
            print("ERROR")
        legendre_pairs.append((s1, s2))
possibilities = set()
with open("99-pairs-11", "r") as f:
    for l in f.readlines():
        seq = list(map(int, l.split()))
        s1, s2 = seq[:len(seq)//2], seq[len(seq)//2:]
        possibilities.add(tuple(Compresion(s1,3)))
        possibilities.add(tuple(Compresion(s2,3)))

f = open("unique_3.txt", 'w')
others = set()
for s1, s2 in legendre_pairs:
    c1 = Compresion(s1, 3)
    c2 = Compresion(s2, 3)
    others.add(tuple(c1))
    others.add(tuple(c2))
    if tuple(c1) in possibilities and tuple(c2) in possibilities:
        f.write(" ".join(map(str, s1)) + "\n")
        f.write(" ".join(map(str, s2)) + "\n")
f.close()
others_normalized = set()
for c in others:
    is_added = False
    for d in range(1, 3):
        for s in range(3):
            if tuple([c[(d*i+s)%3]for i in range(3)]) in others_normalized:
                is_added = True
    if not is_added:
        others_normalized.add(c)
print("longitud:",len(others_normalized))

f = open ("unique_11.txt", 'w')
normalized_33 = set()
normalized_3 = set()
for s1, s2 in legendre_pairs:
        c1 = Compresion(s1, 3)
        c2 = Compresion(s2, 3)
        if not (tuple(c1) in others_normalized and tuple(c2) in others_normalized):
            continue
        is_added = False
        for d in range(1, len(s1)):
            if gcd(d, len(s1)) != 1:
                continue
            s1_decimated = [s1[(d*i)%len(s1)] for i in range(len(s1))]
            s2_decimated = [s2[(d*i)%len(s2)] for i in range(len(s2))]
            if tuple(s1_decimated) in normalized_33 and tuple(s2_decimated) in normalized_33:
                is_added = True
        if not is_added:
            normalized_33.add(tuple(s1))
            normalized_33.add(tuple(s2))
print("longitud:", len(normalized_33))
for s in normalized_33:
    f.write(" ".join(map(str, s)) + "\n")
f.close()
print("longitud:", len(set([tuple(Compresion(s, 3)) for s in normalized_33])))
normalized_11 = set()
g = open ("unique_3.txt", 'w')
with open("99-pairs-11", "r") as f:
    for l in f.readlines():
        seq = list(map(int, l.split()))
        s1, s2 = seq[:len(seq)//2], seq[len(seq)//2:]
        if tuple(Compresion(s1, 3)) in others_normalized and tuple(Compresion(s2, 3)) in others_normalized:
            g.write(" ".join(map(str, s1)) + "\n")
            g.write(" ".join(map(str, s2)) + "\n")
f.close()
            
