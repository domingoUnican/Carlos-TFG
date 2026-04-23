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

psd_values = dict()
pairs_9 = set()
number = 0
with open("9.txt") as f:
    for line in f:
        seq = list(map(int, line.strip().split()))
        L = len(seq)
        psd = PSD(seq)
        psd_values[tuple(psd)] = seq
        if tuple(50 - psd) in psd_values:
            number += 1
            pairs_9.add((tuple(seq), tuple(psd_values[tuple(50 - psd)])))

pairs_33 = set()
print("Number of pairs in 9.txt:", number)
number = 0
with open("33.txt") as f:
    for line in f:
        seq = list(map(int, line.strip().split()))
        L = len(seq)
        psd = PSD(seq)
        psd_values[tuple(psd)] = seq
        if tuple(50 - psd) in psd_values:
            number += 1
            pairs_33.add((tuple(seq), tuple(psd_values[tuple(50 - psd)])))
print("Number of pairs in 9.txt:", len(pairs_9))
print("Number of pairs in 33.txt:", len(pairs_33))
compressions_3 = set()
for psd1, psd2 in pairs_9:
    seq1 = Compresion(psd1, 3)
    seq2 = Compresion(psd2, 3)
    compressions_3.add((tuple(seq1), tuple(seq2)))
print("Number of pairs of compressions with d=3:", len(compressions_3))

candidates = set()
for a,b in compressions_3:
    candidates.add(a)
    candidates.add(b)

def is_less(seq):
    for a in candidates:
        if all(x <= y for x, y in zip(seq, a)):
            return True
    return False

def legendre_pairs(seq, pos):
    if pos ==6:
        for a in candidates:
            temp = Compresion(seq, 3)
            if all(x <= y for x, y in zip(temp, a)):
                seq[6] = a[0]-temp[0]
                seq[7] = a[1]-temp[1]
                seq[8] = a[2]-temp[2]
                if np.max(PSD(seq)) <= 50:
                    yield seq
    if pos <6: 
        for val in range(12):
            temp = list(seq)
            temp[pos] = val
            if is_less(Compresion(temp, 3)):
                yield from legendre_pairs(temp, pos + 1)

psd_values = dict()
pairs_9 = set()
for seq in legendre_pairs([0]*9, 0):
    temp = (PSD(seq))
    psd_values[tuple(temp)] = seq
    if tuple(50 - temp) in psd_values:
        pairs_9.add((tuple(seq), tuple(psd_values[tuple(50 - temp)])))
print("Number of pairs in 9.txt:", len(pairs_9))

candidate_3 = set()
for a,b in pairs_9:
    seq1 = Compresion(a, 3)
    seq2 = Compresion(b, 3)
    candidate_3.add(tuple(seq1))
    candidate_3.add(tuple(seq2))
unico_33 = set()
unico_3 = set()
for a,b in pairs_33:
    seq1 = Compresion(a, 3)
    seq2 = Compresion(b, 3)
    if tuple(seq1) in candidate_3 and tuple(seq2) in candidate_3:
        unico_33.add((tuple(seq1), tuple(seq2)))
        unico_3.add((tuple(seq1), tuple(seq2)))
print("Numer of unique pairs in 33.txt that match candidates from 9.txt:", len(unico_33))
print("number of elements in unico_3:", len(unico_3))
unico_9 = set()
for a,b in pairs_9:
    seq1 = Compresion(a, 3)
    seq2 = Compresion(b, 3)
    if (tuple(seq1), tuple(seq2)) in unico_3:
        unico_9.add((tuple(seq1), tuple(seq2)))

    if (tuple(seq2), tuple(seq1)) in unico_3:
        unico_9.add((tuple(seq2), tuple(seq1)))

print("Number of unique", len(unico_9))
pickle.dump(unico_9, open("unico_9.pkl", "wb"))
pickle.dump(unico_33, open("unico_33.pkl", "wb"))

