def Compresion(seq, d):
    result = []
    T = len(seq) // d
    for j in range(d):
        result.append(sum(seq[j + i*d] for i in range(T)))
    return result
set_3 = set()
with open("9.txt") as f:
    for line in f:
        seq = list(map(int, line.strip().split()))
        set_3.add(tuple(Compresion(seq, 3)))

set_9 = set()
with open("33.txt") as f:
    for line in f:
        seq = list(map(int, line.strip().split()))
        set_9.add(tuple(Compresion(seq, 3)))
        if tuple(Compresion(seq, 3)) in set_3:
            print("Found a match in 33.txt:", seq)

print("Number of unique compressions in 9.txt:", len(set_3))
for comp in set_9:
    print(' '.join(map(str, comp)))

