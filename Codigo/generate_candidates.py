import numpy as np
from itertools import product
from collections import defaultdict

import sys

def PAF(seq):
    """
    This is fast version using the Fast Fourier Transform, 
    but the first element is remove because it is always len(seq).
    """
    temp = np.fft.fft(seq)
    temp *= np.roll(temp[::-1],1) # This is the formula for the PAF
    temp = np.fft.ifft(temp)
    return tuple(np.array(np.rint(temp),int)[1:])


def Compresion(seq, m):
    """
    This is function that compress a sequence, summing the elements
    """
    L = len(seq)
    result = []
    for i in range(m):
        result.append(sum((seq[i+m*j] for j in range(L//m))))
    return tuple(result)

def domina(x,y):
    return all(xi <= yi for xi, yi in zip(x,y))


def camino_mas_largo(candidatos):
    distancias = defaultdict(lambda: 0)
    predecesores = dict()
    max_distancia = 0
    for i in candidatos:
        for j in candidatos:
            if i!=j and domina(i,j):
                if distancias[i]+1 > distancias[j]:
                    distancias[j] = max(distancias[j], distancias[i]+1)
                    predecesores[j] = i
                max_distancia = max(max_distancia, distancias[i] + 1, distancias[j])
    print("Max distance: ", max_distancia)
    for k, v in distancias.items():
        if v == max_distancia:
            resultado = []
            while k in predecesores:
                resultado.append(k)
                k = predecesores[k]
            if len(resultado) != max_distancia:
                print("Error: ", resultado, max_distancia)
            return resultado, set(candidatos) - set(resultado)

    return [], set()

def legendre_sequence(n):
    """
    Generate the Legendre sequence of length n.
    The Legendre sequence is defined as follows:
    - a[i] = 1 if i is a quadratic residue modulo n
    - a[i] = -1 if i is a non-quadratic residue modulo n
    - a[0] = 0

    Parameters
    ----------
    n : int
        Length of the sequence.

    Returns
    -------
    list[int]
        The Legendre sequence of length n.
    """
    legendre_seq = [5] * n
    for i in range(1, n):
        legendre_seq[i] = 3 if pow(i, (n - 1) // 2, n) == 1 else 6
    return legendre_seq


def caminos_grande(candidatos):
    """
    Ordena los candidatos y los divide en una lista de caminos (cadenas),
    donde en cada camino se cumple domina(camino[i], camino[i+1]).

    Parameters
    ----------
    candidatos : iterable[tuple[int, ...]]
        Lista (o conjunto) de candidatos.

    Returns
    -------
    list[list[tuple[int, ...]]]
        Partición en caminos dominantes.
    """
    ordenados = sorted(candidatos)
    caminos = []

    for candidato in ordenados:
        mejor_idx = None
        mejor_final = None

        # Colocar el candidato en el camino "más ajustado" posible
        # para mantener pocos caminos y preservar domina(prev, candidato).
        for idx, camino in enumerate(caminos):
            final = camino[-1]
            if domina(final, candidato):
                if mejor_final is None or final > mejor_final:
                    mejor_final = final
                    mejor_idx = idx

        if mejor_idx is None:
            caminos.append([candidato])
        else:
            caminos[mejor_idx].append(candidato)

    return caminos

if __name__ == "__main__":
    dft_set = dict()
    result = set()
    L, C = int(sys.argv[1]), int(sys.argv[2])
    fich = open(f"candidatos_totales_{L}_{C}.txt", "w")
    for i in product(range(C), repeat=L//C):
        dft_set[PAF(i)] = i
        temp = tuple([C*(L-3)//2 -j for j in PAF(i)])
        if temp in dft_set:
            result.add(i)
            result.add(dft_set[temp])
            fich.write(f'{" ".join([str(j) for j in i])} {" ".join([str(j) for j in dft_set[temp]])}\n')
    fich.close()
    print("Number of candidates: ", len(result))
    caminos = caminos_grande(result)
    for contador, camino in enumerate(caminos):
        print(f"Camino {contador} con longitud {len(camino)}")
        with open(f"candidatos_{contador}_{L}_{C}.txt", "w") as fich:
            for p in camino:
                fich.write(" ".join([str(i) for i in p]) + "\n")
