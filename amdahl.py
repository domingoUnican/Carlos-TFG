import matplotlib.pyplot as plt
import numpy as np

n_cores = np.arange(1, 129)
P = 0.9999
speedup = 1 / ((1 - P) + (P / n_cores))
speedup_ideal = n_cores

plt.figure(figsize=(8, 5))
plt.plot(n_cores, speedup, label='Speedup Real (P=0.9999)', color='blue', linewidth=2.5)
plt.plot(n_cores, speedup_ideal, label='Escalabilidad Ideal (Lineal)', color='red', linestyle='--', alpha=0.7)

plt.title('Escalabilidad Teórica del Algoritmo (Ley de Amdahl)')
plt.xlabel('Número de Núcleos ($N$)')
plt.ylabel('Aceleración ($Speedup$)')
plt.xlim(1, 128)
plt.ylim(1, 128)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.savefig('amdahl_plot.png', dpi=300)