import numpy as np
import matplotlib.pyplot as plt

str_A = "010100011110110110101101011110011110101000110001001110010001011100000110010"
str_B = "000100010100111111010101111001000001011001000100101111100011001001111101110"

A = np.array([1 if x == '0' else -1 for x in str_A])
B = np.array([1 if x == '0' else -1 for x in str_B])

# 1. Espectro continuo (interpolado) para ver el caos de fondo
N_plot = 750
fft_A_cont = np.fft.fft(A, n=N_plot)
fft_B_cont = np.fft.fft(B, n=N_plot)
psd_A_cont = np.abs(fft_A_cont)**2
psd_B_cont = np.abs(fft_B_cont)**2
x_cont = np.linspace(0, 75, N_plot, endpoint=False)

# 2. Espectro discreto (Matemática pura de Legendre, N=75)
fft_A_disc = np.fft.fft(A)
fft_B_disc = np.fft.fft(B)
psd_sum_disc = np.abs(fft_A_disc)**2 + np.abs(fft_B_disc)**2
x_disc = np.arange(75)

plt.figure(figsize=(10, 6))

# Pintamos las ondas continuas
plt.plot(x_cont, psd_A_cont, label='Vector A (Interpolado)', color='#1f77b4', alpha=0.7, linewidth=1)
plt.plot(x_cont, psd_B_cont, label='Vector B (Interpolado)', color='#ff7f0e', alpha=0.7, linewidth=1)

# Pintamos la SUMA DISCRETA (el teorema real)
# Omitimos el k=0 porque por definición da 2, no 152
plt.scatter(x_disc[1:], psd_sum_disc[1:], color='#2ca02c', s=40, zorder=5, 
            label='Suma Discreta Exacta (152)', marker='s')
plt.axhline(y=152, color='#2ca02c', linestyle='-', linewidth=2, zorder=4)

plt.title('Densidad Espectral de Potencia (PSD) de Pares de Legendre (N=75)', fontsize=14, pad=15)
plt.xlabel('Frecuencia (k)', fontsize=12)
plt.ylabel('Amplitud Espectral (Energía)', fontsize=12)

plt.ylim(0, 160)
plt.xlim(0.5, 74.5)

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='lower right', fontsize=11)

plt.tight_layout()
plt.savefig('real_fluctuaciones_psd_image.png', dpi=300)