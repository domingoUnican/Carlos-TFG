import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))
ax.axis('off')

# Nodos: (x, y, texto, color)
nodos = {
    'Raíz': (4, 4, 'Inicio', 'lightblue'),
    'C0_0': (2, 3, 'C0 = -1', 'lightgreen'), 'C0_1': (6, 3, 'C0 = 1', 'lightgreen'),
    'C1_0': (1, 2, 'C1 = -1', 'lightgreen'), 'C1_1': (3, 2, 'C1 = 1', 'salmon'), # Poda
    'C1_2': (5, 2, 'C1 = -1', 'lightgreen'), 'C1_3': (7, 2, 'C1 = 1', 'lightgreen'),
    'C3_0': (0.5, 1, 'C3 = -1', 'lightgreen'), 'C3_1': (1.5, 1, 'C3 = 1', 'lightgreen'),
    'C3_4': (4.5, 1, 'C3 = -1', 'lightgreen'), 'C3_5': (5.5, 1, 'C3 = 1', 'lightgreen'),
    'C3_6': (6.5, 1, 'C3 = -1', 'salmon'),     'C3_7': (7.5, 1, 'C3 = 1', 'salmon') # Poda
}

aristas = [
    ('Raíz', 'C0_0', 'black'), ('Raíz', 'C0_1', 'black'),
    ('C0_0', 'C1_0', 'black'), ('C0_0', 'C1_1', 'red'),
    ('C0_1', 'C1_2', 'black'), ('C0_1', 'C1_3', 'black'),
    ('C1_0', 'C3_0', 'black'), ('C1_0', 'C3_1', 'black'),
    ('C1_2', 'C3_4', 'black'), ('C1_2', 'C3_5', 'black'),
    ('C1_3', 'C3_6', 'red'),   ('C1_3', 'C3_7', 'red')
]

for origen, destino, color in aristas:
    x_vals = [nodos[origen][0], nodos[destino][0]]
    y_vals = [nodos[origen][1], nodos[destino][1]]
    ax.plot(x_vals, y_vals, color=color, linewidth=2, zorder=1)

for k, (x, y, texto, color) in nodos.items():
    ax.scatter(x, y, s=2000, color=color, edgecolor='black', zorder=2)
    ax.text(x, y, texto, ha='center', va='center', fontsize=9, fontweight='bold')

ax.text(3, 1.8, 'PODA POR PSD', color='red', fontweight='bold', ha='center')
ax.text(7, 0.8, 'PODA', color='red', fontweight='bold', ha='center')

plt.tight_layout()
plt.savefig('dfs_tree_l7.png', dpi=300)