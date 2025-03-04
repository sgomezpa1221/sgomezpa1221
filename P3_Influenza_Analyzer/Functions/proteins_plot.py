import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def graficar_mapa_identidad(proteinas, nombres):
    """
    Genera una matriz de identidad comparando la similitud de secuencias de proteínas.

    Parámetros:
    - proteinas (list): Lista de secuencias de proteínas.
    - nombres (list): Lista de nombres de las proteínas.
    """
    num_proteinas = len(proteinas)
    matriz_identidad = np.zeros((num_proteinas, num_proteinas))

    for i in range(num_proteinas):
        for j in range(num_proteinas):
            identidad = sum(a == b for a, b in zip(proteinas[i], proteinas[j])) / min(len(proteinas[i]), len(proteinas[j]))
            matriz_identidad[i, j] = identidad * 100  # Convertir a porcentaje

    # Ajustar el tamaño del gráfico
    plt.figure(figsize=(7, 5))

    # Crear el mapa de calor con nueva estética
    sns.heatmap(
        matriz_identidad,
        annot=True, fmt=".1f", cmap="mako", linewidths=0.5,
        cbar_kws={"shrink": 0.75, "label": "Porcentaje de Identidad"},
        xticklabels=nombres,
        yticklabels=nombres
    )

    # Mejorar etiquetas y título
    plt.title("🔬 Mapa de Identidad de Proteínas (%)", fontsize=14, fontweight='bold', pad=12)
    plt.xlabel("Proteína Comparada", fontsize=12)
    plt.ylabel("Proteína Referencia", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(fontsize=10)

    plt.show()
