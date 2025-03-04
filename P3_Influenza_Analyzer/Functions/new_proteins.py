from Bio.Seq import Seq
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def traducir_secuencia(nucleotidos):
    """
    Traduce una secuencia de nucleótidos a aminoácidos usando el código genético estándar.

    Parámetro:
    - nucleotidos: Secuencia de ADN en formato string.

    Retorna:
    - Secuencia de aminoácidos en formato string.
    """
    secuencia_adn = Seq(nucleotidos)  # Convertir en objeto Seq de Biopython
    proteina = secuencia_adn.translate(to_stop=True)  # Traducir hasta el primer codón de stop
    return str(proteina)


def graficar_heatmap(proteinas, nombres):
    """
    Genera un heatmap comparando las secuencias de proteínas tras las mutaciones.

    Parámetros:
    - proteinas: Lista de secuencias de aminoácidos.
    - nombres: Lista de nombres de cada condición.

    Retorna:
    - Muestra un heatmap con las diferencias en las secuencias.
    """
    # Determinar la longitud máxima de las secuencias
    longitud_max = max(len(p) for p in proteinas)

    # Crear matriz numérica para el heatmap
    matriz = np.full((len(proteinas), longitud_max), np.nan)  # Llenar con NaN para posiciones sin datos

    for i, proteina in enumerate(proteinas):
        for j, aa in enumerate(proteina):
            matriz[i, j] = ord(aa)  # Convertir el aminoácido a valor numérico

    # Crear heatmap
    plt.figure(figsize=(15, 4))
    ax = sns.heatmap(matriz, cmap="coolwarm", linewidths=0.5, xticklabels=10,
                     yticklabels=nombres, cbar_kws={'label': 'Código ASCII de AA'})

    # Personalización
    plt.title("Comparación de Secuencias de Proteínas tras Mutaciones")
    plt.xlabel("Posición en la Secuencia de Aminoácidos")
    plt.ylabel("Condición")
    plt.show()
