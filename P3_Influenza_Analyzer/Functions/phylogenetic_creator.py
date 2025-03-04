import os

def phylogenetic_tree(query, iq_tree):
    """
     Ejecuta IQ-TREE con un archivo que contiene secuencias de genes.
     Parámetros:
     - query: Ruta al archivo PHYLIP (.phy) de entrada.
     - iq_tree: Ruta al ejecutable de IQ-TREE.
    Retorna:
     - None

"""
# Construcción del comando
    command_iq_tree = f'{iq_tree} -s "{query}"'

# Ejecutar IQ-TREE
    print("Ejecutando IQ-TREE...")
    print(command_iq_tree)  # Para depuración
    os.system(command_iq_tree)
    print("IQ-TREE finished")
