import random

def aplicar_transiciones(sequence):
    """Aplica mutaciones de transición (A↔G, C↔T) en una región aleatoria, asegurando un cambio significativo."""
    seq_list = list(sequence)
    n = len(seq_list)

    transiciones = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}

    # Seleccionar un segmento más grande (100-300 nucleótidos)
     #Seleccionar un punto de inicio en cualquier parte de la secuencia excepto la última posición
    start = random.randint(0, n - 1)

    # Determinar un tamaño aleatorio para el segmento (100-300 nucleótidos)
    segment_size = random.randint(100, 500)

    # Asegurar que 'end' no supere el tamaño total de la secuencia
    end = min(start + segment_size, n - 1)

    # Determinar qué porcentaje de nucleótidos mutar (50% - 100%)
    mutar_ratio = random.uniform(0.5, 1.0)

    # Aplicar transiciones
    indices = random.sample(range(start, end), int((end - start) * mutar_ratio))
    for i in indices:
        seq_list[i] = transiciones[seq_list[i]]

    return "".join(seq_list), start, end

def aplicar_translocaciones(sequence):
    """Intercambia un segmento aleatorio con otro distante en la secuencia."""
    seq_list = list(sequence)
    n = len(seq_list)

    # Seleccionar un segmento más grande (100-300 nucleótidos)
    start1 = random.randint(0, n - 500)
    end1 = start1 + random.randint(100, 500)

    # Asegurar que la segunda región esté al menos a 500 nucleótidos de distancia
    start2 = random.randint(max(0, start1 - 1000), min(n - 300, start1 + 1000))
    while abs(start1 - start2) < 500:
        start2 = random.randint(max(0, start1 - 1000), min(n - 300, start1 + 1000))

    end2 = start2 + (end1 - start1)

    # Intercambiar segmentos
    seq_list[start1:end1], seq_list[start2:end2] = seq_list[start2:end2], seq_list[start1:end1]

    return "".join(seq_list), start1, end1

def aplicar_mutaciones(sequence, tipo="ambas"):
    """Aplica transiciones, translocaciones o ambas según el tipo."""
    if tipo == "transicion":
        return aplicar_transiciones(sequence)
    elif tipo == "translocacion":
        return aplicar_translocaciones(sequence)
    elif tipo == "ambas":
        seq_mut, start1, end1 = aplicar_transiciones(sequence)
        seq_mut, start2, end2 = aplicar_translocaciones(seq_mut)
        return seq_mut, (start1, end1), (start2, end2)
    else:
        raise ValueError("Tipo de mutación no válido: use 'transicion', 'translocacion' o 'ambas'")
