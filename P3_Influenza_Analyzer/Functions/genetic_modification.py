def insertar_segmento(secuencia, nuevo_segmento, inicio, fin):
    """
    Inserta un nuevo segmento en una región específica de la secuencia,
    eliminando el fragmento original dentro del rango indicado.

    Parámetros:
    secuencia (str): Secuencia de nucleótidos original.
    nuevo_segmento (str): Nuevo segmento a insertar.
    inicio (int): Posición de inicio del segmento a reemplazar (0 basado).
    fin (int): Posición final del segmento a reemplazar (exclusivo).

    Retorna:
    str: Secuencia modificada.
    """
    if inicio < 0 or fin > len(secuencia) or inicio >= fin:
        raise ValueError("Las posiciones de inicio y fin no son válidas.")

    return secuencia[:inicio] + nuevo_segmento + secuencia[fin:]
