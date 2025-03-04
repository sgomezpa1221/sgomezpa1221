def mostrar_ayuda():
    """Imprime la documentación del paquete de análisis de influenza en la terminal."""
    documentacion = """
# Documentación del Paquete de Análisis de Influenza

## Descripción General
Este paquete permite el análisis de secuencias genéticas del virus de la influenza, incluyendo cálculo de contenido GC y TA, aplicación de mutaciones, traducción a proteínas, creación de árboles filogenéticos y análisis con BLAST.

## Módulos y Funciones

### 1. Conteo de Contenidos GC y TA
**Módulo:** `conteo_gc_ta.py`
- `gc_content(sequence)`: Calcula el porcentaje de nucleótidos GC en una secuencia de ADN.
- `ta_content(sequence)`: Calcula el porcentaje de nucleótidos TA en una secuencia de ADN.

### 2. Aplicación de Mutaciones
**Módulo:** `mutations_insertion.py`
- `aplicar_mutaciones(sequence, tipo="transicion")`: Aplica mutaciones a la secuencia de ADN.

### 3. Inserción de Segmento
**Módulo:** `genetic_modification.py`
- `insertar_segmento(sequence, segmento, inicio, fin)`: Inserta un segmento específico en la secuencia.

### 4. Traducción de Secuencias
**Módulo:** `new_proteins.py`
- `traducir_secuencia(sequence)`: Traduce una secuencia de ADN a proteína.

### 5. Visualización de Datos
**Módulos:** `plot_gc_ta_content.py`, `proteins_plot.py`
- `graficar_scatter_contenidos(nombres, gc_values, ta_values)`: Grafica los contenidos GC y TA.
- `graficar_mapa_identidad(secuencias)`: Muestra el mapa de identidad de proteínas.

### 6. Análisis Filogenético
**Módulo:** `phylogenetic_creator.py`
- `phylogenetic_tree(query_file, iq_tree_app)`: Genera un árbol filogenético.

### 7. Análisis BLAST
**Módulos:** `Blastn_analyze.py`, `blast_database.py`
- `blast_verification(query_file, database, output)`: Realiza una búsqueda BLAST.
- `generar_db_blast(database)`: Crea una base de datos BLAST.

## Uso Recomendado
1. Cargar una secuencia de ADN.
2. Calcular contenido GC y TA.
3. Aplicar mutaciones.
4. Traducir la secuencia a proteína.
5. Comparar con BLAST.
6. Generar un árbol filogenético.

## Contacto
Para dudas o mejoras, contactar al equipo de desarrollo.
    """
    print(documentacion)

# Asignación a la variable ayuda
ayuda = mostrar_ayuda

# Ejemplo de uso
if __name__ == "__main__":
	ayuda()  # Llamada a la variable para mostrar la documentación
