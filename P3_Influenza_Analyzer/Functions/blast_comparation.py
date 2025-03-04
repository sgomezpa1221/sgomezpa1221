import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def verificar_secuencia_blast(query, subject_db, output_file):
    """
    Ejecuta BLASTN con una consulta y una base de datos específicas y filtra los resultados con pident > 70%.

    Parámetros:
    - query: Ruta al archivo FASTA de consulta.
    - subject_db: Ruta a la base de datos BLAST (debe estar en formato BLAST).
    - output_file: Ruta al archivo donde se guardarán los resultados.
    """
    # Formato de salida
    output_format = '"6 qseqid sseqid pident qstart qend length sstart send evalue"'

    # Construcción del comando BLAST
    command_blast = (
        f'blastn -query "{query}" '
        f'-db "{subject_db}" '
        f'-outfmt {output_format} '
        f'> "{output_file}"'
    )

    # Ejecutar BLAST
    print("Ejecutando BLAST...")
    os.system(command_blast)
    print(command_blast)
    print("BLASTN finished")
#
#     # Leer los resultados en un DataFrame
#     columns = ["qseqid", "sseqid", "pident", "qstart", "qend", "length", "sstart", "send", "evalue"]
#     df = pd.read_csv(output_file, sep='\t', names=columns)
#
#     # Filtrar por porcentaje de identidad superior al 70%
#     df_filtered = df[df['pident'] > 70]
#
#     # Mostrar los resultados filtrados
#     print("Resultados filtrados:")
#     print(df_filtered)
#
#     # Crear matriz de correlación de los genes más similares
#     correlation_matrix = df_filtered.pivot(index='qseqid', columns='sseqid', values='pident')
#
#     # Generar el heatmap con Seaborn
#     plt.figure(figsize=(10, 8))
#     sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', linewidths=0.5)
#     plt.title("Matriz de Correlación de Porcentaje de Identidad (>70%)")
#     plt.show()
#
#
# Parámetros de entrada
query_file = "C:\\Users\\lauri\\Documents\\Entrega_03_FDP\\Alcaravan_H_H3N3_AA.fasta"
database = "C:\\Users\\lauri\\Documents\\Entrega_03_FDP\\alphainfluenza_H_y_N_Genbank_data_base.fasta"
output = "Resultado_verificacion_secuencia_BLAST.txt"

# Ejecutar la función
verificar_secuencia_blast(query_file, database, output)
