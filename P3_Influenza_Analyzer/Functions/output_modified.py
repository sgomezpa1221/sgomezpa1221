from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

def guardar_fasta(nombre_archivo, nombres, secuencias):
    """
    Guarda las secuencias en un archivo FASTA compatible con BLAST.

    Parámetros:
    - nombre_archivo (str): Nombre del archivo FASTA de salida.
    - nombres (list): Lista con los nombres de cada secuencia.
    - secuencias (list): Lista con las secuencias en formato string.
    """
    registros = []
    for nombre, secuencia in zip(nombres, secuencias):
        record = SeqRecord(Seq(secuencia), id=nombre, description="|modified sequence|")
        registros.append(record)

    # Guardar el archivo en formato FASTA
    with open(nombre_archivo, "w") as output_handle:
        SeqIO.write(registros, output_handle, "fasta")

    print(f"✅ Archivo {nombre_archivo} guardado con éxito.")
