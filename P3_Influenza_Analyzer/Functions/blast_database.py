import os
def generar_db_blast(query_file, database, output):
    """
    Ejecuta BLASTN para generar una base de datos específicas.

    Parámetros:
    - subject_db: Ruta a la base de datos BLAST (debe estar en formato BLAST).

    Retorna:
    - None
    """

    # Construcción del comando
    command_blast = (
        f'makeblastdb -in "{database}" '
        f'-dbtype nucl'
    )

    # Ejecutar BLAST
    print("Ejecutando BLAST...")
    print(command_blast)  # Para depuración
    os.system(command_blast)
    print("BLASTN finished")
