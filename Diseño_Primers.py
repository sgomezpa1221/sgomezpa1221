"""
Este programa se desarollo con la intención de generar una tabla en formato de archivo.txt
con parametros básicos para el diseño de una PCR convencional a partir de datos .fna
obtenidos del GenBank.
"""

from tabulate import tabulate


def complementary_sequence (sequence):
    """ Esta función retorna la cadena complementaria de ADN a la secuencia.

    La secuencia de retorno se crea reemplazando la base nitrogenada en
    la hebra molde por su respectiva base complementaria, según el apareamiento
    denominado por Watson y Crick
    """
    complementary_sequence = ""
    for nucleotide in sequence:
        if nucleotide == "T":
            complementary_sequence += "A"
        elif nucleotide == "A":
            complementary_sequence += "T"
        elif nucleotide == "G":
            complementary_sequence += "C"
        elif nucleotide == "C":
            complementary_sequence += "G"
        else:
            pass
    return complementary_sequence


def melting_temperature (Temp_melting):
    """ Esta función retorna el valor de la temperatura de melting del primer.

    Cuenta la cantidad de A,T,G y C, en el primer y despúes opera con la formula
    2*(A+T)*4*(G+C)
    """
    primer = Temp_melting
    count_A = primer.count("A")
    count_T = primer.count("T")
    count_G = primer.count("G")
    count_C = primer.count("C")
    Temp_melting = 2*(count_A + count_T) + 4*(count_G + count_C)
    return Temp_melting


gene_sequence = ""
list_gene_names = []
list_genes_sequences = []

with open("D:\Ejercicios_FdeP_CIencias_Biológicas\Bsubtilis_CDS.fna.txt") as file:
    # Reemplazar la anterior línea de código con el archivo a utilizar en formato .fna

    """
    Este bloque de código:
    1. Almacena el nombre del gen o locus en una lista.
    2. Ensambla la secuencia respectiva del gen y guarda en una lista 
    """

    for line in file:
        if line.startswith(">"):
            list_genes_sequences.append(gene_sequence)
            line = line.split() # Separa la línea que contiene la información de cada gen  en el archivo .fna
            gene_or_locus_name = line[1].split("=")[1][: -1] # Extra el nombre del gen o locus
            list_gene_names.append(gene_or_locus_name)
            gene_sequence = ""
        else:
            gene_sequence = gene_sequence + line
            gene_sequence = gene_sequence.rstrip()

list_genes_sequences.append(gene_sequence) # Añade la secuencia del último gen al terminar el ciclo
list_genes_sequences.pop(0) # Elimina el primer valor de la lista que es un espacio vacío

""" 
Punto de control:
El número de elementos en la lista de genes debe 
ser igual al número de secuencias en su respectiva lista
"""

number_gene_names = len(list_gene_names)
number_gene_sequences = len(list_genes_sequences)
# print(number_gene_names)
# print(number_gene_sequences)

control_of_cycle = 0 # Al alcanzar cierto valor interrumpe el siguiente ciclo
size_primer = 22 # Esta es la variable más importante del programa y relaciona todas las condiciones siguientes
main_table = [] # Pondra los datos en forma de tabla [Nombre del gen, Primer f y r_ Tm]

# Este ciclo crea las filas y columnas para una tabla

for gene_name in list_gene_names:
    position_gene_sequences = list_genes_sequences[control_of_cycle]  #Recorre los genes por posición en la lista
    complete_complementary_sequence = (complementary_sequence(position_gene_sequences))
    #  La anterior línea recorre las secuencias por posición en la lista
    number_of_nucleotides = f"{len(complete_complementary_sequence)} pb"
    primer_f = f"5´- {complete_complementary_sequence [0: size_primer]} -3´"
    # La anterior línea toma los primeros 22 nucleotidos de la secuencia complementaria en dirección 3'- 5'
    primer_r = f"3´- {complete_complementary_sequence [-size_primer:]} -5´"
    #La anterior línea toma los primeros 22 nucleotidos de la secuencia complementaria en dirección 5'- 3'
    Tm_primer_f = melting_temperature(primer_f)
    Tm_primer_r = melting_temperature(primer_r)
    Ta_primer_f = Tm_primer_f - 4 # Calcula la temperatura de annealign del primer forward
    Ta_primer_r = Tm_primer_r - 4 # Calcula la temperatura de annealign del primer reverse
    table_renglon = [gene_name,
                     number_of_nucleotides,
                     primer_f, primer_r,
                     Tm_primer_f, Tm_primer_r,
                     Ta_primer_f, Ta_primer_r
                     ]

    main_table.append(table_renglon)
    control_of_cycle += 1
    # if control_of_cycle > 2:
    #     break
#print(main_table)

total_table = tabulate(main_table, headers=["Gene or locus name",
                                            "Number of nucleotides",
                                            "Primer F", "Primer R",
                                            "Tm Primer F", "Tm Primer R",
                                            "Ta Primer F", "Ta Primer R"],
                       colalign=("center", "center",
                                 "center", "center",
                                 "center", "center",
                                 "center", "center" )
                       )

#print(total_table) # Vista previa de la tabla en la consola de Python
"""
Esta bloque de código:
1. Crea un documento en el cual se van a escribir los resultados
2. Genera la tabla final en un archivo .txt
3. Muestra en la consola de Python la finalización del proceso
"""
name_output_file = "Resultados" # Introduzca en esta variable el nombre de su archiivo, especificando el organismo

file1 = open (f"{name_output_file}.txt", "w")
file1.write(total_table)
file1.close()

print(f" El documento ya se creo en tu carpeta de archivos con el siguiente nombre: "
      f"{name_output_file}.txt") # Permite observar la terminación de la ejecución del programa en la consola


