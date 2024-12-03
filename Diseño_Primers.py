"""
This program was developed with the intention of generating a table in .txt file format
with basic parameters for the design of a conventional PCR from .fna data obtained from GenBank.
"""

from tabulate import tabulate


def complementary_sequence (sequence):
    """ This function returns the complementary DNA strand to the sequence.

    The return sequence is created by replacing the nitrogenous base 
    on the template strand with its respective complementary base, 
    according to the pairing called by Watson and Crick
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
    """ This function returns the value of the melting temperature of the primer.
    
    Count the amount of A,T,G and C, in the first and then operate 
    with the formula 2*(A+T)*4*(G+C).
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

with open("D:\Ejercicios_FdeP_CIencias_Biológicas\Bsubtilis_CDS.fna") as file:
    # Replace the above line of code with the file to be used in .fna format

    """
    This code block:
    1. Stores the name of the gene or locus in a list.
    2. Assembles the respective gene sequence and saves it in a list. 
    """

    for line in file:
        if line.startswith(">"):
            list_genes_sequences.append(gene_sequence)
            line = line.split() # Separates the line containing the information for each gene in the .fna file.
            gene_or_locus_name = line[1].split("=")[1][: -1] # Extracts the gene or locus name.
            list_gene_names.append(gene_or_locus_name)
            gene_sequence = ""
        else:
            gene_sequence = gene_sequence + line
            gene_sequence = gene_sequence.rstrip()

list_genes_sequences.append(gene_sequence) # Add the sequence of the last gene at the end of the cycle.
list_genes_sequences.pop(0) # Deletes the first value in the list that is an empty space

""" 
Checkpoint:
The number of items in the gene list must be equal 
to the number of sequences in their respective list
"""

number_gene_names = len(list_gene_names)
number_gene_sequences = len(list_genes_sequences)
# print(number_gene_names)
# print(number_gene_sequences)

control_of_cycle = 0 # Reaching a certain value interrupts the next cycle.
size_primer = 22 # This is the most important variable in the programme and relates all of the following conditions.
main_table = [] # Put the data in the form of a table [Gene name, Primer f and r_ Tm].

# This cycle creates the rows and columns for a table.

for gene_name in list_gene_names:
    position_gene_sequences = list_genes_sequences[control_of_cycle]  # Scroll through the genes by position in the list
    complete_complementary_sequence = (complementary_sequence(position_gene_sequences))
    #  The above line goes through the sequences by position in the list
    number_of_nucleotides = f"{len(complete_complementary_sequence)} pb"
    primer_f = f"5´- {complete_complementary_sequence [0: size_primer]} -3´"
    # The above line takes the first 22 nucleotides of the complementary sequence in the 3‘- 5’ direction.
    primer_r = f"3´- {complete_complementary_sequence [-size_primer:]} -5´"
    # The above line takes the first 22 nucleotides of the complementary sequence in the 5‘- 3’ direction.
    Tm_primer_f = melting_temperature(primer_f)
    Tm_primer_r = melting_temperature(primer_r)
    Ta_primer_f = Tm_primer_f - 4 # Calculate annealign temperature of the first forward
    Ta_primer_r = Tm_primer_r - 4 # Calculate the annealign temperature of the first reverse
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

#print(total_table) # Table preview in the Python console
"""
This code block:
1. Create a document in which the results are to be written.
2. Generates the final table as a .txt file
3. Show in the Python console the end of the process.
"""
name_output_file = "Results" # Enter in this variable the name of your file, specifying the organism

file1 = open (f"{name_output_file}.txt", "w")
file1.write(total_table)
file1.close()

print(f" The document is already created in your file folder with the following name: "
      f"{name_output_file}.txt") # Allows to observe the termination of the programme execution in the console.


