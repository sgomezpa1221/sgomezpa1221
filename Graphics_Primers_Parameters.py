import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate
from matplotlib.backends.backend_pdf import PdfPages

"""
The development of this code aims to establish a starting point for the development of primers and consequently,
aspects related to the process of a conventional PCR for data in .fna format obtained from GenBank.This program allows 
to obtain a basic data report from graphs that intuitively relate the behaviour of the melting and annealing 
temperatures for the first forward and reverse that are obtained considering the first and last 22 nucleotides in the 
sequence of each gene. 
"""


def complementary_sequence (sequence):
    """ This function returns the complementary DNA strand to the sequence. The return sequence is created by replacing
     the nitrogenous base on the template strand with its respective complementary base, according to the pairing named
     by Watson and Crick.
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
    """This function returns the melting temperature value of the primer. It counts the amount of A, T, G and C in the
    primer and then operates on the formula 2*(A+T)*4*(G+C).
    """
    primer = Temp_melting
    count_A = primer.count("A")
    count_T = primer.count("T")
    count_G = primer.count("G")
    count_C = primer.count("C")
    Temp_melting = 2*(count_A + count_T) + 4*(count_G + count_C)
    return Temp_melting


gene_sequence = ""
list_gene_names = [] # Stores the name of all genes coming from the .fna file
list_genes_sequences = [] # Stores the frequency of genes originating from the arhvian .fna

with open("D:\Ejercicios_FdeP_CIencias_Biológicas\Bsubtilis_CDS.fna.txt") as file:
    # Replace the above line of code with the location of the file to be used in .fna format

    """ This code block:
    1. It stores the name of the gene or locus in a list.
    2. Assembles the respective sequence of the gene and stores it in a list. 
    """

    for line in file:
        if line.startswith(">"):
            list_genes_sequences.append(gene_sequence)
            line = line.split() # Separates the line containing the information for each gene in the .fna file
            gene_or_locus_name = line[1].split("=")[1][: -1] # Extract the name of the gene or locus
            list_gene_names.append(gene_or_locus_name)
            gene_sequence = ""
        else:
            gene_sequence = gene_sequence + line
            gene_sequence = gene_sequence.rstrip()


list_genes_sequences.append(gene_sequence) # Add the sequence of the last gene at the end of the cycle.
list_genes_sequences.pop(0) # Deletes the first value in the list that is an empty space

""" 
Checkpoint:
The number of elements in the gene list must be equal to the number of sequences in the respective gene list. 
be equal to the number of sequences in their respective list.
"""

number_gene_names = len(list_gene_names)
number_gene_sequences = len(list_genes_sequences)

control_of_cycle = 0 # It interrupts the next cycle and also allows to scroll through the subsequent lists by position.
size_primer = 22 # This is the most important variable in the programme and relates all of the following conditions

""" 
The following code block creates the row and column vectors for a DataFrame obtained from the sequences of each gene in
the .fna file.
"""

list_number_of_nucleotides = []
list_of_primers_F = []
list_of_primers_R = []
list_of_Tm_primers_F = []
list_of_Tm_primers_R = []
list_of_Ta_primers_F = []
list_of_Ta_primers_R = []

for gene_name in list_gene_names:
    position_gene_sequences = list_genes_sequences[control_of_cycle]  # Scroll through the genes by position in the list
    complete_complementary_sequence = (complementary_sequence(position_gene_sequences))
    # The above line goes through the sequences by position in the list
    list_number_of_nucleotides.append(len(complete_complementary_sequence))

    primer_f = complete_complementary_sequence[0: size_primer]
    # Take the first 22 nucleotides of the complementary sequence in the 3‘- 5’ direction.
    list_of_primers_F.append(primer_f)
    primer_r = complete_complementary_sequence[-size_primer:]
    # Take the first 22 nucleotides of the complementary sequence in the 5‘- 3’ direction.
    list_of_primers_R.append(primer_r)

    # Metlting temperature
    Tm_primer_f = melting_temperature(primer_f)
    list_of_Tm_primers_F.append(Tm_primer_f)
    Tm_primer_r = melting_temperature(primer_r)
    list_of_Tm_primers_R.append(Tm_primer_r)

    # Annealing temperature
    list_of_Ta_primers_F.append(Tm_primer_f - 4)
    list_of_Ta_primers_R.append(Tm_primer_r - 4)

    control_of_cycle += 1

df = pd.DataFrame({"Gene or locus name":list_gene_names,
                   "Number of nucleotides":list_number_of_nucleotides,
                   "Primer F":list_of_primers_F,
                   "Primer R":list_of_primers_R,
                   "Tm Primer F":list_of_Tm_primers_F,
                   "Tm Primer R":list_of_Tm_primers_R,
                   "Ta Primer F":list_of_Ta_primers_F,
                   "Ta Primer R":list_of_Ta_primers_R}
                  )
""" 
The following line of code seeks to display the minimum and maximum values that are considered relevant in the DataFrame
and therefore in the genome of the organism from which the .fna file comes from.
"""
parameters = []
row_max_number_nucleotides = df.loc[df["Number of nucleotides"].idxmax()]
# The previous line Identifies the row of the gene with the highest number of nucleotides
parameters.append(["Higher number of nucleotides", row_max_number_nucleotides["Gene or locus name"],
                   row_max_number_nucleotides["Number of nucleotides"]])
# The previous line adds the most relevant information from the line with the highest number of nucleotides to the list.
row_min_number_nucleotides = df.loc[df["Number of nucleotides"].idxmin()]
# The previous line Identifies the row of the gene with the lowest number of nucleotides
parameters.append(["Lower number of nucleotides", row_min_number_nucleotides["Gene or locus name"],
                   row_min_number_nucleotides["Number of nucleotides"]])
# The previous line adds the most relevant information from the line with the lowest number of nucleotides to the list.
row_max_Tm_primer_F = df.loc[df["Tm Primer F"].idxmax()]
parameters.append(["Máx Tm Primer F", row_max_Tm_primer_F["Gene or locus name"], row_max_Tm_primer_F["Tm Primer F"]])
row_max_Tm_primer_R = df.loc[df["Tm Primer R"].idxmax()]
parameters.append(["Máx Tm Primer R", row_max_Tm_primer_R["Gene or locus name"], row_max_Tm_primer_R["Tm Primer R"]])
row_min_Tm_primer_F = df.loc[df["Tm Primer F"].idxmin()]
parameters.append(["Mín Tm Primer F", row_min_Tm_primer_F["Gene or locus name"], row_min_Tm_primer_F["Tm Primer F"]])
row_min_Tm_primer_R = df.loc[df["Tm Primer R"].idxmin()]
parameters.append(["Mín Tm Primer R", row_min_Tm_primer_R["Gene or locus name"], row_min_Tm_primer_R["Tm Primer R"]])

"""
This code block contains two ways of generating the data report in a PDF file, in the first one each graph is displayed 
independently, while in the second one, all graphs are presented on a single page.
"""

# Generate PDF with all larger charts on multiple pages
with PdfPages('output_graphics.pdf') as pdf:
# Histogram
    sns.histplot(df['Number of nucleotides'],
                 color = '#006d77',
                 kde = True, #Shows the distribution line
                 bins = 'fd' # Define the intervals according to the Freedman Diaconis rule.
                 )
    plt.title("DNA Sequence Length Distribution")
    plt.xlabel("Number of Nucleotides [log]")
    plt.xscale('log')  # Escala logarítmica en el eje X
    plt.ylabel("Frequency")
    # plt.show()
    pdf.savefig()
    plt.close()

# Boxplot
    temperature_columns = ["Tm Primer F", "Tm Primer R", "Ta Primer F", "Ta Primer R"]

    df_long = df.melt(value_vars = temperature_columns, var_name = "Temperature", value_name = "Values")

    plt.figure(figsize = (10, 6))
    sns.boxplot(x = "Temperature",
                y = "Values",
                hue = "Temperature",
                palette = "viridis",
                data = df_long,
                legend = False
                )
    sns.stripplot(df_long,
                x ="Temperature",
                y = "Values",
                size = 1,
                color = "0.2"
                )

    plt.title("Distribution of Temperatures", fontsize=16)
    plt.ylabel("Values (°C)", fontsize=14)
    # plt.show()
    pdf.savefig()
    plt.close()

# Heatmap
    columns_of_interest = df[["Number of nucleotides", "Tm Primer F", "Tm Primer R", "Ta Primer F", "Ta Primer R"]]
    correlation_matrix = columns_of_interest.corr()

    plt.figure(figsize=(6, 4))
    sns.heatmap(correlation_matrix,
                annot=True,
                fmt=".2f",
                cmap="coolwarm",
                linewidths=0.5,
                square=True,
                cbar_kws={"label": "Correlation Coefficient"}
                )

    plt.title("Correlation Matrix: Gene Length and Temperatures", fontsize=12)
    # plt.show()
    pdf.savefig()
    plt.close()

# Scatterplot
    sns.scatterplot(y=df['Number of nucleotides'], x=df['Tm Primer F'], label="Tm primer F", alpha=0.5)
    sns.scatterplot(y=df['Number of nucleotides'], x=df['Tm Primer R'], label="Tm primer R", alpha=0.5)
    sns.scatterplot(y=df['Number of nucleotides'], x=df['Ta Primer R'], label="Ta primer F", alpha=0.5)
    sns.scatterplot(y=df['Number of nucleotides'], x=df['Ta Primer R'], label="Ta primer R", alpha=0.1)
    plt.title("Relationship between Number of Nucleotides and Temperature in PCR")
    plt.ylabel("Number of Nucleotides")
    plt.xlabel("Temperature")
    plt.legend()
    # plt.show()
    pdf.savefig()
    plt.close()

# Lineplot
    plt.figure(figsize=(10, 6))
    # Tm Primer F graph
    sns.lineplot(y='Number of nucleotides',
                 x='Tm Primer F',
                 data=df,
                 marker='o',
                 label='Tm Primer F',
                 color='b'
                 )
    # Tm Primer R graph
    sns.lineplot(y='Number of nucleotides',
                 x='Tm Primer R',
                 data=df,
                 marker='o',
                 label='Tm Primer R',
                 color='r'
                 )

    plt.title('Number of Nucleotides vs Temperatures of primers F and R', fontsize=16)
    plt.xlabel('Number', fontsize=12)
    plt.ylabel('Temperature (°C)', fontsize=12)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    # plt.show()
    pdf.savefig()
    plt.close()

# Relevant parameter report table
    plt.figure(figsize=(10, 4))
    plt.axis('off')
    table_data = [["Parameter", "Name of gene or locus", "Value"]] + parameters
    table = plt.table(cellText=table_data,
                      colLabels=None,
                      loc='center',
                      cellLoc='center',
                      colColours=["#f2f2f2", "#f9f9f9", "#f2f2f2"])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

    plt.title("Summary of Key Parameters", fontsize=16, y=1.08)
    pdf.savefig()
    plt.close()

# Generate PDF with all small graphics on a single page
with PdfPages('single_page_report.pdf') as pdf:
    # Create a figure with subplots (organised in 3 rows and 3 columns)
    fig, axes = plt.subplots(3, 3, figsize=(18, 12))
    fig.suptitle("Summary of Graphs and Parameters", fontsize=16)

    # Histograma
    sns.histplot(df['Number of nucleotides'], kde=True, bins='fd', ax=axes[0, 0], color='skyblue')
    axes[0, 0].set_title("Sequence Length Distribution")
    axes[0, 0].set_xlabel("Number of Nucleotides")
    axes[0, 0].set_ylabel("Frecuency")

    # Boxplot
    temperature_columns = ["Tm Primer F", "Tm Primer R", "Ta Primer F", "Ta Primer R"]
    df_long = df.melt(value_vars=temperature_columns, var_name="Temperature", value_name="Values")
    sns.boxplot(x="Temperature", y="Values", data=df_long, ax=axes[0, 1], palette="viridis")
    axes[0, 1].set_title("Distribution of Temperatures")
    axes[0, 1].set_ylabel("Values (°C)")

    # Heatmap
    columns_of_interest = ["Number of nucleotides", "Tm Primer F", "Tm Primer R", "Ta Primer F", "Ta Primer R"]
    correlation_matrix = df[columns_of_interest].corr()
    sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap="coolwarm", ax=axes[0, 2])
    axes[0, 2].set_title("Correlation matrix")

    # Scatterplot
    sns.scatterplot(x=df['Tm Primer F'], y=df['Number of nucleotides'], ax=axes[1, 0], alpha=0.5, label="Tm Primer F")
    sns.scatterplot(x=df['Tm Primer R'], y=df['Number of nucleotides'], ax=axes[1, 0], alpha=0.5, label="Tm Primer R")
    sns.scatterplot(x=df['Ta Primer F'], y=df['Number of nucleotides'], ax=axes[1, 0], alpha=0.5, label="Ta Primer F")
    sns.scatterplot(x=df['Ta Primer R'], y=df['Number of nucleotides'], ax=axes[1, 0], alpha=0.1, label="Ta Primer R")
    axes[1, 0].set_title("Relationship between Length and Tm")
    axes[1, 0].set_xlabel("Temperature (°C)")
    axes[1, 0].set_ylabel("Number of Nucleotides")
    axes[1, 0].legend()

    # Lineplot
    sns.lineplot(x='Tm Primer F', y='Number of nucleotides', data=df, ax=axes[1, 1], label="Tm Primer F")
    sns.lineplot(x='Tm Primer R', y='Number of nucleotides', data=df, ax=axes[1, 1], label="Tm Primer R")
    axes[1, 1].set_title("Number of Nucleotides vs Tm")
    axes[1, 1].set_xlabel("Melting Temperature (Tm)")
    axes[1, 1].set_ylabel("Number of Nucleotides")
    axes[1, 1].legend()

    # Table
    axes[2, 0].axis('off')
    table_data = [["Parameter", "Name of gene or locus", "Value"]] + parameters
    table = axes[2, 0].table(cellText=table_data, colLabels=None, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)
    axes[2, 0].set_title("Summary of keys parameters")

    axes[1, 2].axis('off')
    axes[2, 1].axis('off')
    axes[2, 2].axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig(fig)
    plt.close()
