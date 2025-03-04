""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Librerias
import pandas as pd
import random
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import os

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Funciones
from Funciones.conteo_gc_ta import gc_content, ta_content
from Funciones.mutations_insertion import aplicar_mutaciones
from Funciones.genetic_modification import insertar_segmento
from Funciones.plot_gc_ta_content import graficar_scatter_contenidos
from Funciones.new_proteins import traducir_secuencia, graficar_heatmap
from Funciones.proteins_plot import graficar_mapa_identidad
from Funciones.output_modified import guardar_fasta
from Funciones.Blastn_analyze import blast_verification
from Funciones.phylogenetic_creator import phylogenetic_tree
from Funciones.blast_database import generar_db_blast
from Funciones.ayuda_influenza_analyzer import mostrar_ayuda
"""
# """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# """" Prueba Secuencia Conocida"""
#
# """
#     - Se realizo el analisis de blast y filogenia con una secuencia de hemaglutinina de el Alcaravan
#     una especie de pajaro.
#     - Estas funciones se realizan solo una vez.
# """
# """
#     1. Creacion de una database de las secuencias de las 18 tipos de Hemaglutininas conocidas
#     de influeza A
# """
# # # Parámetros de entrada
# # database = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Base_de_datos\\alphainfluenza_H_y_N_Genbank_data_base.fasta"
# # generar_db_blast(database)  # Llamar a la función
# # #
# # # """
# # """
# #     2. Creacion del archivo iqtree
# # """
# # # Parámetros de entrada
# # query_file = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Hemaglutinina.phy"
# # iq_tree_app = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\iqtree-1.6.12-Windows\\bin\\iqtree.exe"
# # phylogenetic_tree(query_file, iq_tree_app) # Llamar a la función
#
#
# """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # """
# #     - Procesamiento la secuencia conocida para realizar los analisis -
# # """
# #
# # # Ruta al archivo nucleotidos
# # Nucleotides_file = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Validación\\Hemaglutinina\\Alcaraván_H_H3N3_AA.fasta"
# # # Ruta al archivo de interes (aminoacidos)
# # AA_file = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Validación\\Hemaglutinina\\Alcaraván_H_H3N3_PRO.fasta"
# #
# # # Lectura y alineamiento de las secuencias (de Nucleotidos y de Proteinas)
# #
# # sequence = ""
# # counter = 0
# # with open(Nucleotides_file, 'r') as file:
# #     next(file)  # Saltar la primera línea
# #     for line in file:
# #         counter += 1
# #         sequence += line.strip()
# #
# # protein = ""
# # count = 0
# # with open(AA_file, 'r') as file:
# #     next(file)  # Saltar la primera línea
# #     for line in file:
# #         count += 1
# #         protein += line.strip()
# #
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # # Aplica la funcion de los contenidos de GC y TA
# # gc_ori, ta_ori = gc_content(sequence),ta_content(sequence)
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # #  Aplica transiciones
# # seq_tr, st_tr, end_tr = aplicar_mutaciones(sequence, tipo="transicion")
# # gc_tr,ta_tr = gc_content(seq_tr), ta_content(seq_tr)
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # # Aplica translocaciones
# # seq_tl, st_tl, end_tl = aplicar_mutaciones(sequence, tipo="translocacion") # Mutación
# # gc_tl, ta_tl= gc_content(seq_tl),ta_content(seq_tl)
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # # Aplica ambas mutaciones (transiciones + translocaciones)
# # seq_both, st_both, end_both = aplicar_mutaciones(sequence, tipo="ambas") # Mutación
# # gc_both, ta_both = gc_content(seq_both), ta_content(seq_both)
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # # Inserción de una secuencia especifica
# # inicio_insert = 533
# # fin_insert = 611
# # nuevo_segmento = "ATGCGTACGTTAGC"
# # seq_mod = insertar_segmento(sequence, nuevo_segmento, inicio_insert, fin_insert)
# # gc_mod, ta_mod  = gc_content(seq_mod), ta_content(seq_mod)
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # # Traducción de las secuencias
# # o_pro = traducir_secuencia(sequence)
# # tr_pro = traducir_secuencia(seq_tr)
# # tl_pro = traducir_secuencia(seq_tl)
# # both_pro = traducir_secuencia(seq_both)
# # mod_pro = traducir_secuencia(seq_mod)
# #
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # #Resultados Numéricos
# # print(f"\n🔹 Antes de la mutación:")
# # print(f"Contenido GC : {gc_ori:.2f}%")
# # print(f"Contenido TA : {ta_ori:.2f}%")
# #
# # print(f"\n🔹 Solo transiciones:")
# # print(f"Contenido GC : {gc_tr:.2f}%")
# # print(f"Contenido TA : {ta_tr:.2f}%")
# # print(f"Segmento mutado entre posiciones {st_tr} y {end_tr}")
# # # print(seq_transiciones)
# #
# # print(f"\n🔹 Solo translocaciones:")
# # print(f"Contenido GC : {gc_tl:.2f}%")
# # print(f"Contenido TA : {ta_tl:.2f}%")
# # print(f"Segmento mutado entre posiciones {st_tl} y {end_tl}")
# # # print(seq_translocaciones)
# #
# # print(f"\n🔹 Ambas mutaciones (transiciones + translocaciones):")
# # print(f"Contenido GC : {gc_both:.2f}%")
# # print(f"Contenido TA : {ta_both:.2f}%")
# # print(f"Segmento mutado entre posiciones {st_both} y {end_both}")
# # # print(seq_ambas)
# #
# # print(f"\n🔹 Modificación genética en posiciones {inicio_insert} - {fin_insert}:")
# # print(f"Contenido GC : {gc_mod:.2f}%")
# # print(f"Contenido TA : {ta_mod:.2f}%")
# # #print(seq_modificada)
# #
# # # Comparación de las secuencias de proteínas
# # print(f"\n🔹 Secuencia de aminoácidos ORIGINAL:\n{o_pro}\n")
# # print(f"🔹 Tras transiciones:\n{tr_pro}\n")
# # print(f"🔹 Tras translocaciones:\n{tl_pro}\n")
# # print(f"🔹 Tras ambas mutaciones:\n{both_pro}\n")
# # print(f"🔹 Tras modificación genética:\n{mod_pro}\n")
# #
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# # # Graficación
# # #------------------------------------- Contenidos GC y TA--------------------------------------------#
# # graficar_scatter_contenidos(["Ori", "Tr", "Tl", "Both", "M-G"],
# #                               [gc_ori, gc_tr, gc_tl, gc_both, gc_mod],
# #                               [ta_ori, ta_tr, ta_tl, ta_both, ta_mod])
# # #-------------------------------------- Nuevas Proteinas --------------------------------------------#
# # proteinas = [
# #       o_pro,
# #       tr_pro,
# #       tl_pro,
# #       both_pro,
# #     mod_pro
# # ]
# # nombres = ["Ori", "Trans", "Transl", "Both", "M_G"] # "M_G"
# # graficar_heatmap(proteinas, nombres)
# # graficar_mapa_identidad(proteinas, nombres)
# # """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# #
# # """
# #     Llamado a la función para guardar en archivo FASTA
# #     en el cual se guardan todas las secuencias generadas, incluida la original
# #     solo se ejecuta una vez. Luego de usado "silenciarlo"
# # """
# # secuencias = [sequence, seq_tr, seq_tl, seq_both, seq_mod]
# # guardar_fasta("Funciones/secuencias_mutadas.fasta", nombres, secuencias)
#
# """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# """ !!!!!!!!! Ejecutar las siguientes lineas de codigo una a la vez !!!!!!!!!!!!!!!
# Una vez creado el archivo con las secuencias del proyecto, se realiza blast de la siguiente manera:
#     1. Verificación del porcentaje de identidad
#     2. Filtrar unicamente secuencias con un porcentaje mayor al 90%
# """
#
#
#
# #query_file = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Funciones\\secuencias_mutadas.fasta"
# #database = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Base_de_datos\\alphainfluenza_H_y_N_Genbank_data_base.fasta"
# #output = "Resultado_verificación_secuencia_BLASTTT.txt"
# # # Llamar a la función
# #generar_db_blast(query_file, database, output)
#
#
# #     2. Verificacion de que todas las secuencias de Hemaglutinina coinciden al 100%
# # """
# #query_file = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Validación\\Neuraminidasa\\Pato_N_H5N1_AA.fasta"
# query_file = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Validación\\IA.fasta"
# database = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\Base_de_datos\\alphainfluenza_H_y_N_Genbank_data_base.fasta"
# output = "Resultado_verificación_secuencia_BLASTT.txt"
# blast_verification(query_file, database, output) # Llamar a la función
# #
#
# """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# """
# Despues de confrimar el porcentaje de identidad:
#     1. Se incluyen estas secuencias en un archivo fasta que contiene las secuencias de las 18 hemaglutininas o
#     de las 10 neuramidasas segun sea el caso.
#     2. Con esa totalidad de secuencias se realiza un analisis de filogenia, generando un archivo .iqtree
#     que posteriormente sera analizado en la aplicacion que lleva el mismo nombre.
# """
#
#
# # Parámetros de entrada
# # query_file = "C:\\Users\\LENOVO\Desktop\\Entrega_03_Influenza_Analyzer\\Hemaglutininaaaaz.PHY"
# # iq_tree_app = "C:\\Users\\LENOVO\\Desktop\\Entrega_03_Influenza_Analyzer\\iqtree-1.6.12-Windows\\bin\\iqtree.exe"
#
# # Llamar a la función
# # phylogenetic_tree(query_file, iq_tree_app)

print(mostrar_ayuda())
