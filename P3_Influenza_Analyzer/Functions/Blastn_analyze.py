import subprocess

def blast_verification(query, subject_db, output_file):
    output_format = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
    command_blast = [
         "blastn",
         "-query", query,
         "-db", subject_db,
         "-outfmt", output_format,
         "-out", output_file
     ]
    print("Ejecutando BLAST...")
    subprocess.run(command_blast, shell=True)
    print("BLASTN finalizado.")
