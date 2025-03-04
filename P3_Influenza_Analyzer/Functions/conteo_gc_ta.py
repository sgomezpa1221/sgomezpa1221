def gc_content(sequence):
  """Calcula el porcentaje de contenido GC en una secuencia de nucleótidos."""
  sequence = sequence.upper()
  sequence =  sequence.rstrip("\n")
  G_count = sequence.count("G")
  C_count = sequence.count("C")
  gc_content = (G_count + C_count) / len(sequence) * 100
  return gc_content

def ta_content(sequence):
  """Calcula el porcentaje de contenido GC en una secuencia de nucleótidos."""
  sequence = sequence.upper()
  sequence =  sequence.rstrip("\n")
  T_count = sequence.count("T")
  A_count = sequence.count("A")
  gc_content = (T_count + A_count) / len(sequence) * 100
  return gc_content
