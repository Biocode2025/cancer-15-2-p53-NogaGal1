import random
RNA_codon_table = {}

#function that replaces a nucleotid
def Mutate_DNA(seq):
  Mutated_seq = ""
  nuc_list = ["A", "T", "G", "C"]
  seq = seq.upper()
  
  #replace a nucleotid in the given sequence with a random nucleotid
  rand_index = random.randrange(0, len(seq))
  rand_nuc = seq[rand_index]
  nuc_list.remove(rand_nuc)
  rand_base = random.randrange(0, len(nuc_list))
  new_nuc = nuc_list[rand_base]
  
  #create & combine the parts
  fst_part = seq[0 : rand_index]
  snd_part = seq[rand_index + 1 :]
  Mutated_seq = fst_part + new_nuc + snd_part
  
  return Mutated_seq

#function that adds a nucleotid 
def Insert_DNA(seq):
  nuc_list = ["A", "T", "G", "C"]
  seq = seq.upper()
  
  #select a random index and a random nucleotid
  rand_index = random.randrange(0, len(seq)) 
  new_nuc = random.choice(nuc_list)     
    
  #create & combine the parts
  fst_part = seq[:rand_index]
  snd_part = seq[rand_index:]
  mutated_seq = fst_part + new_nuc + snd_part

  return mutated_seq

#function that deletes a nucleotid
def Delete_DNA(seq):
  seq = seq.upper()
  
  rand_index = random.randrange(0, len(seq)) 
  
  #create & combine the parts
  fst_part = seq[:rand_index]
  snd_part = seq[rand_index + 1:] #ignores a nucleotid (deletes it)

  mutated_seq = fst_part + snd_part

  return mutated_seq

#compare the sequences, return number of diffrences
def Comp_seq(a, b):
  count = 0
  for x, y in zip(a, b):
    if x != y:
      count += 1
  return count

#translates codons to amino acids
def Read_dict():
  global RNA_codon_table
  with open("codon_AA.txt") as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      codon, aa = line.split()  
      RNA_codon_table[codon] = aa

#transcribes DNA to RNA 
def DNA_RNA_Cod(dna_seq):
  rna_seq = dna_seq.upper().replace("T", "U")
  return rna_seq

#translates RNA to protein
def RNA_prot(rna_seq):
  protein = ""
  for i in range(0, len(rna_seq)-2, 3):
    codon = rna_seq[i:i+3]
    aa = RNA_codon_table[codon]
    protein += aa
  return protein
