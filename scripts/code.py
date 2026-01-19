import random
RNA_codon_table = {}


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