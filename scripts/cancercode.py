import random
import os

RNA_codon_table = {}

# Set up paths
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(script_dir, '..', 'data'))
results_dir = os.path.normpath(os.path.join(script_dir, '..', 'results'))

####functions####

# function that replaces a nucleotide
def Mutate_DNA(seq):
    if len(seq) == 0:
        return seq
    nuc_list = ["A", "T", "G", "C"]
    seq = seq.upper()

    rand_index = random.randrange(len(seq))
    rand_nuc = seq[rand_index]
    nuc_list.remove(rand_nuc)
    new_nuc = random.choice(nuc_list)

    return seq[:rand_index] + new_nuc + seq[rand_index + 1:]

# function that adds a nucleotide
def Insert_DNA(seq):
    if len(seq) == 0:
        return seq
    nuc_list = ["A", "T", "G", "C"]
    seq = seq.upper()

    rand_index = random.randrange(len(seq) + 1)
    new_nuc = random.choice(nuc_list)

    return seq[:rand_index] + new_nuc + seq[rand_index:]

# function that deletes a nucleotide
def Delete_DNA(seq):
    if len(seq) == 0:
        return seq
    seq = seq.upper()

    rand_index = random.randrange(len(seq))
    return seq[:rand_index] + seq[rand_index + 1:]

# compare sequences 
def Comp_seq(a, b):
    count = abs(len(a) - len(b))
    for x, y in zip(a, b):
        if x != y:
            count += 1
    return count

# reads codon table
def Read_dict():
    global RNA_codon_table
    with open(os.path.join(data_dir, "codon_AA.txt")) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                codon = parts[0].strip()
                aa = parts[1].strip()
                RNA_codon_table[codon] = aa

# transcribes DNA to RNA
def DNA_RNA_Cod(dna_seq):
    return dna_seq.upper().replace("T", "U")

# translates RNA to protein 
def RNA_prot(rna_seq):
    protein = ""
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        if codon not in RNA_codon_table:
            break
        aa = RNA_codon_table[codon]
        if aa == "*":  # stop codon
            break
        protein += aa
    return protein

#### main code #####

def main():
    Read_dict()

    # read p53 DNA
    orgnl_seq = ""
    with open(os.path.join(data_dir, "p53_sequence.fa"), "r") as f:
        for line in f:
            line = line.strip().upper()
            if line.startswith(">"):
                continue
            orgnl_seq += line

    mutated_seq = orgnl_seq
    mutations = 3

    for _ in range(mutations):
        random_func = random.randrange(1, 4)
        if random_func == 1:
            mutated_seq = Mutate_DNA(mutated_seq)
        elif random_func == 2:
            mutated_seq = Insert_DNA(mutated_seq)
        else:
            mutated_seq = Delete_DNA(mutated_seq)

    diff_num = Comp_seq(orgnl_seq, mutated_seq)
    percent_diff = diff_num / max(len(orgnl_seq), len(mutated_seq)) * 100

    p53_orig = RNA_prot(DNA_RNA_Cod(orgnl_seq))
    mut_p53_prot = RNA_prot(DNA_RNA_Cod(mutated_seq))
    p53_diff = Comp_seq(p53_orig, mut_p53_prot)

    # write output
    with open(os.path.join(results_dir, "mutated_p53.fasta"), "w") as f:

        f.write(">original p53 protein\n")
        f.write(p53_orig + "\n\n")

        f.write(">mutated p53 protein\n")
        f.write(mut_p53_prot + "\n\n")


main()