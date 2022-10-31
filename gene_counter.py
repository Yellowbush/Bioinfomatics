# Jason Sy
# B Bio 340: Intro to Bioinformatics
# 10/20/2022
#
# Parse Fasta file and count nucleotides

import sys

infile = ('', 'r')

#In prompt input fasta filename
ff = sys.argv[1]

#Function to parse through fasta
def parser(infile):
    seqs = []
    headers = []
    with open(infile) as f:
        sequence = ""
        header = None
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:-1])
                if header:
                    seqs.append([sequence])
                sequence = ""
                header = line[1:]
            else:
                sequence += line.rstrip()
        seqs.append([sequence])
    return headers, seqs

headers, seqs = parser(ff)

flat_seqs = [item for sublist in seqs for item in sublist]


#Function to count nucleotides in fasta
def count_nucelotides(sequence):
    nuc_dict = {
        "g" : "Guanine ",
        "c" : "Cytosine ",
        "a" : "Adenine ",
        "t" : "Thymine "}
    firstthree = sequence[:3]
    g = sequence.upper().count('G') 
    c = sequence.upper().count('C')
    a = sequence.upper().count('A')
    t = sequence.upper().count('T')

    gc = (g + c)
   
    total = (g + c + a + t)
    
    gc_content = (gc/total * 100)

    codons = (total//3)

  
    name_of_first_three = nuc_dict[firstthree[0]]
    name_of_first_three += nuc_dict[firstthree[1]]
    name_of_first_three += nuc_dict[firstthree[2]]


    return '\n Total = {}\n Codons = {} \n First Three Nucleotides = {}\n Name of First Three Nucleotides = {} \n GCs = {} %\n G = {}, C = {}, A = {}, T = {}.'.format(total, codons, firstthree, name_of_first_three, round(gc_content, 2), g, c, a, t, )

for header, seq in zip(headers, flat_seqs):
    print(header, count_nucelotides(seq))

