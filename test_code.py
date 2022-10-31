# Jason Sy
# B Bio 340: Intro to Bioinformatics
# 10/30/2022
#
# Convert nucleotide to protein, and test mutation type and similarity


from itertools import dropwhile
import sys
from difflib import SequenceMatcher


#Test sequence for similarity, rounds to 2 decimal places
def similar(original, mutation):
    return round(SequenceMatcher(None, original, mutation).ratio() * 100, 2)




#Function to get mutation type
def get_mutation_type(original, mutation):
    mutation_type = ""
    
    #If protein doesn't change mutation is silent and therefore unlikely to be deleterious
    if original == mutation:
        mutation_type = "Silent : Not Deleterious"
        return mutation_type
    #If protein is missing an amino acid, mutation is frameshift and likely to be deleterious
    if len(original) > len(mutation):
        mutation_type = "Frameshift Deletion : Likely Deleterious"
        return mutation_type
    #If protein has amino acid inserted, mutation is frameshift and likely to be deleterious
    if len(original) < len(mutation):
        mutation_type = "Frameshift Insertion : Likely Deleterious"
        return mutation_type
    #If protein amino acids change mutation is missense and likely to be deleterious
    else:
        mutation_type = "Missense : Likely Deleterious"
        return mutation_type



#Function for to convert RNA to Amino Acids
def codon_table(sequence):
    #Dictionary for Codons
    codons = {
        "UUU":"Phe", "UUC":"Phe", "UUA":"Leu", "UUG":"Leu",
        "UCU":"Ser", "UCC":"Ser", "UCA":"Ser", "UCG":"Ser",
        "UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
        "UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
        "CUU":"Leu", "CUC":"Leu", "CUA":"Leu", "CUG":"Leu",
        "CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
        "CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
        "CGU":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
        "AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
        "ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
        "AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
        "AGU":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
        "GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
        "GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
        "GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
        "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly",}

    aa_profile = ""

    #Break RNA into List of 3 nucleotides
    rna = [sequence[i:i+3] for i in range(0, len(sequence), 3)]


    #For loop to begin defining protein with AUG as the start codon
    for codon in dropwhile(lambda k: k != 'AUG', (rna)):
        profile = codon
        codon.find("AUG")
        #Break loop when approached with stop Codon
        if codon in codons.keys():
            if codon == "UAG":
                aa_profile += "STOP"
                break
            if codon == "UAA":
                aa_profile += "STOP"
                break
            if codon == "UGA":
                aa_profile += "STOP"
                break
            profile = codons[codon]
        aa_profile += profile + "-"
        
    return aa_profile




#TESTING for Silent Mutation
print("Mutation Test for Silent Mutation")
original_sequence = "AUGUUUUCUUAUUGUCUUCCUCAUCGUUGA"
mutated_sequence = "AUGUUCUCCUACUGCCUACCCCACCGGUAA"
print("     Nucleotides of original sequence:       " + original_sequence)
print("     Nucleotides of mutated sequence:        " + mutated_sequence)
print("     Amino Profile of orignigal sequence:    " + codon_table(original_sequence))
print("     Amino Profile of mutated sequence:      " + codon_table(mutated_sequence))
print("     Mutation Type is:                       " + get_mutation_type(codon_table(original_sequence), codon_table(mutated_sequence)))
print("     Protein Similarity Percentage:          {}".format(similar(codon_table(original_sequence), codon_table(mutated_sequence))) + " %")
print()
print()

#TESTING for Frameshift Deletion Mutation
print("Mutation Test for Frameshift Deletion Mutation")
original_sequence2 = "AUGUUUUCUUAUUGUCUUCCUCAUCGUUGA"
mutated_sequence2 = "AUGUCUUAUUGUCUUCCUCAUCGUUGA"
print("     Nucleotides of original sequence:       " + original_sequence2)
print("     Nucleotides of mutated sequence:        " + mutated_sequence2)
print("     Amino Profile of orignigal sequence:    " + codon_table(original_sequence2))
print("     Amino Profile of mutated sequence:      " + codon_table(mutated_sequence2))
print("     Mutation Type is:                       " + get_mutation_type(codon_table(original_sequence2), codon_table(mutated_sequence2)))
print("     Protein Similarity Percentage:          {}".format(similar(codon_table(original_sequence2), codon_table(mutated_sequence2))) + " %")
print()
print()

#TESTING for Frameshift Insertion Mutation
print("Mutation Test for Frameshift Insertion Mutation")
original_sequence3 = "AUGUUUUCUUAUUGUCUUCCUCAUCGUUGA"
mutated_sequence3 = "AUGUUAUUUUCUUAUUGUCUUCCUCAUCGUUGA"
print("     Nucleotides of original sequence:       " + original_sequence3)
print("     Nucleotides of mutated sequence:        " + mutated_sequence3)
print("     Amino Profile of orignigal sequence:    " + codon_table(original_sequence3))
print("     Amino Profile of mutated sequence:      " + codon_table(mutated_sequence3))
print("     Mutation Type is:                       " + get_mutation_type(codon_table(original_sequence3), codon_table(mutated_sequence3)))
print("     Protein Similarity Percentage:          {}".format(similar(codon_table(original_sequence3), codon_table(mutated_sequence3)))  + " %")
print()
print()

#TESTING for Missense Mutation
print("Mutation Test for Missense Mutation")
original_sequence4 = "AUGUUUUCUUAUUGUCUUCCUCAUCGUUGA"
mutated_sequence4 = "AUGCUCUCUUAUUGUCUUCCUCAUCGUUGA"
print("     Nucleotides of original sequence:       " + original_sequence4)
print("     Nucleotides of mutated sequence:        " + mutated_sequence4)
print("     Amino Profile of orignigal sequence:    " + codon_table(original_sequence4))
print("     Amino Profile of mutated sequence:      " + codon_table(mutated_sequence4))
print("     Mutation Type is:                       " + get_mutation_type(codon_table(original_sequence4), codon_table(mutated_sequence4)))
print("     Protein Similarity Percentage:          {}".format(similar(codon_table(original_sequence4), codon_table(mutated_sequence4))) + " %")
print()
print()
