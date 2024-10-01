# --- QUESTION ONE ---
def translate_dna_to_mrna_to_amino(dna_sequence):

    seq = dna_sequence.upper()
    print("Input DNA: " + seq + "\n")

    # checks if DNA sequence is valid
    valid = seq.count("A") + seq.count("C") + seq.count("G") + seq.count("T")
    if valid == len(seq):
        pass
    else:
        print("DNA sequence invalid!")
        return False

    # check if multiple of 3
    if valid % 3 == 0:
        pass
    else:
        print("DNA sequence is not a multiple of 3!")
        return False
    
    # computes the reverse complement in the DNA sequence
    comp = ""

    for c in seq:
        if c == "A":
            comp += "T"
        elif c == "T":
            comp += "A"
        elif c == "G":
            comp += "C"
        elif c == "C":
            comp += "G"
    
    print("Complement = " + comp)

    # find the mRNA
    new_seq = comp.replace("T", "U")
    print("mRNA = " + new_seq)

    # translate into amino acid

    # codon table
    codon_table = {
        'AUG': ('Methionine', 'M'), 'UUU': ('Phenylalanine', 'F'), 'UUC': ('Phenylalanine', 'F'),'UUA': ('Leucine', 'L'), 'UUG': ('Leucine', 'L'), 'UCU': ('Serine', 'S'), 'UCC': ('Serine', 'S'),
        'UCA': ('Serine', 'S'), 'UCG': ('Serine', 'S'), 'UAU': ('Tyrosine', 'Y'), 'UAC': ('Tyrosine', 'Y'),'UGU': ('Cysteine', 'C'), 'UGC': ('Cysteine', 'C'), 'UGG': ('Tryptophan', 'W'), 'CUU': ('Leucine', 'L'),
        'CUC': ('Leucine', 'L'), 'CUA': ('Leucine', 'L'), 'CUG': ('Leucine', 'L'), 'CCU': ('Proline', 'P'),'CCC': ('Proline', 'P'), 'CCA': ('Proline', 'P'), 'CCG': ('Proline', 'P'), 'CAU': ('Histidine', 'H'),'CAC': ('Histidine', 'H'), 'CAA': ('Glutamine', 'Q'), 'CAG': ('Glutamine', 'Q'), 'CGU': ('Arginine', 'R'),
        'CGC': ('Arginine', 'R'), 'CGA': ('Arginine', 'R'), 'CGG': ('Arginine', 'R'), 'AUU': ('Isoleucine', 'I'),'AUC': ('Isoleucine', 'I'), 'AUA': ('Isoleucine', 'I'), 'GUU': ('Valine', 'V'), 'GUC': ('Valine', 'V'),
        'GUA': ('Valine', 'V'), 'GUG': ('Valine', 'V'), 'GCU': ('Alanine', 'A'), 'GCC': ('Alanine', 'A'),'GCA': ('Alanine', 'A'), 'GCG': ('Alanine', 'A'), 'GAU': ('Aspartic acid', 'D'), 'GAC': ('Aspartic acid', 'D'),
        'GAA': ('Glutamic acid', 'E'), 'GAG': ('Glutamic acid', 'E'), 'GGU': ('Glycine', 'G'), 'GGC': ('Glycine', 'G'),'GGA': ('Glycine', 'G'), 'GGG': ('Glycine', 'G'), 'UAA': ('Stop', None), 'UAG': ('Stop', None),
        'UGA': ('Stop', None), 'ACU': ('Threonine', 'T'), 'ACC': ('Threonine', 'T'), 'ACA': ('Threonine', 'T'),'ACG': ('Threonine', 'T'), 'AAU': ('Asparagine', 'N'), 'AAC': ('Asparagine', 'N'), 'AAA': ('Lysine', 'K'),
        'AAG': ('Lysine', 'K'), 'AGU': ('Serine', 'S'), 'AGC': ('Serine', 'S'), 'AGA': ('Arginine', 'R'),'AGG': ('Arginine', 'R')
    }

    # protein letter table to show the letters in the amino acid
    protein_letters = {
        "A": "Ala","C": "Cys","D": "Asp","E": "Glu","F": "Phe","G": "Gly","H": "His","I": "Ile","K": "Lys","L": "Leu","M": 
        "Met","N": "Asn","P": "Pro","Q": "Gln","R": "Arg","S": "Ser","T": "Thr","V": "Val","W": "Trp","Y": "Tyr",
    }

    # separate mRNA every 3 letters
    mrna_split = [(new_seq[i:i+3]) for i in range(0, len(new_seq), 3)]
    # print(mrna_split)

    # replace the mRNA with the corresponding colon table
    amino_acid_list = []
    for x in mrna_split:
     
        # searching from codon table
        amino_acid_list.append(codon_table.get(x, x))
     
    # print("Aminoacid = " + str(amino_acid_list))

    # replace the codon table with shorter names from the protein letter table for format accuracy    
    shorter_amino_acid_list = []

    for elements in amino_acid_list:
        for items in elements:
            shorter_amino_acid_list.append(protein_letters.get(items))
            while(None in shorter_amino_acid_list):
                shorter_amino_acid_list.remove(None)
                
    # print(shorter_amino_acid_list)

    # just for formatting
    final_result = ""
    for a in range(0, len(shorter_amino_acid_list)):
        final_result += shorter_amino_acid_list[a] + " (" + amino_acid_list[a][1] + ")" + " - "
    final_result = final_result[:-2]
    print("Aminoacid = " + final_result)



# --- QUESTION TWO ---
def freq_of_RNA_codon_in_aminoacid(amino_acid):
    protein_letters = {
        "A": ["GCC", "GCA", "GCG", "GCU"],"C": ["UGU", "UGC"],"D": ["GAU", "GAC"],"E": ["GAA", "GAG"],
        "F": ["UUU", "UUC"],"G": ["GGU", "GGC", "GGA", "GGG"],"H": ["CAU", "CAC"],"I": ["AU", "AUC", "AUA"],
        "K": ["AAA", "AAG"],"L": ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"],"M": ["AUG"],"N": ["AAU", "AAC"],
        "P": ["CCU", "CCC", "CCA", "CCG"],"Q": ["CAA", "CAG"],"R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],"T": ["ACU", "ACC", "ACA", "ACG"],"V": ["GUU", "GUC", "GUA", "GUG"],
        "W": ["UGG"],"Y": ["UAU", "UAC"],
    }

    amino_acid = amino_acid.upper()

    # checks if there are 3 aminoacids
    if len(amino_acid) <= 3:
        pass
    else:
        print("There can only be 3 amino acids maximum!")

    # checks validity
    for letter in amino_acid:
        if letter in protein_letters:
            pass
        else:
            print("Amino acid is invalid!")

    # create all possible combinations of codons manually :(
    codon_combinations = []
    if len(amino_acid) == 1:
        for first_letter in protein_letters[amino_acid[0]]:
            codon_combinations.append([first_letter])
    elif len(amino_acid) == 2:
        for first_letter in protein_letters[amino_acid[0]]:
            for second_letter in protein_letters[amino_acid[1]]:
                codon_combinations.append([first_letter, second_letter])
    elif len(amino_acid) == 3:
        for first_letter in protein_letters[amino_acid[0]]:
            for second_letter in protein_letters[amino_acid[1]]:
                for third_letter in protein_letters[amino_acid[2]]:
                    codon_combinations.append([first_letter, second_letter, third_letter])

    # count the frequency of each codon in the combination

    for combination in codon_combinations:
        codon_frequency = {}
        codon_sequence = "".join(combination)
        print("mRNA = " + codon_sequence)

        for codon in combination:
            if codon in codon_frequency:
                codon_frequency[codon] += 1
            else:
                codon_frequency[codon] = 1
        
        # print the final count
        for key, value in codon_frequency.items():
            print(str(key) + " = " + str(value))

# user interface
def user_interface():
    while True:
        answer = input("""Pick an Option:
1. Translate a DNA sequence -> mRNA (using “U” instead of ”T”) -> into an aminoacid sequence (protein)
2. Provides the frequency of each RNA codon encoding a given aminoacid, in a DNA sequence
3. Exit the program\n""")
        if answer == "1":
            answer_one = input("Insert the DNA sequence here: ")
            translate_dna_to_mrna_to_amino(answer_one)
            return False
        elif answer == "2":
            answer_two = input("Insert amino acid here: ")
            freq_of_RNA_codon_in_aminoacid(answer_two)
            return False
        elif answer == "3":
            print("Thank you for using the program!")
            return False
        else:
            print("That option doesn't exist here, please try running the program again.\n")
            return False

user_interface()