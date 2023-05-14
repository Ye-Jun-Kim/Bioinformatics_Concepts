
#This is the scientifically known replication of origin for Vibrio Cholerae
#for the purpose of this exercise, it will be a hypothesized ori
vibrio_cholerae_ori = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCAC\
ACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTC\
ATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTG\
TTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATC\
ATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGC\
CAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTT\
GATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCA\
TGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"

#Finding the length of the nucelotides making up the ori of Vibrio Cholerae
#print(len(vibrio_cholerae_ori))
#Vibrio Cholerae's Replication of origin is 539 base pairs long, out of a total 


#creating a function called PatternCount to count how many times a certain pattern is observed in a given strip of DNA
def PatternCount(text,pattern):
    count = 0
    for i in range(len(text)-len(pattern)+1):
        if text[i:i+len(pattern)] == pattern:
            count = count + 1
    return count


#Create a function called FrequencyMap to createa a dictionary of certain k-mers
def FrequencyMap(text,k_num):
    #first, need to create a dictionary of all the k-mers after a sepcified k-value
    freq = {}
    n = len(text)
    #initial for loop to go through each k-mer in the whole text once
    for i in range(n-k_num+1):
        pattern = text[i:i+k_num]
        count = 0
        #for loop within to take the initial k-mer and run through the text counting matching pairs. The final count is added onto the dictionary
        for num in range(n-k_num+1):
            if text[num:num+k_num] == pattern:
                count = count + 1
                freq[pattern] = count
    
    return freq



#Createa a function FrequentWords to output a list of most frequently used k-mers
def FrequentWords(text,k_num):
    words = []
    freq_dict = FrequencyMap(text,k_num)
    m = max(freq_dict.values())
    for key in freq_dict:
        if freq_dict[key] == m:
            words.append(key)

        words.sort()
        
    return words

#Testing to see if the functions above work as intended for the Vibrae cholerae oriC for a 10-mer
#print(FrequentWords(vibrio_cholerae_ori,10))


#DNA pairs Adenine(A) to Tyrosine(T) and Cytosine(C) to Guanine(G)
#DNA strands are always read from the 5' to 3' direction
#In order to analyze the complementary strand of a given order of nucelotides, two actions must be taken:
#1) the order of the nucleotides must be reversed to account for the flip of the 5' - 3' direction
#2) A-T and C-G and vise versa

#Create function to address the first problem, reversing to address the 5'-3' direction of DNA
def Reverse(text):
    reverse_text = ''
    for char in text:
        reverse_text = char + reverse_text

    return reverse_text


#Create a fcuntion to match the base pairs with its complementary counterpart
def Complement(text):
    comp = ''
    for char in text:
        if char == 'A':
            comp = comp + 'T'
        elif char == 'T':
            comp = comp + 'A'
        elif char == 'C':
            comp =comp + 'G'
        elif char == 'G':
            comp = comp + 'C'
    
    return comp
        
#Create a function to combine the past two functions to return a completed complementary strand of DNA
def ReverseComplement(text):
    new_text = Reverse(text)
    new_text = Complement(new_text)

    return new_text


#In order to make sure that the DnaA box is only contained within the replication of origin, you needs to
#check if it exists within the entire genome.
#Create a function to count the position index of a given pattern in the entire genome and compile into a list
def PatternMatching(pattern,genome):
    positions = []
    for i in range(len(genome)-len(pattern)+1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(i)
        
    return positions

#Checking to see PatternMatching on the entire vibrio cholera genome imported into a txt file
vibrio_cholerae_genome = open('vibrio_cholerae_genome.txt')
vibrio_cholerae_str = vibrio_cholerae_genome.read()
#pattern_of_interest = 'CTTGATCAT'
#result = PatternMatching(pattern_of_interest,vibrio_cholerae_str)
#print(result)
#print(len(result))

#WARNING:As it calls upon PatternCount, as well as symbol being one letter, running this code on an entire genome will
#either take foreer, or crash the computer
def SymbolArray(genome, symbol):
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(extended_genome[i:i+(n//2)],symbol)
    return array


#Intimidating at first, but the very basic explanation of the following code is that instead of calling upon PatternCount
#in each single letter of the genome, it references the last letter and the first of the next letter. If the last letter 
#checked was a matching letter, it subtracts 1, and if the next letter is a matching letter, it adds 1. The array is already
#predetermined/initialized by PatternCount BEFORE the for loop.
def FasterSymbolArray(genome, symbol):
    array = {}
    n = len(genome)
    ExtendedGenome = genome + genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(genome[0:n//2],symbol)

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array


#Creating a function to count the difference between C and G on a given genome sequence
#The values will correspond to the position on the genome, and will be given a value based on the count of C and G in the sequence
def SkewArray(genome):
    genome = genome.upper()
    skew = []
    skew.append(0)
    n = len(genome)
    if genome[0] == 'C':
        skew.append(-1)
    elif genome[0] == 'G':
        skew.append(1)
    else:
        skew.append(0)
      
    for i in range(1,n):
        if genome[i] == 'C':
            skew.append(skew[i]-1)
        elif genome[i] == 'G':
            skew.append(skew[i]+1)
        else:
            skew.append(skew[i])
    
    return skew


#The SkewArray creates an array of skew values. As we are most interested at the point where the laading strand
#ends and the lagging strand begins, we can make a function to find the minimum value of the Skew array.
def MinimumSkew(genome):
    var = SkewArray(genome)
    position = []
    num = 0
    min_var = min(var)

    for i in range(len(var)):
        if var[i] == min_var:
            position.append(num)
        num = num + 1
        
    return position


#Hamming Distance is the number of mismatches between 2 given k-mers. 
def HammingDistance(p,q):
    mismatch_count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            mismatch_count = mismatch_count + 1
    
    return mismatch_count

#Creating a function to index all positions in the genome where the sections resembles a given pattern,
#with a HammingDistance of less than or equal to a value 'd'
def ApproximatePatternMatching(pattern,genome,d):
    position = []
    for i in range(len(genome)-len(pattern)+1):
        mismatch_count = HammingDistance(pattern,genome[i:i+len(pattern)+1])

        if mismatch_count <= d:
            position.append(i)
    
    return position

#Function to run the previous ApproximatePatternMatching function, and count how many times it counted  
def ApproximatePatternCount(pattern,genome,d):
    position = ApproximatePatternMatching(pattern,genome,d)
    count = len(position)

    return count 
