#Genetic motifs are proteins or DNA sequences with a specific structure or pattern

#In order to visualize a clear trend of a motif given a set of genetic sequences, we can create a dictionary of lists
#counting each respective nucleotide in the given sequences.

#We can start by creating a base function to count the number of nucleotides in each str given an X amount of motifs,
#assuming that the given motifs are of the same length of characters.
#(Argument should be a LIST of STRINGS, and returns a dictionary of lists)
def Count(motif):
    #The following series of steps initializes a dictionary of lists. Dictionary contains 4 lists, each corresponding to the 
    #count of the nucleotides counted at each given index of the motifs.
    count = {}
    k = len(motif[0])
    for nuc in 'ATCG':
        count[nuc] = []
        for i in range(k):
            count[nuc].append(0)
    
    #Goes through each letter in the motif to count the number of nucleotides at each index
    t = len(motif)
    for i in range(t):
        for num in range(k):
            nuc = motif[i][num]
            count[nuc][num] = count[nuc][num] + 1

    return count


#We can build a function to return a new dictioanry with values derived from the previous function, 
#but with a proportional value out of the total number of strings/motifs given to the function.
#It allows user to quantify per capita of how many nucleotides there are aside from other nucleotides.
#(Argument is a LIST OF STRINGS; function returns a dictionary)
def Profile(motif):
    k = len(motif[0])
    table = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
    for i in range(k):
        for j in range(len(motif)):
            nuc = motif[j][i]
            table[nuc][i] += 1
    return table


#We can continue on to build a function returning a STRING of the highest called nucelotides from the given motifs
#(the nucleotide with the highest proportion)
#(Argument motif is a list of strings; returns a string)
def Consensus(motif):
    table = Count(motif)
    consensus_str = ""
    for i in range(len(motif[0])):
        m = 0
        frequentSymbol = ""
        for symbol in "ATCG":
            if table[symbol][i] > m:
                m = table[symbol][i]
                frequentSymbol = symbol
        consensus_str = consensus_str + frequentSymbol
    
    return consensus_str

#Now we elaborate from the last function to return an evaluated score of the consensus string
#(Argument is a list of strings; returns an integer)
def Score(motif):
    consensus_str = Consensus(motif)
    score_val = 0

    for i in motif:
        num = 0
        for letter in i:
            if letter != consensus_str[num]:
                score_val = score_val + 1
            num = num + 1
    
    return score_val 
    

#Fun extra exercise to return the same function as Consensus, but return a List instead of a String
#(Argument is a list of strings; returns a list of strings)
def ConsensusList(motif):
    table = Profile(motif)
    consensus_tab = []
    nuc='ATCG'
    for x in range(len(table[nuc[0]])):
        consensus_tab.append('X')

    for i in range(len(nuc)):
        for j in range(len(table[nuc[i]])):
            if table[nuc[i]][j] > table[nuc[i-1]][j]:
                consensus_tab[j] = nuc[i]
    
    return consensus_tab


#A probability function can be created to calculate the probability of a certain sequence/profile happening.
#Since we already have a function to generate the evaluated score of each consensus, we just need to multiply our
#probability value by the subsequent probaiblity chances to get the total probability chance.
#(Arguments are a string,dictionary of lists; returns a float)
def Pr(text,profile):
    prob = 1

    for i in range(len(text)):
        prob = prob * profile[text[i]][i]
    
    return prob 


#We can now use the previous functions to create a new function that returns the most likely profile occurring in a 
#following sequence of 4 x k matrix profile
#Aguments are a string,integer,dictionary of lists; returns a string
def ProfileMostProbablePattern(text,k,profile):
    temp = 0
    result = ""

    for i in range(len(text)-k+1):
        if Pr(text[i:i+k],profile) > temp:
            temp = Pr(text[i:i+k],profile)
            result = text[i:i+k]
    
    #if all specified k-mers returns 0, the first k-mer identified will be given by default
    if len(result) == 0:
        result = text[0:k]

    return result
           
#The following function goes through a given strings of DNA, and returns lists of the best motifs identified through
#the scoring functions written above.
#(Arguments are a list of strings,k-mer integer, t-integer of len of the list: returns a lists of strings)
def GreedyMotifSearch(Dna,k,t):
    best_motifs = []
    
    for i in range(0, t):
        best_motifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        
        if Score(Motifs) < Score(best_motifs):
            best_motifs = Motifs
    
    return best_motifs

#Although a more complicated function, GreedyMotifSearch actually trades efficiency for speed. It uses the first identified 
#motif as a reference to the rest of the motifs. The problem with this is outlined in the 'DosR.py' file.

def CountWithPseudocounts(motif):
    #initializes a dictionary with values as 0 for all locations
    count = {}
    k = len(motif[0])
    for nuc in 'ATCG':
        count[nuc] = []
        for i in range(k):
            count[nuc].append(0)
    
    #loops through all dictiionary values and adds 1 count if a certain nucleotide is read
    t = len(motif)
    for i in range(t):
        for num in range(k):
            nuc = motif[i][num]
            count[nuc][num] = count[nuc][num] + 1

    #now that the initial dictionary is complete, we can loop through again to replace 0s with pseudocounts      
    for nuc in 'ATCG':
        if 0 in count[nuc]:
            for nuc in 'ATCG':
                for i in range(len(count[nuc])):
                    count[nuc][i] += 1
                            

    return count

#Profiles the given motifs with pseudocounts in place, updated in the CountWithPseudocounts function above
def ProfileWithPseudocounts(motif):
    count_dict = CountWithPseudocounts(motif)
    k = len(motif[0])
    new_dict = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
 
    
    for i in range(k):
        col_sum = 0
        for nuc in count_dict:
            col_sum = col_sum + count_dict[nuc][i]
        
        for nuc in new_dict:
            new_dict[nuc][i] = count_dict[nuc][i]/col_sum

        
    
    return new_dict


#Update our GreedyMotifSearch function, replacing the Profiling with the ProfileWithPseudocounts,
#in order to update the values with pseudocounts in place of 0s
#What results is a more accurate and helpful probability result of likely motifs
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    best_motifs = []
    
    for i in range(0, t):
        best_motifs.append(Dna[i][0:k])

    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        
        if Score(Motifs) < Score(best_motifs):
            best_motifs = Motifs
    
    return best_motifs


#Uses the ProfileMostProbablePattern and Pr functions written above in order to return the most probable 4-mer
#from a given a list of Dna sequences
def Motifs(Profile,Dna,k):
    motifs = []
    t = len(Dna)
    for i in range(t):
        motif = ProfileMostProbablePattern(Dna[i], k, Profile)
        motifs.append(motif)
    return motifs


test_result = ['CCA',\
'CCT',\
'CCT',\
'TTG']

test_dna = ['AAGCCAAA',\
'AATCCTGG',\
'GCTACTTG',\
'ATGTTTTG']
print(Motifs(Profile(test_result), test_dna,3))
