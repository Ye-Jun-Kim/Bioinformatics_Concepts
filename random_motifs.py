#continues on from motifs.py, but implements the random library to iterate randomness into the motif search
#initiating the best motif search through randomness might seem counterproductive at first,
#but it works due to the fact that strings in DNA are not random and include implanted motifs.

#Up until so far, we have built functions that reference the first given motif as the beginning of the search,
#but nothing guarantees that the first given motif is a good point of reference. 
#Therefore, it is in the interest of scientists to explore the results of implementing randomness to the motif search
#in order to generate a realistic result.

#First we import the necessary random library and the functions from the motif file
import random 
import motifs as m


#Built a function to select random motifs, given a list of DNA sequences. The result is a list of random k-mer motifs 
#selected from a randomized value.
def RandomMotifs(Dna, k, t):
    random_motifs = []
    t = len(Dna)
    l = len(Dna[0])
    for i in range(t):
        ran_num = random.randint(0,l-k)
        random_motifs.append(Dna[i][ran_num:ran_num+k])
    return random_motifs

#Takes the previous function and generates the best motif depending on the Score function built in the motifs file
#usually ran many, many times in order to generate the best result after many sample sizes.
def RandomizedMotifSearch(Dna,k,t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = m.ProfileWithPseudocounts(M)
        M = m.Motifs(Profile, Dna,k)
        if m.Score(M) < m.Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 


#In order to normalize results, we can create a function to normalize all probabilities reletive to 1, by 
#taking a dictionary of probabilities similar to a result of running the Count and Profile functions from motifs
def Normalize(Probabilities):
    prob_sum = sum(Probabilities.values())
    new_dict = {}

    for key,value in Probabilities.items():
        new_dict[key] = value/prob_sum
    
    return new_dict

def WeightedDie(Probabilities):
    count = 0
    n = random.uniform(0,1)
    for keys,values in Probabilities.items():
        count = count+values
        if n < count:
            return keys
        

def ProfileGeneratedString(text,profile,k):
    n = len(text)
    probabilities = {}
    for i in range(0,n-k+1):
        probabilities[text[i:i+k]] = m.Pr(text[i:i+k], profile)

    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

def GibbsSampler(Dna, k, t, N):
    best_motifs = []
    motifs = RandomMotifs(Dna, k, t)
    best_motifs = motifs
    for j in range(1,N):
        i = random.randint(0,t-1)
        reduce_motifs = []
        for j in range(0,t):
            if j != i:
                reduce_motifs.append(motifs[j])
        Profile = m.ProfileWithPseudocounts(reduce_motifs)
        motif_i = ProfileGeneratedString(Dna[i], Profile, k)
        motifs[i] = motif_i
        if m.Score(motifs) < m.Score(best_motifs):
                best_motifs=motifs
    return best_motifs




test_dict = {'A':0.45,'B':0.63, 'C':0.09,'D':0.27,'E':0.36}
print(Normalize(test_dict))