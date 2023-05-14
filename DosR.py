#DosR stands for Dormancy Surival Regulator, a transcription factor that controls gene expression under hypoxic conditions
#It is the reason why bacteria like Mycobacterium Tuberculosis Bacterium(MTB) that causes Tuberculosis can stay dormant
#inside the human body for so long.

#The DosR dataset consists of genes in the bacteria that heavilty change experssion when in a hypoxic environment.

#We can use the functions built in the motif.py file in order to narrow down on what the motif could be for the DosR genes.

import motifs

#readlines() is used because the file is a long text, but we require a list of sequences in order to use our functions.
with open('DosR_dataset.txt','r') as file:
    data = file.readlines()

#In order to clean up the data, we need to get rid of the built in '\n' codes.
clean_data = []
for seq in data:
    seq_data = seq.strip().replace('\n','')
    clean_data.append(seq_data)
   

#We have the functions and the text file available. 
#We can use the GreedyMotifSearch function to create a list of strings containing the most probable motifs from the DosR set
#Arguments will be the DosR dataset, with a k-mer length of 15

greedy_result = motifs.GreedyMotifSearch(clean_data,15,len(clean_data))

#Through the results collected by GreedyMotifSearch, we can analyze its shortcomings
#GreedyMotifSearch sacrifices efficiency for speed: it references the rest of the motif analysis by the first motif observed
#in the sequence. There are many problems with this: we dont know if the first observed motif is a good motif or not, and
#there are so many 0 possibilities when analyzing all motifs that the probability has a suspiciously high number of 0s.