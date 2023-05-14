#Another method of visualizing where the ori is located is through visualizing the difference of Guanine and Cytosine
#Due to the likelihood of Cytosine mutating in the leading strand, and the liklihood of this causing an imbalance of G-C pairr
#An increasing difference of the number of G vs C signifies a potential leading strand, while
#A decreasing difference of the number of G vs C signifies a potential lagging stand.

import vibrio_cholerae as vib
import matplotlib.pyplot as plt

with open('e_coli_genome.txt') as file:
    genome = file.read()

skew = vib.SkewArray(genome)
plt.plot(skew)
plt.xlabel('Genome Position')
plt.ylabel('Skew')
plt.show()

#Similarylt to visualizing the difference of Cytosine in the E Coli genome, the data shows a clear trend of 
#a peak value at 1600000 and a valley at 4000000. 

#The skew should reach a minimum at the point the lagging strand ends and the leading strand begins.
#This minimum value is achieved at around the 4000000 position, promoting the idea that this is where the ori is for 
#E.Coli genome.