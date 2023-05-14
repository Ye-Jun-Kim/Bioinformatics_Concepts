#Locating the ori of E Coli through data visualization

import vibrio_cholerae as vib
import matplotlib.pyplot as plt

with open('e_coli_genome.txt') as file:
    e_coli = file.read()

array = vib.FasterSymbolArray(e_coli, "C")

plt.plot(*zip(*sorted(array.items())))
plt.show()

#The plot shows that there is a maximum number of Cytosine around position 1600000, which likely notes that the lagging strand
#begins replicating there. The minimum number of Cytosine happens around position 4000000, which likely notes that 
#the leading strand begins replicating there. Because the ori is located where the lagging strand transitions into 
#the leading strand, we can deduce that the replication origin is around position 4000000 for E.Coli