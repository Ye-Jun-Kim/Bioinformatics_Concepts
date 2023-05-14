#This file imports the vibrio_cholerae file in order to use its functions for analyzing the known ori of Thermotaga Petrophilia

import vibrio_cholerae as vib 

#The known ori of other species can be analyzed in order to compare the ori of the species being studied.

#Thermotoga Petrophilia is a bacterium that thrives in extreme temperatures. It also has a known DnaA box:

#Given the DnaA box as a string, capitalizing all letters would make the functions from vibrao_cholerae more convenient to use
thermotoga_petrophilia_ori = "aactctatacctcctttttgtcgaatttgtgtgatttatagagaaaatcttattaactga\
aactaaaatggtaggtttggtggtaggttttgtgtacattttgtagtatctgatttttaa\
ttacataccgtatattgtattaaattgacgaacaattgcatggaattgaatatatgcaaa\
acaaacctaccaccaaactctgtattgaccattttaggacaacttcagggtggtaggttt\
ctgaagctctcatcaatagactattttagtctttacaaacaatattaccgttcagattca\
agattctacaacgctgttttaatgggcgttgcagaaaacttaccacctaaaatccagtat\
ccaagccgatttcagagaaacctaccacttacctaccacttacctaccacccgggtggta\
agttgcagacattattaaaaacctcatcagaagcttgttcaaaaatttcaatactcgaaa\
cctaccacctgcgtcccctattatttactactactaataatagcagtataattgatctga".upper()


#Using the PatternCount function from vibrio_cholerae,
#We can analyze if the Thermotaga Petrophilia DnaA box shares the same significant 9-mers as Vibrio Cholerae.
count_1 = vib.PatternCount(thermotoga_petrophilia_ori, 'ATGATCAAG')
count_2 = vib.PatternCount(thermotoga_petrophilia_ori,'CTTGATCAT')
print(count_1, count_2)

#we can see from the printed results that these 9-mers do not exist in the ori of Thermotoga Petrophilia
#It might be obvious, but we can conclude from this finding that different genomes have different DnaA boxes
