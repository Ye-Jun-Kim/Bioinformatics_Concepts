p = 'GGGCCGTTGGT'
q = 'GGACCGTTGAC'

import motifs 

test_data = {'A':[0.4,0.3,0.0,0.1,0.0,0.9],
'C':[0.2,0.3,0.0,0.4,0.0,0.1],
'G':[0.1,0.3,1.0,0.1,0.5,0.0],
'T':[0.3,0.1,0.0,0.4,0.5,0.0]}

test_result = motifs.Pr('TCGGTA',test_data)

print(test_result)
