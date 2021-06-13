import itertools
from collections import defaultdict
import random
import numpy as np


fname = 'sample_sliding_window_fst.txt'


def randomPairs(combinations,pairs=5):
    tmp = combinations
    P = []
    S = []
    i = 1
    while i < int(pairs+1):
        t1 = random.sample(tmp,1)[0]
        if (t1[0] not in P) & (t1[1] not in P):
            #print(t1)
            S.append('_'.join(t1))
            tmp = [ele for ele in tmp if ele != t1]
            P.append(t1[0])
            P.append(t1[1])
            i+=1
        else:
            continue 
    return(S)

if __name__ == '__main__':

    north = ['p2','p3','p8','p9','p10']
    south = ['p1','p4','p5','p6','p7']
    ns = list(itertools.product(south,north))

    # Reading file
    fh = open(fname,'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    header = k[0]
    D = {}
    body = list(zip(*k[1:]))
    for i in range(len(header)):
        D[header[i]] = body[i]
    D.keys()

    # Calculating mean fst between random nonredundant pairs of north & south, resampled 10 times
    W = []
    for win in range(len(D['Chromosome'])):
        #print(win)
        boot = []
        for n in range(10):
            b = []
            pairs = randomPairs(ns)
            for pair in pairs:
                b.append(float(D[pair][win].replace('na','NaN')))
            row_mean = np.mean(b)
            boot.append(row_mean)
        boot_mean = np.mean(boot)
        W.append(boot_mean)

    
    # Adding column north_south in the output
    nname = fname.replace('.txt','')
    
    wh = open(nname+'_mean.txt','w')
    wh.write('\t'.join(header)+'\tnorth_south\n')
    for row in range(len(k[1:])):
        nrow = k[1:][row] + [W[row]]
        #nk.append(nrow)
        wh.write('\t'.join(k[1:][row]) + '\t' + str(round(W[row],8)) + '\n')
    wh.flush()
    wh.close()  
        
        
