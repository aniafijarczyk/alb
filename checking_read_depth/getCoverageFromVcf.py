import gzip
import numpy as np

fname = "test.vcf.gz"
faifile = "Agla_Btl03082013.genome.fa.fai"

def getTag(string_,tag_):
    dp = int(string_.split(':')[tag_])
    return(dp)

def getDP4(info_):
    info = info_.split(';')
    dp4 = [i.split('=')[1].split(',') for i in info if i[0:3] == 'DP4'][0]
    sum_dp4 = sum([int(j) for j in dp4])
    return(sum_dp4)

def findCumulative(faifile_):
    #fai = pd.read_csv(faifile_, sep = "\t", header = None, names = ["chrom","length","x1","x2","x3"])
    fh = open(faifile_,"r")
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    D = {a[0]:int(a[1]) for a in k}
    nums = list(D.values())
    N = []
    i = 0
    for ele in nums:
        N.append(i)
        i+=ele
    #ND = {a:[c,b] for a,b,c in zip(list(D.keys()),N,nums)}
    #dN = pd.DataFrame({'chrom':list(D.keys()),'length':nums,'addition':N})
    dN = list(zip(list(D.keys()),nums,N))
    D = {chrom:[leng,cum] for chrom,leng,cum in dN}
    return(D)

if __name__ == '__main__':

    # reading file
    fh = gzip.open(fname,"rt")
    linie = fh.readlines()
    k = [ele.split() for ele in linie if ele[0] != "#"]

    # check DP field
    tags = k[0][8].split(':')
    tags_ind = [j for j in range(len(tags)) if tags[j] == 'DP'][0]

    # check number of inds
    gts = k[0][9:]
    gts_num = len(gts)

    # getting DP4
    S = {}
    for snp in k:
        info_dp4 = getDP4(snp[7])
        ind_dp = []
        for ind in snp[9:]:
            dpi = getTag(ind,tags_ind)
            ind_dp.append(dpi)
        S[snp[0]+'_'+snp[1]] = [info_dp4] + ind_dp


    # getting scaffold lengths
    df = findCumulative(faifile)

    # writing output to file
    wh = open(fname.split('.')[0]+'_cov.tab','w')
    header = ['chrom','pos','pos_cum','length','DP4','\t'.join(['pool_'+str(i) for i in range(gts_num)])]
    positions = list(S.keys())
    subset = positions[0::1000]
    wh.write('\t'.join(header)+'\n')
    for ele in subset:
        inds = '\t'.join([str(ind) for ind in S[ele]])
        chrom = ele.split('_')[0]
        pos = ele.split('_')[1]
        n = chrom, pos, str(int(pos)+df[chrom][1]), str(df[chrom][0]), inds
        f = '\t'.join(n) + '\n'
        wh.write(f)
    wh.flush()
    wh.close()

    # calculating stats
    dp = list(zip(*S.values()))
    print("mean = "+str(np.mean(dp[0])))
    print("median = "+str(np.median(dp[0])))
    print("min = "+str(min(dp[0])))
    print("max = "+str(max(dp[0])))
    print("std = "+str(np.std(dp[0])))
    print("mean*2 = "+str(np.mean(dp[0])*2))
    print("perc99 = "+str(np.percentile(dp[0],99)))
    print("mean+4*sqrt(mean) = "+str(np.mean(dp[0]) + 4*np.sqrt(np.mean(dp[0]))))





