### SNP calling in pool-seq with bcftools

using samtools v1.9, bcftools v1.9, vcftools v0.1.17

#### Files

```
bam=/home/anna/Dropbox/Pool-seq/data/pool1.sorted.markdup.Scaffold1part.bam
fasta=/home/anna/Dropbox/Pool-seq/data/albgenome_Scaffold1.fa
```

#### SNP calling

```
samtools calmd -Aru -@2 ${bam} ${fasta} | bcftools mpileup -C50 -f ${fasta} -q4 \
-a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP -Ou - | bcftools call -mv -f gq -Oz -o pool1.vcf.gz -
tabix -p vcf pool1.vcf.gz
```

#### Annotation of bad quality snps
```
zcat pool1.vcf.gz | grep -v "^#" | awk '{print $1"\t"$2"\t"$2"\tsnp_"NR}' > pool1_annotations
bgzip pool1_annotations
tabix -s 1 -b 2 -e 3 pool1_annotations.gz
```
```
export PERL5LIB=/home/anna/github/vcftools/src/perl
bcftools view pool1.vcf.gz | vcf-annotate -f +/d=15/q=20/Q=10/-W/w=10 \
--fill-HWE --fill-type -n -a pool1_annotations.gz -c CHROM,FROM,TO,INFO/SNP_ID \
-d key=INFO,ID=SNP_ID,Number=1,Type=Integer,Description='SnpList' | bgzip -c > pool1_annotated.vcf.gz
```

```
StrandBias FLOAT Min P-value for strand bias (given PV4) [0.0001]
BaseQualBias FLOAT Min P-value for baseQ bias [1e-100]
MapQualBias FLOAT Min P-value for mapQ bias [0]
EndDistBias FLOAT Min P-value for end distance bias [0.0001]
a = MinAB INT Minimum number of alternate bases [2]
c = SnpCluster INT1,INT2 Filters clusters of 'INT1' or more SNPs within a run of 'INT2' bases []
D = MaxDP INT Maximum read depth [10000000]
d = MinDP INT Minimum read depth [15]
q = MinMQ INT Minimum RMS mapping quality for SNPs [20]
Q = Qual INT Minimum value of the QUAL field [10]
W = GapWin INT Window size for filtering adjacent gaps []
w = SnpGap INT SNP within INT bp around a gap to be filtered [10]
```
#### Filtering

```
bcftools view -i 'TYPE=="snp" & FORMAT/AD[:1]>=4 & FORMAT/AD[:0]>=4 & (DP4[0]+DP4[1]+DP4[2]+DP4[3])>=15 & (DP4[0]+DP4[1]+DP4[2]+DP4[3])<=150 & FILTER=="PASS"' \
pool1_annotated.vcf.gz -Oz -o pool1.f.vcf.gz
```

#### Final filters

