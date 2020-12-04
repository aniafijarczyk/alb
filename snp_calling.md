### SNP calling in pool-seq with bcftools

using samtools v1.9, bcftools v1.9, vcftools v0.1.17

#### Files

```
bams=*.bam
fasta=albgenome.fa
```

#### Read realignment
```
for bam in $bams
  do
  nameis=$(echo $bam | sed 's/.bam//g')
  echo $nameis
  samtools calmd -bAr -@2 ${bam} ${fasta} > ${nameis}_calmd.bam
  done
```

#### SNP calling
```
ls *_calmd.bam > bam_list
bcftools mpileup -C50 -f ${fasta} -q4 -a FORMAT/AD,FORMAT/ADR,FORMAT/ADF,FORMAT/DP -Ou -b bams_list | \
bcftools call -cv -f gq -Oz -o snp_bcftools.vcf.gz -
tabix -p vcf snp_bcftools.vcf.gz
```

#### Annotation of bad quality snps
```
zcat snp_bcftools.vcf.gz | grep -v "^#" | awk '{print $1"\t"$2"\t"$2"\tsnp_"NR}' > snp_bcftools_annotations
bgzip snp_bcftools_annotations
tabix -s 1 -b 2 -e 3 snp_bcftools_annotations.gz
```
```
export PERL5LIB=/home/anna/github/vcftools/src/perl
bcftools view snp_bcftools.vcf.gz | vcf-annotate -f +/d=15/q=20/Q=10/-W/w=10 \
--fill-HWE --fill-type -n -a snp_bcftools_annotations.gz -c CHROM,FROM,TO,INFO/SNP_ID \
-d key=INFO,ID=SNP_ID,Number=1,Type=Integer,Description='SnpList' | bgzip -c > snp_bcftools_annotated.vcf.gz
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
bcftools view -i 'TYPE=="snp" & FORMAT/AD[:1]>=4 & FORMAT/AD[:0]>=4 & FORMAT/DP[:]>=10 & \
(DP4[0]+DP4[1]+DP4[2]+DP4[3])>=40 & (DP4[0]+DP4[1]+DP4[2]+DP4[3])<=1000 & FILTER=="PASS"' \
snp_bcftools_annotated.vcf.gz -Oz -o snp_bcftools.f.vcf.gz
```
```
min depth per pool = 10
max depth across pools = 1000
min number of reads supporting a variant to call a snp = 4
```
