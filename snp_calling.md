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

#### Checking coverage
A table with coverage for each pool and contig length and script for plots
```
getCoverageFromVcf.py
plotCoverage.R
```
#### Annotation of bad quality snps
```
zcat snp_bcftools.vcf.gz | grep -v "^#" | awk '{print $1"\t"$2"\t"$2"\tsnp_"NR}' > snp_bcftools_annotations
bgzip snp_bcftools_annotations
tabix -s 1 -b 2 -e 3 snp_bcftools_annotations.gz
```
```
export PERL5LIB=/home/anna/github/vcftools/src/perl
bcftools view snp_bcftools.vcf.gz | vcf-annotate -f +/2=1e-100/3=1e-100/q=20/Q=20/-W/w=10 \
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


- excluding SNP with StrandBias (min P-value 0.0001)
- excluding SNP with BaseQualBias (min P-value 1e-100)
- excluding SNP with MapQualBias (min P-value 1e-100)
- excluding SNP with EndDistBias (min P-value 0.0001)
- min MapQ = 20
- min QUAL = 20
- no variant within 10bp from a gap
- only snps
- max coverage  = (c + 4*sqrt(c))*n , where c is average pool read depth and n is number of pools
- at least 4 reads supporting a variant (DP4{2}+DP4{3}>=4)
- at least 1 forward and 1 reverse read supporting a variant (DP4{2}>1 & DP4{3}>1)
- MAF = 0.01 (similar to min 4 reads) (DP4{2}+DP4{3})/DP4{0}+DP4{1}+DP4{2}+DP4{3})
- masking pools with DP < 10
- selecting only heterozygotes (?)

Keeping SNPs that PASS above filters + max DP = 705 + maf 0.01, min supporting variants = 4 (at least 1 forward and 1 reverse)
```
bcftools view -i 'QUAL>=20 & TYPE=="snp" & (DP4[0]+DP4[1]+DP4[2]+DP4[3])<=705 & (DP4[2]>=1) & (DP4[3]>=1) & (DP4[2]+DP4[3])>=4 & ((DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3]))>0.01 & FILTER=="PASS"' \
snp_bcftools_annotated.vcf.gz -Oz -o snp_bcftools_annotated.f.vcf.gz
```
Masking pools with DP < 10
```
vcftools --gzvcf snp_bcftools_annotated.f.vcf.gz --minDP 10 --recode --recode-INFO-all --out snp_bcftools_annotated.f2
bgzip snp_bcftools_annotated.f2.recode.vcf
vcftools --gzvcf snp_bcftools_annotated.f2.recode.vcf.gz --missing-indv --out snp_bcftools_annotated.f2.recode
vcftools --gzvcf snp_bcftools_annotated.f2.recode.vcf.gz --missing-site --out snp_bcftools_annotated.f2.recode

```

#### Plotting AF across pools

```
while read p
  do
  bcftools query -s $p -f '%CHROM\t%POS[\t%DP\t%AD{0}\t%AD{1}\t%GT\n]' $vcf -o snp_bcftools_annotated.f2.recode_${p}.tab
  done <pools.txt
```
Script for plotting AF for a sample of 10000 SNPs from each pool -> plotAF.R

