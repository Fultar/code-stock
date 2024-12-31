# 1. plink过滤低质量snps
```
## vcftools过滤低质量SNPs
vcftools --vcf /home/yzhou/big/vcf_100_sample/Coilia_ZY_Big_Small_all_samples_raw_filtered_snp_hard.vcf  --minQ 30  --minDP 3  --min-meanDP 3 --recode --recode-INFO-all --out Coilia_ZY_Big_Small_all_samples_vcftoolsfilted
#After filtering, kept 2468370 out of a possible 5061188 Sites

## 过滤SNP缺失率高于5%的SNP
plink --vcf Coilia_ZY_Big_Small_all_samples_vcftoolsfilted.recode.vcf --geno 0.05 --make-bed --out Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno --chr-set 24 --double-id
#4109465 variants removed due to missing genotype data (--geno).
#317703 variants and 100 samples pass filters and QC.

## maf过滤
plink --bfile Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno --maf 0.01 --make-bed --out Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf --chr-set 24 --double-id
#208286 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#109417 variants and 100 samples pass filters and QC.

## 哈温平衡过滤
plink --bfile Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf --hwe 1e-4 --make-bed --out Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf_hwe --chr-set 24 --double-id
#--hwe: 385 variants removed due to Hardy-Weinberg exact test.#23811 variants and 100 samples pass filters and QC.
#109032 variants and 100 samples pass filters and QC.

```

# 2.模型构建
```
#### 体长
#输出亲缘矩阵
gemma-0.98.5 -bfile ../Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf_hwe -p coilia_genotype_length.gemma.txt -gk 2 -o length_kinship

#LMM模型
gemma-0.98.5 -bfile ../Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf_hwe -k output/length_kinship.sXX.txt -lmm 1 -p coilia_genotype_length.gemma.txt -c corvariates_weight_pca.gemma.txt -o Coilia_length_LMM_corv


#### 体重
gemma-0.98.5 -bfile ../Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf_hwe -p coilia_genotype_weight.gemma.txt -gk 2 -o weight_kinship

gemma-0.98.5 -bfile ../Coilia_ZY_Big_Small_all_samples_vcftoolsfilted_geno_maf_hwe -k output/weight_kinship.sXX.txt -lmm 1 -p coilia_genotype_weight.gemma.txt -c corvariates_length_pca.gemma.txt -o Coilia_weight_LMM_corv1

```
