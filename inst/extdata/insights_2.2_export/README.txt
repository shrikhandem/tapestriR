This folder contains data exported from the Tapestri Insights 2.2 application.
The following files are included:

AF.csv - variant allele frequency
         (# of reads with evidence for mutation / total # of reads * 100)

DP.csv - read depth
    The Read Depth per Variant metric is the filtered depth, at the cell level. This gives the number
    of filtered reads that support each of the reported alleles.

GQ.csv - quality scores
    The genotype quality score represents the Phred-scaled confidence that the genotype assignment (GT)
    is correct, derived from the genotype normalized likelihoods of possible genotypes (PL). Specifically,
    the quality score is the difference between the PL of the second most likely genotype and the PL
    of the most likely genotype. The values of the PLs are normalized so that the most likely PL is always 0.
    So, the quality score ends up being equal to the second smallest PL, unless that PL is greater than 99.
    
    In GATK, the value of quality score is capped at 99 because larger values are not more informative,
    but they take up more space in the file. So if the second most likely PL is greater than 99, we still
    assign a quality score of 99.
    
    Basically, the GQ gives you the difference between the likelihoods of the two most likely genotypes.
    If it is low, you can tell there is not much confidence in the genotype, i.e., there was not enough
    evidence to confidently choose one genotype over another.

NGT.csv - genotype call converted to categorical
          (0-reference, 1-heterozygous mutation, 2-homozygous mutation, 3-unknown)
          Calls are made by GATK/HaplotypeCaller.

Variants.csv - variant metadata
