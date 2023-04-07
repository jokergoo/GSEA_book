# GSEA in Genomics

## Overview

Genomics and epigenomics studies always generate may lists of genomic regions
of interest. For example, single-nucleotide variants (SNVs) or small insertion and deletions
(InDels) from whole genome sequencing or exon sequencing, peak regions of a certain
chromatin modifications or TFBS from ChIP-Sequencing data, or differetially methylated
regions fro whole genome bisulfite sequencing. Once the regions are obtained, the next
step is always functional intepretation. In this chapter, I will introduce the method
that specifically applies on genomic regions for functional enrichment analysis.


## Problem of ORA-based method

Since the functional intepretation is based on genes, a natural thought is to first annotate genomic
regions to genes, which generate a gene list, then to apply ORA on the gene list. The annotation to genes
can be done by looking for the nearest genes of every region. If the regions are believed to affect gene
transciprtion, they can also be annotated to the nearest gene TSS. In case some regions are too far from genes,
a cutoff of the distance from regions to the nearest genes can be set.

This conversion is very straightforward, however, let's recall the null assumption of ORA. In ORA, genes are treated as balls and all genes have the same probability and independently to be picked. However, in the mapping from regions to genes, the mapping from regions to genes are not balanced. Let's go back to the regions, actually, we assume, in the null hypotheiss, genomic regions are distrbuted uniformly in the genome, when converting to genes, due to somes regions locate in the gene desert and some reginos locates close to genes, this results in genes are not assigned wth the same probability to be picked. For exampe, if there is a clusters of regions close to a gene, then this genes is more likely to be picked. thus genes are with different weights. All these violates the null assumption of ORA. And such inproper comvertion will produce false positive results.


## The GREAT method

Recall the straighgy we proposed is to convert regions to genes, since it is not proper in the region to gene step, the tool GREAT proposes to directly constrcut region sets, which is a list of regions that associated to a biological term, then apply enrichment tests of input regions on the region set. This solves xxx. In the following sections, I will first intruduce how to construct such region sets and the enrichment tests.

<img src="figures/great_workflow.png">

### Construct functional genomic domains

For a gene set, GREAT directly constructs such “region sets” (or genomic domains) that associate with individual biological functions. In GREAT, there are three modes of linking regions to genes


<img src="06-in-genomics_files/figure-html/unnamed-chunk-1-1.png" width="576" style="display: block; margin: auto;" />

1. Basal plus extension

2. Two nearest genes

3. single nearest gene

### The binomial model



<img src="06-in-genomics_files/figure-html/unnamed-chunk-2-1.png" width="576" style="display: block; margin: auto;" />

The enrichment test is applied as follows. For a specific biological term represented as a gene set, denote the fraction of its associated functional domains in the genome (Fig. 1C) as p, the total number of input regions as N, the observed number of input regions that fall in the associated domains as n and the corresponding random variable as X, then X follows a binomial distribution: X ~ B(p, N) and the p-value of the enrichment is calculated as Pr(X ≥ n)

Note in GREAT, input regions are treated as single points and their middle points are used.


--|------|---------|
  |region |  gene
p | width(x)/widht(genome) | N_geneset/n_genome
N | N_region | n_DE gene
x | region in domain | DE gene in gene set


### implementation

We first extend the gene TSS. In the following code, we take human genome versio hg19 as an example.
In Bioconductor, there are also a list of "TxDb" packages which are standard Bioconductor
packages for transcripts information for a specific genome build. Here the package
**TxDb.Hsapiens.UCSC.hg19.knownGene** follows the naming skema as package_type.organism.source.genome.gene_id_type.

All the TxDb packages are built with the core facility package **GenomicFeatures**, which provide
functions to extract various transcipt level sof information. `genes()` extract the genomic coordinates of
all genes.


```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
g = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
g
# GRanges object with 23056 ranges and 1 metadata column:
#        seqnames            ranges strand |     gene_id
#           <Rle>         <IRanges>  <Rle> | <character>
#   1       chr19 58858172-58874214      - |           1
#   10       chr8 18248755-18258723      + |          10
#   ...       ...               ...    ... .         ...
#   9994     chr6 90539619-90584155      + |        9994
#   9997    chr22 50961997-50964905      - |        9997
#   -------
#   seqinfo: 93 sequences (1 circular) from hg19 genome
```

The gene variable `g` not only contains chromosomes, but also contains small contigs. This won't affect
the analysis too much, but to make the analysis clean. Let's restrict `g` to the chromosomes.

Note `g` also has a `gene_id` column which we will used to map to the gene set


```r
g = g[seqnames(g) %in% paste0("chr", c(1:22, "X", "Y"))]
```

I also reset the seqlevels


```r
all_chr = paste0("chr", c(1:22, "X", "Y"))
seqlevels(g) = all_chr
```

To get the positions of TSS which is a one-base region, `promoters()` regions by setting upstream extensions to zero
and downstream extension to 1. Note the gene variable `g` is stranded, so the direction of upstream
and downstream is automatically picked based on the strand in `g`.


```r
tss = promoters(g, upstream = 0, downstream = 1)
tss
# GRanges object with 23033 ranges and 1 metadata column:
#        seqnames    ranges strand |     gene_id
#           <Rle> <IRanges>  <Rle> | <character>
#   1       chr19  58874214      - |           1
#   10       chr8  18248755      + |          10
#   ...       ...       ...    ... .         ...
#   9994     chr6  90539619      + |        9994
#   9997    chr22  50964905      - |        9997
#   -------
#   seqinfo: 24 sequences from hg19 genome
```

To extend TSS is a two-step process. First construct the basal domains, which extend 5kb upstream and 1kb downstream.
THis can be donw also with `promoters()` on the  TSS.


```r
basal_domain = promoters(tss, upstream = 5000, downstream = 1000)
```

Construct the exteded domain is a little bit complex because we need to check whether the extension
reached neighbour gene's basal domains. SInce we need to compare to a gene's left and right genes, let's
first sort the coordinate of `basal_domain`. Also after obtained the basal domains, teh strandness is not
used anymore, here we set the strand of `basal_domains` to "\*" to mark the regions are not stranded.

```r
strand(basal_domain) = "*"
basal_domain = sort(basal_domain)
```

To illustrate the process of extension, and easily explain the idea, from here we assume `basal_domain`
only for one single chromosome. Note if you want to imeplement xx, teh next chunk of code should be applied 
to every chromsome.


The thought is we first extend basal domains by 1mb both upstream and downstream, next
we compare the neighbouring TSS, if it overlaps to the left gene's basal domain, the start position
the left basal domain is used as the start of the extension. Similar, if it overlaps to the right gene's basal
domain, the end positino ...

Step 1, we directly extend 1MB to the left and to the right. Note to make the code short, I assume the right extension
does not reach teh end of hte chromosome and the start extension does not reach the start of the chromosome.
In real case, you must additionally apply `left = pmax(left, 1)` and `right = pmin(right, chr_len)`.

```r
extension = 1e6
# from here we assume `basal_domain
left = start(basal_domain) - extension
right = end(basal_domain) + extension
```

To compare the left, we start form teh second region, and to compre the right, we start from the first region,
but only till the `n-1`th region.

```r
n = length(basal_domain)
left[2:n] = pmax(left[2:n], end(basal_domain)[1:(n-1)])
right[1:(n-1)] = pmin(right[1:(n-1)], start(basal_domain)[2:n])
```

Now `left` and `right` contains borders of xxx. since each region corrspond to the same
in `basal_domain`. We copy basal_domain, and directly change its start and end positions. 
Note the meta columns  in `basal_domains` is als kept in `extended_tss`.

```r
extended_tss = basal_domain
start(extended_tss) = left
end(extended_tss) = right
```

Step 2, construct the genomic domain that correspond to a gene set. Let's random select
100 genes as a gene set. SInce in `extended_tss` there is already a `gene_id` column, 
we can select the extened TSS for the genes in the gene set

```r
gs = sample(g$gene_id, 100)

region_set = extended_tss[extended_tss$gene_id %in% gs]
```

The extneded Tss may overlap, here it is necessary to apply `reduce()` to get a "flat"
version.

```r
region_set = reduce(region_set)
```

Now `region_set` is teh regions in teh genome that associate to a certain gene set.

Next we apply the binomial test with the region set.


```r
library(rGREAT)
gr = randomRegions(genome = "hg19")
```

In GREAT, only the middle point of the input regions are used. We first construct a `GRanges`
object which contains the middle points of regions.

```r
gr_mid = gr
start(gr_mid) = end(gr_mid) = mid(gr)
```

We overlap `gr_mid` to the region set. This can be done with the standard overlapping function
`findOverlaps()`. the object returned by `findOverlaps()` is a xx object. since each region in 
`region_set` does not overlap other regions and `gr_mid` only contains single base regions, 
the number of hits in the overlapping is the number of regions that fall in the genomic domain.

```r
ov = findOverlaps(gr_mid, region_set)

n_total = length(gr)
n_hits = length(ov)
```

To get the probabiliyt of xxx, we need to divide the totla length of xx to the total length of genome.
If the gene information is from teh TxDb package, the chromosome lengths information
is already in teh `seqlengths` xx and can be obtained with the `seqlengths()` function. Recall
we adjusted the seqlevels to teh normal chromsomes, teh chromosome lengths are adjusted automatically.
If you don't have such inforamtion, you have to look teh chromosme lengths for the
corect genome version. UCSC and NCBI genome are two useful resource to xx.

```r
prop = sum(width(region_set))/sum(seqlengths(extended_tss))
```

With all variables ready, we can simply apply `pbinom()` to calcualte the p-value.

```r
p = pbinom(n_hits - 1, n_total, prop, lower.tail = TRUE)
```

The implementation is a little bit complex, the complex part is to construct the region set.

### The hypergeometric model when dealing with background

GREAT also supports background. Here the "background" is more like a universal sets or a super sets of input regions. 
To provide the background, it has a special format. Assume the background regions is a set of 
regions [xi, yi] with index set A.., the input regions can only be a list of intervals (xj, yj) where
the the index set B is a subset of A. In this case...



<img src="06-in-genomics_files/figure-html/unnamed-chunk-9-1.png" width="576" style="display: block; margin: auto;" />

Such type of background/universal set is normally difficult to obtain. In Section x, I will introduce another
background model where background is not a unversal set while more like a mask on the genome where the enrichment
is only applied. I will also demonstrate the hypergenomic xx is a special case of the xx.

### The web-based tool

GREAT is implemented as a web-based tool. The use is very straightforward. 

### The rGREAT package

GREAT web service is useful for xxx, but xxx. The rGREAT package can make
GREAT analysis automatic by automaticlaly submitting genomic regions to GREAT
web service and retrieve results from there on the fly. The several functions
introduced later is to silmulate the process of xx teh website.

I used a TFBS bind sites of the TFBS (JUN) from UCSC database.... I first read it into a data frame
where the first three columns are the chromsomes, the start positinos oand the end 
positions of every peak regions of JUN.


```r
bed = read.table("data/tb_encTfChipPkENCFF708LCH_A549_JUN_hg19.bed")
head(bed)
#     V1      V2      V3 V4   V5 V6
# 1 chr1  936221  936347  . 1000  .
# 2 chr1 1280252 1280492  .  480  .
# 3 chr1 2143913 2144153  .  982  .
# 4 chr1 4052186 4052426  .  544  .
# 5 chr1 4440698 4440938  .  837  .
# 6 chr1 7359438 7359678  .  493  .
```

#### Submit regions

Just like the web interface, the only mandatory argumnet is the input regions. The input regions can
also be in teh `GRanges` object. The two foramts are baically give identical information.

`submitGreatJob()` automatically construct all necessary parameters and submit to the GREAT web service. Just like
you clicked "submit" button in the web interface. After the job is successfually submitted, the returned object
`job` contains the settings of the xxx.


```r
job = submitGreatJob(bed)
job
# Submit time: 2023-04-07 09:34:07 
#   Note the results may only be avaiable on GREAT server for 24 hours.
# Version: 4.0.4 
# Species: hg19 
# Inputs: 1726 regions
# Mode: Basal plus extension 
#   Proximal: 5 kb upstream, 1 kb downstream,
#   plus Distal: up to 1000 kb
# Include curated regulatory domains
# 
# Enrichment tables for following ontologies have been downloaded:
#   None
```

A specific genome version for an organism should be set. By default it is hg19 for human genome.
GREAT has several versions. The newest version is 4.0.4 released on 2019. It supports both hg19 and hg38. Other genomoes
includes... Note it is very important to select teh correct genome version for the organism. In Section
xx I will demonstrate the consequence of selecting the wrong genome version.

The genome version can be set via the `species` argumetn

```r
submitGreatJob(bed, species = "hg38")  # assume bed is from hg38
submitGreatJob(bed, species = "mm10")  # assume bed is from mm10
```

`submitGreatJob()` supports all historical versions of GREAT. Normally we use the newest version since the annotation
is the newest. But when GREAT updates from version 3 to 4, many gene set collections are removed. So, if you want to
use more gene set collections, you may need to swtich to older verison. This can be set via the `version` argument

```r
submitGreatJob(bed, version = "3.0")
```

If by any chance, you also have the "background" universe set of the regions, you can set it with the `bg` argument.
But remeber the format of the two xxx.

```r
bed_subset = bed[sample(nrow(bed), 500), ]
submitGreatJob(bed_subset, bg = bed)
```

#### obtain enrichment tables

Once the regions are successfully submitted. GREAT performs the analysis and temporary results will be on the
GREAT website for some time, which other functions can automaticaly retrieve the result from there.
GERAT applies analysis on all supported gene set collections. The function `getEnrichmentTables()` downloads
the complete result table for specific gene set collections.

Calling `getEnrichmentTables()` without argumetn downloads the enrichment tables for the three GO collectios, i.e.
BP, MF and CC.


```r
tbl = getEnrichmentTables(job)
names(tbl)
# [1] "GO Molecular Function" "GO Biological Process" "GO Cellular Component"
```

Once the enrichment table is downloaded, it is cached. Repeatedly call `getEnrichmentTables()` will directly use
the table from cache and will not download it again. This is good because the results actally are temporary saved
on the GREAT server. Once the tables are cached in the `job` object. It can be xxx

`tbl` is just a list of data frames, where each data frmae is xxx which is identical to the result if you submit
the regions directly from GREAT web service (of course it is because all calculations are performed on GREAT).



```r
head(tbl[[1]])
#           ID                           name Binom_Genome_Fraction
# 1 GO:0005515                protein binding            0.75150510
# 2 GO:0005488                        binding            0.86944070
# 3 GO:0019899                 enzyme binding            0.17198290
# 4 GO:0019900                 kinase binding            0.06564039
# 5 GO:0045296               cadherin binding            0.02779121
# 6 GO:0050839 cell adhesion molecule binding            0.04901075
#   Binom_Expected Binom_Observed_Region_Hits Binom_Fold_Enrichment
# 1     1297.09800                       1500              1.156428
# 2     1500.65500                       1628              1.084860
# 3      296.84260                        460              1.549643
# 4      113.29530                        206              1.818257
# 5       47.96764                        113              2.355755
# 6       84.59255                        162              1.915062
#   Binom_Region_Set_Coverage Binom_Raw_PValue Binom_Adjp_BH Hyper_Total_Genes
# 1                0.86906140     8.378648e-34  3.534952e-30             11129
# 2                0.94322130     4.881720e-24  1.029799e-20             14268
# 3                0.26651220     7.001972e-23  9.847107e-20              1885
# 4                0.11935110     2.376577e-16  2.378356e-13               647
# 5                0.06546929     2.818625e-16  2.378356e-13               304
# 6                0.09385863     8.118197e-15  5.708446e-12               461
#   Hyper_Expected Hyper_Observed_Gene_Hits Hyper_Fold_Enrichment
# 1     1360.75100                     1579              1.160389
# 2     1744.55900                     1873              1.073624
# 3      230.48030                      346              1.501213
# 4       79.10917                      134              1.693862
# 5       37.17031                       79              2.125352
# 6       56.36681                      114              2.022467
#   Hyper_Gene_Set_Coverage Hyper_Term_Gene_Coverage Hyper_Raw_PValue
# 1              0.69620810                0.1418816     2.492149e-24
# 2              0.82583770                0.1312728     1.339393e-12
# 3              0.15255730                0.1835544     2.887885e-16
# 4              0.05908289                0.2071097     3.499515e-10
# 5              0.03483245                0.2598684     3.394731e-11
# 6              0.05026455                0.2472885     6.033715e-14
#   Hyper_Adjp_BH
# 1  1.051438e-20
# 2  1.412725e-09
# 3  6.091993e-13
# 4  2.109208e-07
# 5  2.864474e-08
# 6  8.485415e-11
```


The columns on the enrichment are self-explorary. The columns are:

- ID
- name
- Binom_Genome_Fraction
...


All supported gene set collections can be obtained by `availableOntologies(job)`:


```r
availableOntologies(job)
# [1] "GO Molecular Function"     "GO Biological Process"    
# [3] "GO Cellular Component"     "Mouse Phenotype"          
# [5] "Mouse Phenotype Single KO" "Human Phenotype"          
# [7] "Ensembl Genes"
```

and the ontologies can be set in `getEnrichmentTables()`


```r
getEnrichmentTables(job, ontology = c("GO Biological Process", "Ensembl Genes"))
```


#### get region-gene associations

In the GREAT process, each extened TSS is associated to a gene. For a list of input regions,
we can obtain the associations between regions and genes, i.e. whether the region (middle point)
falls in the extended TSS of a gene. The function `getRegionGeneAssociations()` returns an
`GRanges` object of the associations between regions and genes. Note teh number of regions
in the object may be less than the total number of input xxx because some regions are locate
in teh gene desert...

In teh meta column of xxx. You may see the column of `annotated_genes` and `dist_to_TSS` are in
`Characterlist` and `IntegerList`...


```r
getRegionGeneAssociations(job)
# GRanges object with 1720 ranges and 2 metadata columns:
#          seqnames              ranges strand | annotated_genes    dist_to_TSS
#             <Rle>           <IRanges>  <Rle> | <CharacterList>  <IntegerList>
#   [1]        chr1       936221-936347      * |            HES4           -732
#   [2]        chr1     1280252-1280492      * |     DVL1,TAS1R3     4120,13678
#   ...         ...                 ...    ... .             ...            ...
#   [1719]     chrX 149369317-149369557      * |   MAMLD1,MAGEA8 -162249,359496
#   [1720]     chrX 153237211-153237451      * |   TMEM187,HCFC1       -447,-73
#   -------
#   seqinfo: 23 sequences from an unspecified genome; no seqlengths
```


The distributions can b


```r
plotRegionGeneAssociations(job)
```

<img src="06-in-genomics_files/figure-html/unnamed-chunk-16-1.png" width="1152" style="display: block; margin: auto;" />

THe associations can also be obtains for a specific gene set, just specigy the gene set collection
and the corresponding gene set ID:

```r
getRegionGeneAssociations(job, ontology = "GO Molecular Function",
    term_id = "GO:0004984")
```


#### visualizations and reports


Similar as differentiall expression analysis, the enrichment table can also be visualized as a
volcano plot. Recall in the DE analysis, in the volcano, on x-axis is teh log2 fold change
and on y-axis is the log10(p-value) or log10(FDR). Similar, we can calcualte the log2 enrichment rate as

$$log2(obs/exp) $$

where `obs` is the observed number of regions in the region set, and exp is calcualted as $p*n$. 
The expected values is already in the enrichment table.

Function `plotVolcano()` makes a volcano plot for a specific gene set collection. In
the following example, we use hte "GO Biologucal Process".


```r
plotVolcano(job, ontology = "GO Biological Process")
```

<img src="06-in-genomics_files/figure-html/unnamed-chunk-17-1.png" width="672" style="display: block; margin: auto;" />

THe volcano plot is actually "one-sided" because we only look at the over-representation. Similar
as the differnetial expression analysis, we can set two cutoffs, one for the log2 enrichment ratio,
and one for the log10 pvalue, i.e. one from biologucal aspect and one from statstical aspect. 

In teh plot, we can see it is composed by several curves. Each curve correspond to a same region hits.
we also see xxx which we will explain in more detail in Section x.

One last function to introuce is the `shinyReport()` function which can create xxx

```r
shinyReport(job)
```

## Local GREAT analysis

As an online tool, all annotation resources are only controlled by GREAT
developers; they are not extensible by users. The current version (4.0.4) of
GREAT only supports human and mouse, and it only supports seven gene set
collections which have not yet been updated to the most recent ones.

Bioconductor has a rich and up-to-date annotation resource, also as I have
demonstrated in sectionx, impleement the GREAT algorithm from scratch with
Bioconductor facility package is not difficult, thus, in the **rGREAT**
package, I also implement a "local GREAT" function, which greatly ...

On one hand, with local GREAT, it is possible to perform enrichment analysis
on any organism and with any type of gene set collection; and on the other
hand, Bioconductor annotation packages are well maintained and updated, which
ensures that a local GREAT analysis always uses the most up-to-date annotation
data. Local GREAT by default supports many gene set collections and more than
600 organisms, and, more importantly, local GREAT allows self-provided gene
sets and organisms from users.

To demonstrate its use, I first convert the object `bed` to a `GRanges` object
because ...


```r
gr = GRanges(seqnames = bed[, 1], ranges = IRanges(bed[, 2], bed[, 3]))
gr
# GRanges object with 1726 ranges and 0 metadata columns:
#          seqnames            ranges strand
#             <Rle>         <IRanges>  <Rle>
#   [1]        chr1     936221-936347      *
#   [2]        chr1   1280252-1280492      *
#   ...         ...               ...    ...
#   [1725]    chr22 43011250-43011490      *
#   [1726]    chr22 50364455-50364695      *
#   -------
#   seqinfo: 23 sequences from an unspecified genome; no seqlengths
```

The core function `great()` performs GREAT analysis locally. It used gene set and gene 
all from Biocondutor packages. The first three argumetns are teh input regions,
the name of teh gene set collection, and the name/version of hte genome.

In teh following example, the regions correspond to human genome hg19 and we use GO BP gene sets.
The returned object prints parametes.



```r
res = great(gr, "GO:BP", "hg19")
res
# 1726 regions are associated to 2268 genes' extended TSSs.
#   TSS source: TxDb.Hsapiens.UCSC.hg19.knownGene
#   Genome: hg19
#   OrgDb: org.Hs.eg.db
#   Gene sets: GO:BP
#   Background: whole genome excluding gaps
# Mode: Basal plus extension
#   Proximal: 5000 bp upstream, 1000 bp downstream,
#   plus Distal: up to 1000000 bp
```

### use predefined-gene sets and TSS sources

### use self-defined extended TSS

### Use self-defined gene sets

### Use both self-defined gene sets and TSS

### The binomial background model



```r

p5 = grid.grabExpr({
grid.newpage()
grid.lines(c(0.02, 0.98), c(0.5, 0.5))
grid.text("A", x = unit(0.02, "npc"), y = unit(0.5, "npc") + unit(4, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 14))

grid.segments(pos[2] - 0.05/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[5] - 0.05/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))
grid.segments(pos[6] - 0.01/2 - 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[7] + 0.01/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))
grid.segments(pos[9] - 0.05/2 - 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[10] + 0.01/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))
grid.segments(pos[12] - 0.05/2 - 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[12] + 0.01/2 + 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))

grid.rect(0.13, 0.5, width = 0.12, height = unit(1, "cm"), gp = gpar(fill = "#FF000020", col = "#FF000080"))
grid.rect(0.36, 0.5, width = 0.28, height = unit(1, "cm"), gp = gpar(fill = "#FF000020", col = "#FF000080"))
grid.rect(0.68, 0.5, width = 0.12, height = unit(1, "cm"), gp = gpar(fill = "#FF000020", col = "#FF000080"))

grid.text(TeX("$\\textit{p}_2 = 0.619$"), x = 0.9, y = unit(0.5, "npc") + unit(3, "mm"), just = "bottom", gp = gpar(fontsize = 10))
})


p6 = grid.grabExpr({
grid.newpage()
grid.lines(c(0.02, 0.98), c(0.5, 0.5))
grid.text("B", x = unit(0.02, "npc"), y = unit(0.5, "npc") + unit(4, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 14))

grid.segments(pos[2] - 0.05/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[5] - 0.05/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))
grid.segments(pos[6] - 0.01/2 - 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[7] + 0.01/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))
grid.segments(pos[9] - 0.05/2 - 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[10] + 0.01/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))
grid.segments(pos[12] - 0.05/2 - 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), pos[12] + 0.01/2 + 0.1/2, unit(0.5, "npc") - unit(1.5, "mm"), gp = gpar(col = "#4DAF4A", lwd = 2))


grid.points(x, unit(rep(0.5, 20), "npc") + unit(2, "mm"), default.unit = "npc", pch = 25, gp = gpar(col = "blue", fill = "blue", cex = 0.6))

grid.rect(0.13, 0.5, width = 0.12, height = unit(1, "cm"), gp = gpar(fill = "#FF000020", col = "#FF000080"))
grid.rect(0.36, 0.5, width = 0.28, height = unit(1, "cm"), gp = gpar(fill = "#FF000020", col = "#FF000080"))
grid.rect(0.68, 0.5, width = 0.12, height = unit(1, "cm"), gp = gpar(fill = "#FF000020", col = "#FF000080"))
grid.text(TeX("$Pr_{binom}(\\textit{X}_2 \\geq 8 | \\textit{N}_2 = 15, \\textit{p}_2 = 0.619)$"), x = 0.75, y = unit(0.5, "npc") + unit(5.5, "mm"), just = "bottom", gp = gpar(fontsize = 10))

})

library(ComplexHeatmap)
lgd = Legend(labels = c("Genomic domains associated to a gene set", "Input regions", "Background regions"),
    graphics = list(
        function(x, y, w, h) {
        	grid.segments(x-w*0.5, y, x+w*0.5, y, gp = gpar(col = "#4DAF4A", lwd = 2))
        },
        function(x, y, w, h) {
        	grid.points(x, y, pch = 25, gp = gpar(col = "blue", fill = "blue", cex = 0.6))
        },
        function(x, y, w, h) {
        	grid.rect(x, y, w, h, gp = gpar(fill = "#FF000020", col = "#FF000080"))
        }
    ), nrow = 1, gap = unit(5, "mm")
)


library(cowplot)
lgd@grob$vp$y = unit(1, "npc")

print(
	plot_grid(p5, p6, lgd@grob, ncol = 1, scale = 0.95)
)
```

<img src="06-in-genomics_files/figure-html/unnamed-chunk-20-1.png" width="672" style="display: block; margin: auto;" />

### with more organisms

### Genome version is important

## Other thoughts

### large sample size

### compare between different sets of genomic regions

Distance to TSS

### Use other statistics and permutate

Once we have a region set They have a same problem where larger xxx tend to give more significant p-values.

### use covariates

## Summary


