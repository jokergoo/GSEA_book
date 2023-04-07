# Clustering and simplifying GSEA results


## Overview

Functional enrichment is becoming a standard analysis in xxx

## measures of similarities

### overlap-based 



Show which terms tend to show high similarity, parent-child, or siblinds

### semantic measures

#### IC


```r
go_id_1 = "GO:0032888"
go_id_2 = "GO:0051256"
```


```r
library(org.Hs.eg.db)
tb = toTable(org.Hs.egGO2ALLEGS)
tb = tb[tb$Ontology == "BP", , drop = FALSE]
go_gene_sets = split(tb$gene_id, tb$go_id)

library(GO.db)
offspring = as.list(GOBPOFFSPRING)
```


```r
n = sum(sapply(c(go_id_1, offspring[[go_id_1]]), function(id) {
	length(go_gene_sets[[id]])
}), na.rm = TRUE)
n_bg = sum(sapply(go_gene_sets, length))
-log(n/n_bg)
# [1] 13.91922
```


```r
calc_IC = function(go_id) {
	n = sum(sapply(c(go_id, offspring[[go_id]]), function(id) {
		length(go_gene_sets[[id]])
	}), na.rm = TRUE)
	n_bg = sum(sapply(go_gene_sets, length))
	-log(n/n_bg)
}
IC_1 = calc_IC(go_id_1)
IC_2 = calc_IC(go_id_2)
```



```r
CA = intersect(GOBPANCESTOR[[go_id_1]], GOBPANCESTOR[[go_id_2]])
CA = setdiff(CA, "all")
```


```r
IC_MICA = max(sapply(CA, calc_IC))
```


```r
2*IC_MICA/(IC_1 + IC_2)
# [1] 0.8660491
```

## Enrichment map

## David

## simplifyEnrichment
