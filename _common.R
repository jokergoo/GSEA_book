# example R options set globally
options(width = 60)

# example chunk options set globally
knitr::opts_chunk$set(
  comment = "#",
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.pos = "htp"
  # out.width = "1\\linewidth"
)

options("width" = 80)

convert_to_entrez_id = function(x, from) {
  if(is.matrix(x)) {
    map = cola:::guess_id_mapping(rownames(x))
    if(is.function(map)) {
        new_rn = map(rownames(x))
    } else {
        new_rn = map[rownames(x)]
    }
    l = is.na(new_rn)

    x = x[!l, , drop = FALSE]
    new_rn = new_rn[!l]

    x2 = do.call(rbind, tapply(1:nrow(x), new_rn, function(ind) {
      colMeans(x[ind, , drop = FALSE])
    }))
    return(x2)

  } else if(is.numeric(x)) {
    map = cola:::guess_id_mapping(names(x))
    x2 = s
    if(is.function(map)) {
        names(x2) = map(names(x))
    } else {
        names(x2) = map[names(x)]
    }
    x2 = x2[!is.na(names(x2))]
    x2 = tapply(x2, names(x2), mean)
    return(x2)
  } else {
      map = cola:::guess_id_mapping(x)
      if(is.function(map)) {
          x2 = map(x)
      } else {
          x2 = map[x]
      }
      x2 = x2[!is.na(x2)]
      x2 = unique(x2)
      return(x2)
  }
}


format.character = function(x, ...) {
  w = getOption("width")
  ifelse(nchar(x) > w, paste0(substr(x, 0, w-2),  ".."), x)
}

str = function(object, ...) {
  utils::str(object, strict.width = "cut", ...)
}


assignInNamespace("get_showHeadLines", ns = "S4Vectors",
    value = function() S4Vectors:::.get_showLines(2L, "showHeadLines"))
assignInNamespace("get_showTailLines", ns = "S4Vectors",
    value = function() S4Vectors:::.get_showLines(2L, "showTailLines"))


head.enrichResult = function(x, ...) {
    tb = x@result
    rownames(tb) = NULL
    tb$ID = ifelse(nchar(tb$ID) > 20, paste0(substr(tb$ID, 0, 20),  ".."), tb$ID)
    tb$Description = ifelse(nchar(tb$Description) > 20, paste0(substr(tb$Description, 0, 20),  ".."), tb$Description)
    tb$geneID = ifelse(nchar(tb$geneID) > 20, paste0(substr(tb$geneID, 0, 20),  ".."), tb$geneID)
    print(utils::head(tb, n = 4, ...))
}


gs_list2dataframe = function(list) {
    data.frame(gene_set = rep(names(list), times = sapply(list, length)),
               gene = unlist(list))
}


gs_dataframe2list = function(df) {
    n1 = length(unique(df[, 1]))
    n2 = length(unique(df[, 2]))
    if(n1 < n2) {
        split(df[, 2], df[, 1])
    } else {
        split(df[, 1], df[, 2])
    }
}

# print.character = function(x, ...) {
#     w = getOption("width")
#     nm = names(x)
#     x = ifelse(nchar(x) > w, paste0(substr(x, 0, w-2),  ".."), x)
#     names(x) = nm
#     x
# }

msigdb_base_url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
msigdb_env = new.env(parent = emptyenv())
msigdb_env$all_versions = NULL
msigdb_env$file_list = list()
msigdb_env$all_collections = list()
msigdb_env$gene_sets = list()


list_msigdb_versions = function() {
    if(is.null(msigdb_env$all_versions)) {
        all_versions = get_file_list(msigdb_base_url)
        msigdb_env$all_versions = all_versions
    } else {
        all_versions = msigdb_env$all_versions
    }
    all_versions
}

choose_msigdb_version = function() {
    all_versions = list_msigdb_versions()
    ind = menu(all_versions, title = "Choose an MSigDB version:")
    all_versions[ind]
}

list_msigdb_collections = function(version = NULL) {
    if(is.null(version)) {
        version = choose_msigdb_version()
    }

    if(is.null(msigdb_env$file_list[[version]])) {
        files = get_file_list(paste0(msigdb_base_url, "/", version))
        files = files[grep("\\.gmt$", files)]
        files = files[!grepl("^msigdb", files)]
        msigdb_env$file_list[[version]] = files
    } else {
        files = msigdb_env$file_list[[version]]
    }
    collections = unique(gsub(paste0(".v", version, ".*$"), "", files))
    msigdb_env$all_collections[[version]] = collections
    collections
}

choose_msigdb_collection = function(version) {
    
    collections = list_msigdb_collections(version)
    ind = menu(collections, title = paste0("Choose a gene set collection for version ", version, ":"))
    collections[ind]
}

get_file_list = function(url) {
    con = url(url)
    on.exit(close(con))
    ln = readLines(con)
    ind = grep("^<tr><td", ln)
    rows = ln[ind]
    rows = gsub("</td><td[^>]*>", ";", rows)
    rows = gsub("<.*?>", "", rows)
    files = sapply(strsplit(rows, ";"), function(x) x[2])[-1]
    gsub("/", "", files)
}


get_msigdb = function(version = choose_msigdb_version(), 
    collection = choose_msigdb_collection(version),
    gene_id_type = c("entrez", "symbols"), as_table = FALSE) {
    
    version = force(version)
    gene_id_type = match.arg(gene_id_type)

    if(is.null(msigdb_env$all_versions)) {
        list_msigdb_versions()
    }
    
    if(!version %in% msigdb_env$all_versions) {
        i = grep(version, msigdb_env$all_versions, ignore.case = TRUE)
        if(length(i) == 1) {
            version = msigdb_env$all_versions[i]
        } else {
            message(paste0("Cannot find version '", version, "', please select a valid value."))
            version = choose_msigdb_version()
        }
    }

    collection = force(collection)
    if(is.null(msigdb_env$all_collections[[version]])) {
        list_msigdb_collections(version)
    }
    if(!collection %in% msigdb_env$all_collections[[version]]) {
        i = grep(collection, msigdb_env$all_collections[[version]], ignore.case = TRUE)
        if(length(i) == 1) {
            collection = msigdb_env$all_collections[[version]][i]
        } else {
            message(paste0("Cannot find collection '", collection, "', please select a valid value."))
            collection = choose_msigdb_collection(version)
        }
    }
    
    url = paste0(msigdb_base_url, "/", version, "/", collection, ".v", version, ".", gene_id_type, ".gmt")
    basename = basename(url)
    if(is.null(msigdb_env$gene_sets[[basename]])) {
        con = url(url)
        on.exit(close(con))
        ln = readLines(con)
        ln = strsplit(ln, "\t")
        gs = lapply(ln, function(x) x[-(1:2)])
        names(gs) = sapply(ln, function(x) x[1])
        msigdb_env$gene_sets[[basename]] = gs
    } else {
        gs = msigdb_env$gene_sets[[basename]]
    }
   
    if(as_table) {
        df = data.frame(gene_set = rep(names(gs), times = sapply(gs, length)),
                        gene = unlist(gs))
        rownames(df) = NULL
        df
    } else {
        gs
    }
}


library(Rcpp)
sourceCpp(code = '
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
List intersectList(List input, StringVector vec) {

  int n = input.size();
  List out(n);

  std::unordered_set<String> seen;
  seen.insert(vec.begin(), vec.end());

  for (int i = 0; i < n; i++) {
      
      StringVector sp = as< StringVector > (input[i]);
      LogicalVector l(sp.size());

      std::unordered_set<String> seen2;

      for(int j = 0; j < sp.size(); j ++) {
        l[j] = seen.find(sp[j]) != seen.end() && seen2.insert(sp[j]).second;
      }

      out[i] = sp[l];
  }

  return out;
}
')


ora = function(genes, gene_sets, universe) {
    # 1
    if(is.data.frame(gene_sets)) {
        gene_sets = gs_dataframe2list(gene_sets)
    }
    # 2
    if(missing(universe)) {
        universe = unique(unlist(gene_sets))
    } else {
        universe = unique(universe)
    }
    # 3
    gs_names = names(gene_sets)
    genes = intersect(genes, universe)
    # gene_sets = lapply(gene_sets, function(x) intersect(x, universe))
    gene_sets = intersectList(gene_sets, universe)

    n_universe = length(universe)
    n_genes = length(genes)
    #4
    # x = sapply(gene_sets, function(x) length(intersect(x, genes)))
    x = sapply(intersectList(gene_sets, genes), length)
    m = sapply(gene_sets, length)
    n = n_universe - m
    k = n_genes
    p = phyper(x - 1, m, n, k, lower.tail = FALSE)
    # 5
    data.frame(gene_set = gs_names, 
               hits = x,
               gene_set_size = m,
               ratio_in_gene_set = x/m,
               ratio_in_genes = x/n_genes,
               enrichment_score = x*n_universe/m/n_genes,
               p_value  = p,
               p_adjust = p.adjust(p, "BH"))
}


get_go_gene_sets_from_orgdb = function(orgdb, ontology = "BP") {
    if(any(!ontology %in% c("BP", "CC", "MF"))) {
        stop("`ontology` can only be one of 'BP', 'CC' and 'MF'.")
    }
    suppressMessages(tb <- select(orgdb, keys = ontology,
        keytype = "ONTOLOGYALL", columns = c("GOALL", "ENTREZID")))
    gs = split(tb$ENTREZID, tb$GOALL)
    lapply(gs, unique)
}

get_all_pc_genes_from_orgdb = function(orgdb) {
    suppressMessages(tb <- select(orgdb, keys = "protein-coding", keytype = "GENETYPE", columns = "ENTREZID"))
    tb$ENTREZID
}
