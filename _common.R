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

convert_to_entrez_id = function(x) {
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


head = function(x, ...) {
    utils::head(x, n = 4, ...)
}

# print.character = function(x, ...) {
#     w = getOption("width")
#     nm = names(x)
#     x = ifelse(nchar(x) > w, paste0(substr(x, 0, w-2),  ".."), x)
#     names(x) = nm
#     x
# }
