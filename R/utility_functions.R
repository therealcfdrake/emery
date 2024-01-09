
name_thing <-
  function(thing = "", n = 1){

    paste0(thing, stringr::str_pad(1:n, width = nchar(n), side = "left", pad = "0"))

  }

define_disease_state <-
  function(D = NULL, n_obs = NULL, prev = NULL){

    if(is.null(D)){
      if(!is.numeric(prev)|!is.numeric(n_obs)){
        stop("Invalid disease state specified", call. = FALSE)
      }else{
        pos <- round(n_obs * prev, 0)
        neg <- n_obs - pos
        D <- c(rep(1, pos), rep(0, neg))
      }
    }else{
      n_obs <- length(D)
      pos <- sum(D)
      neg <- sum(1 - D)
      prev <- mean(D)
    }

    return(
      list(
        D = D,
        n_obs = n_obs,
        prev = prev,
        pos = pos,
        neg = neg
      )
    )

  }


