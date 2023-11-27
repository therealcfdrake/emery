
name_thing <- function(thing = "", n = 1){

  paste0(thing, stringr::str_pad(1:n, width = nchar(n), side = "left", pad = "0"))

}
