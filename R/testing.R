
#### Convert ACC data to .rda format using the save function inside the "data" folder.
#### Then simply load the data back from that same folder.



fred <- function(x) UseMethod("fred")
fred(1:4)


fred.default <- function(x) {
  cat("in fred.default\n")
  invisible(NULL)
}
fred(1:4)
