parse.dd <- function (.path) {
  .parse.file <- function (.path) {
    suppressWarnings(df <- read.table(file = gzfile(.path), 
                                      header = T, 
                                      sep = "\t", 
                                      skip = 0, 
                                      strip.white = T, 
                                      comment.char = "", 
                                      quote = ""))
    df
  }
  
  res <- list()
  if (dir.exists(.path)) {
    res <- lapply(list.files(.path, full.names = T), .parse.file)
    names(res) <- list.files(.path)
  } else if (file.exists(.path)) {
    res <- .parse.file(.path)
  } else {
    cat('Can\'t find folder or file:\t"', .path, '"', sep = '', end = '\n')
  }
  res
}