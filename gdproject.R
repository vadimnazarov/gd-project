parse.dd <- function (.path) {
  .parse.file <- function (.path) {
    suppressWarnings(df <- read.table(file = gzfile(.path), 
                                      header = T, 
                                      sep = "\t", 
                                      skip = 0, 
                                      strip.white = T, 
                                      comment.char = "", 
                                      quote = ""))
    names(df)[1:11] <- c("Read.count", 
                      "Read.proportion",
                      "CDR3.nucleotide.sequence",
                      "CDR3.amino.acid.sequence",
                      "V.gene", "D.gene", "J.gene",
                      "V.end", "D5.end", "D3.end", "J.start")
    df$Umi.count <- df$Read.count
    df$Umi.proportion <- df$Read.proportion
    
    logic <- df$V.end != -1 & df$D5.new1 != -1
    df$VD1.insertions <- -1
    df$VD1.insertions[logic] <- df$D5.new1 - df$VD1.insertions
    
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