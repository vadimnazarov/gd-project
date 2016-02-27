parse.dd <- function (.path) {
  .parse.file <- function (.path) {
    print(.path)
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
    
    .add.one <- function(.df, .col) {
      .df[[.col]][.df[[.col]] != -1] <- .df[[.col]][.df[[.col]] != -1] + 1
      .df
    }
    df <- .add.one(df, "V.end")
    df <- .add.one(df, "D5.new1")
    df <- .add.one(df, "D5.new2")
    df <- .add.one(df, "D5.new3")
    df <- .add.one(df, "D3.new1")
    df <- .add.one(df, "D3.new2")
    df <- .add.one(df, "D3.new3")
    df <- .add.one(df, "J.start")
    
    df$Total.insertions <- -1
    for (i in 1:nrow(df)) {
      acc <- 0
      if (df$D5.new1[i] != -1) {
        if (df$V.end[i] >= df$D5.new1[i]) {
          acc <- -11
        } else {
          acc <- df$D5.new1[i] - df$V.end[i] - 1
          
          if (df$D5.new2[i] != -1) {
            if (df$D3.new1[i] >= df$D5.new2[i]) {
              acc <- -12
            } else {
              acc <- acc + df$D5.new2[i] - df$D3.new1[i] - 1
              
              if (df$D5.new3[i] != -1) {
                if (df$D3.new2[i] >= df$D5.new3[i] || df$J.start[i] <= df$D3.new3[i]) {
                  acc <- -13
                } else {
                  acc <- acc + df$D5.new3[i] - df$D3.new2[i] - 1 + df$J.start[i] - df$D3.new3[i] - 1
                }
              } else {
                if (df$J.start[i] <= df$D3.new2[i]) {
                  acc <- -14
                } else {
                  acc <- acc + df$J.start[i] - df$D3.new2[i] - 1
                }
              }
            }
          } else {
            if (df$J.start[i] <= df$D3.new1[i]) {
              acc <- -15
            } else {
              acc <- acc + df$J.start[i] - df$D3.new1[i] - 1
            }
          } 
        }
      } else {
        acc <- df$J.start[i] - df$V.end[i] - 1
        if (acc < 0) {
          acc <- 0
        }
      }
      
      df$Total.insertions[i] <- acc
    }
    

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


vis.div <- function (.data) {
  if (!has.class(.data, "list")) {
    .data <- list(Data = .data)
  }
  
  df2 <- lapply(.data, function (x) {
    y <- as.data.frame(group_by(x, D.new) %>% summarise(Count = n()))
    y$Freq <- y$Count / sum(y$Count)
    y
  })
  df <- do.call(rbind, lapply(1:length(.data), 
                          function (i) { 
                            cbind(Subject = names(.data)[i], 
                                  df2[[i]]) })
                          )
  ggplot() + 
    geom_bar(aes(x = D.new, y = Freq, fill = Subject), colour = "black",
             data = df, stat = "identity", position = "dodge") +
    theme_linedraw() + 
    theme(axis.text.x  = element_text(angle=90)) + 
    ggtitle("Length D1,D2,D3 - 3; D2D3 - 5")
}