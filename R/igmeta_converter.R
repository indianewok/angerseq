igmeta_converter<-function(df){
  df$barcode<-gsub(pattern = "_.+", replacement = "", x = df$sequence_id)
  df<-df[which(df$productive == TRUE & !is.na(df$cdr3)),]
  df<-split.data.frame(df, f = df$barcode)
  reads<-sapply(df, function(x){nrow(x)})
  df<-df[which(reads >= 4)]
  progress<-txtProgressBar(min = 0, max=length(df), style = 3, char = "~")
  for(i in 1:length(df)){
    df[[i]]<-left_join(unique(df[[i]][,c("v_call","d_call","j_call","cdr3","cdr3_aa",
                                         "junction","junction_aa","locus",
                                         "barcode","sequence_alignment",
                                         "sequence_alignment_aa")]),
                       data.frame(table(df[[i]]$cdr3)), by = c("cdr3" = "Var1"))
    names(df)[[i]]<-unique(df[[i]]$barcode)
    df[[i]]<-lapply(unique(df[[i]]$locus), function(x){
      df[[i]][which(df[[i]]$locus == x),][which.max(df[[i]][which(df[[i]]$locus == x),"Freq"]),]
    })
    for(j in 1:length(df[[i]])){
      colnames(df[[i]][[j]])[-which(colnames(df[[i]][[j]]) == "barcode")]<-paste0(
        colnames(df[[i]][[j]][-which(colnames(df[[i]][[j]]) == "barcode")]), "_",df[[i]][[j]]$locus)
    }
    if(length(df[[i]]) > 1){
      df[[i]]<-df[[i]][which(sapply(df[[i]], function(x){nrow(x)}) > 0)]
      df[[i]]<-data.frame(unlist(df[[i]], recursive = FALSE))
    }
    setTxtProgressBar(progress, i)
  }
  rm(progress,i,j)
  df[which(sapply(df, function(x){class(x)}) == "list")]<-unlist(df[which(
    sapply(df, function(x){class(x)}) == "list")],recursive = FALSE)
  df<-rbindlist(df, fill = TRUE)
}
