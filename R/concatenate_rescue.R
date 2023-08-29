concatenate_rescue<-function(conc_seq){
  mtx<-str_locate(pattern = c(adapters), string = conc_seq)
  if(any(is.na(mtx))){
    mtx<-adapters[which(is.na(mtx), arr.ind = TRUE)[1,1]] %>%
      afind(x = conc_seq, pattern = ., method = "lv", nthread = 1) %>%
      unlist(.) %>% .[3] %>%
      str_locate(string = conc_seq, pattern = .) %>%
      replace_na(data = mtx, replace = .)
  }
  mtx<-mtx[order(mtx[,1], decreasing = FALSE),]
  if(mtx[3,1]-mtx[2,2] <= 2){
    seqs<-list(c(substr(conc_seq, start = mtx[1,1], stop = mtx[2,2]), substr(conc_seq, start = mtx[3,1], stop = mtx[4,2])))
    return(seqs)
  }
  else{
    seqs<-substr(x = conc_seq, start = mtx[2,2]+1, stop = mtx[3,1]-1) %>% str_split(string = conc_seq, pattern = ., simplify = FALSE)
    return(seqs)
  }
}
