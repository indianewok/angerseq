erred_in_anger<-function(anger_frame){
  concatenated_indx<-sapply(concatenates, function(x){grep(x = anger_frame$error_code, pattern = x, value = FALSE)}, simplify = TRUE) %>%
    unlist(.) %>%
    unique(.)
  concatenated<-anger_frame[concatenated_indx,]
  splitseqs<-pbmclapply(concatenated$fastq_files, function(x){concatenate_rescue(x)}, mc.cores = 8)
  splitseqs<-unlist(splitseqs, recursive = FALSE)
  concatenated<-concatenated[-which(sapply(splitseqs, function(x){length(x)}) > 2),]
  splitseqs<-splitseqs[-which(sapply(splitseqs, function(x){length(x)}) > 2)]

  for(i in 1:length(splitseqs)){
    splitseqs[[i]]<-data.frame(splitseqs[[i]], concatenated$id[i])
    colnames(splitseqs[[i]])<-c("fastq_files","id")
    splitseqs[[i]]$qc<-str_locate_all(string = concatenated$fastq_files[i], pattern = splitseqs[[i]]$fastq_files) %>%
      sapply(., function(x){str_sub(string = concatenated$qc[i], x)})
  }
  splitseqs<-rbindlist(splitseqs) %>% .[,c(2,3,1)]
  splitseqs<-left_join(splitseqs, concatenated[,c("r1_dist","r1rev_dist","polydt_dist","revdt_dist","id","error_code")], by= "id")
  rm(concatenated)

  splitseqs<-anger_distance(anger_frame = splitseqs, primer_1 = primer_1, primer_2) %>%
    anger_pbarc(whitelist_path = whitelist_path, anger_frame = ., adapter_distance = 5, primer_1 = primer_1) %>%
    stripping_in_anger(.)
  splitseqs$putative_bcs[which(sapply(splitseqs$cleaned_seq, function(x){str_length(x)})==0)]<-NA
  splitseqs$cleaned_seq[which(sapply(splitseqs$cleaned_seq, function(x){str_length(x)})==0)]<-NA

  anger_frame<-anger_frame[-match(unique(splitseqs$id), table = anger_frame$id)]
  splitseqs$id<-paste0(splitseqs$id, " splitid=",seq(nrow(splitseqs)))
  anger_frame<-rbind(anger_frame, splitseqs)

  anger_frame$putative_bc_log<-!is.na(anger_frame$putative_bcs)
  anger_frame$putative_bcs[which(sapply(anger_frame$putative_bcs, function(x){str_length(x)}) == 0)]<-NA
  anger_frame$qc_updated<-mapply(x = anger_frame$cleaned_seq, y = anger_frame$fastq_files, z = anger_frame$qc,
                                 function(x,y,z){
                                   str_locate(pattern = x, string = y) %>% str_sub(string = z, .)
                                 })
  anger_frame<-anger_frame[which(!is.na(anger_frame$putative_bcs) & !is.na(anger_frame$cleaned_seq)),]
  lazy_bcs<-anger_frame[which(anger_frame$cr_barcode == FALSE),"putative_bcs"]
  true_bcs<-anger_frame[which(anger_frame$cr_barcode == TRUE),"putative_bcs"]
  output<-list(anger_frame, lazy_bcs, true_bcs)
  names(output)<-c("anger_frame","lazy_bcs","true_bcs")
  return(output)
}
