stripping_in_anger<-function(anger_frame){

  rc_primer_1<-revcomp(primer_1)
  rc_primer_2<-revcomp(primer_2)
  anger_frame[,c("polydt_loc","revdt_loc","polyt","polya")]<-NA
  print("Finding adapter sequence locations...")
  anger_frame$polyt<-str_locate_all(string = anger_frame$fastq_files, pattern = "T{12,}+") %>%
    sapply(., function(x){tail(as.vector(x), n = 1)}) %>%
    replace(., !sapply(., length),NA) %>%
    unlist(.)
  anger_frame$polya<-str_locate_all(string = anger_frame$fastq_files, pattern = "A{12,}+") %>%
    sapply(., function(x){head(as.vector(x), n = 1)}) %>%
    replace(., !sapply(., length),NA) %>%
    unlist(.)
  anger_frame$polydt_loc<-afind(x = anger_frame$fastq_files, pattern = primer_2, method = "lv")$location[,1]
  anger_frame$revdt_loc<-afind(x = anger_frame$fastq_files, pattern = rc_primer_2, method = "lv")$location[,1]
  forwards<-which(anger_frame$r1_dist<anger_frame$r1rev_dist & anger_frame$polydt_dist<anger_frame$revdt_dist & !is.na(anger_frame$polyt))
  reverse<-which(anger_frame$r1_dist>anger_frame$r1rev_dist & anger_frame$polydt_dist>anger_frame$revdt_dist & !is.na(anger_frame$polya))
  print("Sequences located! Removing adapters...")
  anger_frame$cleaned_seq<-NA
  anger_frame$cleaned_seq[forwards]<-mapply(x = anger_frame$fastq_files[forwards],
                                            y = anger_frame$polyt[forwards],
                                            z = anger_frame$polydt_loc[forwards],function(x,y,z){gsub(
                                              pattern = substr(x, start = 1, stop = y),
                                              replacement = "", x = x, fixed = TRUE) %>%
                                                gsub(pattern = substr(x, start = z, stop = str_length(x)),
                                                     replacement = "", x = ., fixed = TRUE)
                                            })
  anger_frame$cleaned_seq[reverse]<-mapply(x = anger_frame$fastq_files[reverse],
                                           y = anger_frame$polya[reverse],
                                           z = anger_frame$revdt_loc[reverse],function(x,y,z){gsub(
                                             pattern = substr(x, start = 1, stop = z+str_length(rc_primer_2)-1),
                                             replacement = "", x = x, fixed = TRUE) %>%
                                               gsub(pattern = substr(x, start = y, stop = str_length(x)),
                                                    replacement = "", x = ., fixed = TRUE)
                                           })
  anger_frame$error_code<-paste(anger_frame$r1_dist, anger_frame$r1rev_dist, anger_frame$polydt_dist, anger_frame$revdt_dist, sep = "_")
  return(anger_frame)
}
