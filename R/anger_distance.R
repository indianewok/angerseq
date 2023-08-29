anger_distance<-function(anger_frame, primer_1, primer_2){
  rc_primer_1<-revcomp(primer_1)
  rc_primer_2<-revcomp(primer_2)
  if(any(class(anger_frame) == "list")){
    progress<-progress::progress_bar$new(total = length(anger_frame), format = "  progress [:bar] :percent eta: :eta",)
    for(i in 1:length(anger_frame)){
      anger_frame[[i]][,c("r1_dist","r1rev_dist","polydt_dist","revdt_dist")]<-NA
      anger_frame[[i]]$r1_dist<-afind(x = anger_frame[[i]]$fastq_files, pattern = primer_1, method = "lv")$distance[,1]
      anger_frame[[i]]$r1rev_dist<-afind(x = anger_frame[[i]]$fastq_files, pattern = rc_primer_1, method = "lv")$distance[,1]
      anger_frame[[i]]$polydt_dist<-afind(x = anger_frame[[i]]$fastq_files, pattern = primer_2, method = "lv")$distance[,1]
      anger_frame[[i]]$revdt_dist<-afind(x = anger_frame[[i]]$fastq_files, pattern = rc_primer_2, method = "lv")$distance[,1]
      progress$tick()
    }
  }
  else{
    print("Finding angry distances on a single frame...")
    anger_frame[,c("r1_dist","r1rev_dist","polydt_dist","revdt_dist")]<-NA
    anger_frame$r1_dist<-afind(x = anger_frame$fastq_files, pattern = primer_1, method = "lv")$distance[,1]
    anger_frame$r1rev_dist<-afind(x = anger_frame$fastq_files, pattern = rc_primer_1, method = "lv")$distance[,1]
    anger_frame$polydt_dist<-afind(x = anger_frame$fastq_files, pattern = primer_2, method = "lv")$distance[,1]
    anger_frame$revdt_dist<-afind(x = anger_frame$fastq_files, pattern = rc_primer_2, method = "lv")$distance[,1]
  }
  return(anger_frame)
}
