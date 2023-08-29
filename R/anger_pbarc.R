anger_pbarc<-function(whitelist_path, anger_frame, adapter_distance, primer_1){
  progress<-progress::progress_bar$new(total = length(anger_frame),
                                       format = "  progress [:bar] :percent eta: :eta",)
  progress$tick(0)
  rc_primer_1<-revcomp(primer_1)
  whitelist<-data.table::fread(input = whitelist_path, header = FALSE)
  if(any(class(anger_frame)=="list")){
    print("Barcoding a list of frames...")
    for(i in 1:length(anger_frame)){
      anger_frame[[i]]$putative_bcs<-NA
      anger_frame[[i]]$cr_barcode<-NA
      fadapted<-which(anger_frame[[i]]$r1_dist <= adapter_distance)
      {
        anger_frame[[i]]$putative_bcs[fadapted]<-substr(
          x = anger_frame[[i]]$fastq_files[fadapted],
          start = str_length(
            afind(x = anger_frame[[i]]$fastq_files[fadapted],
                  pattern = primer_1, method = "lv")$match[,1])+
            afind(x = anger_frame[[i]]$fastq_files[fadapted],
                  pattern = primer_1, method = "lv")$location[,1],
          stop = str_length(
            afind(x = anger_frame[[i]]$fastq_files[fadapted],
                  pattern = primer_1, method = "lv")$match[,1])+
            afind(x = anger_frame[[i]]$fastq_files[fadapted],
                  pattern = primer_1, method = "lv")$location[,1]+15)
      }#forward-primer based barcode extraction
      radapted<-which(anger_frame[[i]]$r1rev_dist <= adapter_distance & is.na(anger_frame[[i]]$putative_bcs))
      {
        anger_frame[[i]]$putative_bcs[radapted]<-
          revcomp(substr(anger_frame[[i]]$fastq_files[radapted],
                         start = afind(x = anger_frame[[i]]$fastq_files[radapted],
                                       pattern = rc_primer_1, method = "lv")$location[,1]-16,
                         stop = afind(x = anger_frame[[i]]$fastq_files[radapted],
                                      pattern = rc_primer_1, method = "lv")$location[,1]-1))
      }#reverse-barcode based barode extraction
      anger_frame[[i]]$cr_barcode<-anger_frame[[i]]$putative_bcs %in% whitelist$V1
      progress$tick()
    }
    return(anger_frame)
  }
  else{
    print("Barcoding a single frame...")
    anger_frame$putative_bcs<-NA
    anger_frame$cr_barcode<-NA
    fadapted<-which(anger_frame$r1_dist <= adapter_distance)
    {
      anger_frame$putative_bcs[fadapted]<-substr(
        x = anger_frame$fastq_files[fadapted],
        start = str_length(
          afind(x = anger_frame$fastq_files[fadapted],
                pattern = primer_1, method = "lv")$match[,1])+
          afind(x = anger_frame$fastq_files[fadapted],
                pattern = primer_1, method = "lv")$location[,1],
        stop = str_length(
          afind(x = anger_frame$fastq_files[fadapted],
                pattern = primer_1, method = "lv")$match[,1])+
          afind(x = anger_frame$fastq_files[fadapted],
                pattern = primer_1, method = "lv")$location[,1]+15)
    }#forward-primer based barcode extraction
    radapted<-which(anger_frame$r1rev_dist <= adapter_distance & is.na(anger_frame$putative_bcs))
    {
      anger_frame$putative_bcs[radapted]<-
        revcomp(substr(anger_frame$fastq_files[radapted],
                       start = afind(x = anger_frame$fastq_files[radapted],
                                     pattern = rc_primer_1, method = "lv")$location[,1]-16,
                       stop = afind(x = anger_frame$fastq_files[radapted],
                                    pattern = rc_primer_1, method = "lv")$location[,1]-1))
    }#reverse-barcode based barode extraction
    anger_frame$cr_barcode<-anger_frame$putative_bcs %in% whitelist$V1
    return(anger_frame)
  }
}
