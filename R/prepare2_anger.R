prepare2_anger<-function(file_path, path_to_seqkit, path_to_output){
  fastq_files<-list.files(path = paste0(file_path), full.names = TRUE, pattern = ".gz")
  if(length(fastq_files) == 1){
    seqkit<-path_to_seqkit
    output_directory<-paste0(path_to_output, basename(dirname(fastq_files)), "_split")
    if(file.exists(output_directory)==TRUE){
      output_directory<-paste0(output_directory,"mem")
    }
    dir.create(path = output_directory)
    {
      seqkit_stats_args<-c("stats",fastq_files,"-T","-a")
      stats<-processx::run(seqkit, seqkit_stats_args, echo_cmd =FALSE, echo = FALSE, spinner = TRUE)
      stats<-strsplit(x = stats$stdout, split = "\n") %>%
        lapply(., function(x){strsplit(x, split = "\t")}) %>%
        unlist(.,recursive = FALSE) %>%
        data.frame(.[[2]], row.names = .[[1]]) %>%
        .[,c(1,3)]
      colnames(stats)<-NULL
      stats[,1]<-NULL
    }#global statistics
    {
      seqkit_split2_args<-c("split2", fastq_files, "--by-part","1000","--by-part-prefix",
                            paste0(basename(dirname(fastq_files)), "_split_"),
                            "-O",output_directory,
                            "-e",".gz")
      split<-processx::run(seqkit, seqkit_split2_args, echo_cmd =FALSE, echo = FALSE, spinner = TRUE)
      file_path<-output_directory
      fastq_files<-list.files(path = paste0(file_path), full.names = TRUE, pattern = ".gz")
      output<-list(fastq_files, stats)
      return(output)
    }#splitting the files
  }
}
