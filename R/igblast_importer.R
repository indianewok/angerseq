igblast_importer<-function(file_path){
  df<-fread(file_path, stringsAsFactor = TRUE, na.strings = "")
  df<-utils::type.convert(df, as.is = TRUE)
  return(df)
}
