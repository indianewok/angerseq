revcomp<-Vectorize(function(x){
  str_to_upper(c2s(rev(comp(s2c(x)))))
})
