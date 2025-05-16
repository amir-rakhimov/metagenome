df<-read.table("output/mag_assembly/2D10_cov.tsv",
               comment.char = "",
               header = T,
               sep = "\t")
df.tpm<-df%>%
  as_tibble()%>%
  mutate(RPK=(Plus_reads+Minus_reads)/(Length/1000),
         RPK_sum=sum(RPK),
         TPM=(RPK/RPK_sum)*10^6)

k141_112109 
k141_39531

