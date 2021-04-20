awk -F ',' '{ print $1 }' heart_lung_skin_samples.txt > samplenames.txt


awk '
FNR==NR{
  for(i=1;i<=NF;i++){
    arr[$i]
  }
  next
}
FNR==3{
  for(i=1;i<=NF;i++){
    if($i in arr){
      valArr[i]
      header=(header?header OFS:"")$i
    }
  }
  print header
  next
}
{
  val=""
  for(i=1;i<=NF;i++){
    if(i in valArr){
       val=(val?val OFS:"")$i
    }
  }
  print val
}
' samplenames.txt FS="\t" OFS="\t" GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct

#for i in "${arr[@]}"
#do
#    arrind=($(sed -n 3p GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct | awk -F '\t' -v a=$i '{ for(j=1;j<=NF;j++) { if($j==a) { print j }}}'))
#    for j in "${arrind[@]}"
#    do
#	awk -F '\t' -v a=$j '{ for(j=1;j<=NF;j++) { if(j=a) {print $j}}}' GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct
#    done
#done
