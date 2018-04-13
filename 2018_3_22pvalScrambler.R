# Switching up 

netwas.input=read.table("NetWasPariMidInput.txt")

for(ii in 1:10){
xx=cbind(as.character(netwas.input[,1]),as.character(netwas.input[,1]),netwas.input[sample(1:nrow(netwas.input)),c(3,3,(3:ncol(netwas.input)))])
write.table(xx,file=paste0("NetwasScramble",ii,".txt"),col.names=FALSE,row.names=FALSE)
}



