#Run the code below step by step.
#From gab krouk (gkrouk@gmail.com)
####Import the Omaley dataset. This is long! 3.68 million lines!
read.table("DAP.seq.database.Dec.2022.txt", header=T)->DAP.Seq.Network
dim(DAP.Seq.Network)



####Import TF annotations showing up in the final table

read.table("TF.descriptions.txt", header=FALSE, sep='',row.names=1)->TF.descriptions
dim(TF.descriptions)
TF.descriptions[1:3,]
row.names(TF.descriptions)



#### makes a simple funtion that reads the gene list of interrest that you will put in gene.list.txt.
i<-function()
{
	
read.table('gene.list.txt', sep='\t', header=F)->a
as.character(a[1:length(t(a)),])->b
	
	return(b)
	
	
	
	
	}

i()


#########################################################################################################
# Install necessay packages 
#########################################################################################################

install.packages('parallel')
library('parallel')


#########################################################################################################
# The function below makes some sort of memoization to accelerate the process/
# it needs to be ran onces to generates "DAP.Seq.Network.Fast.list"
#########################################################################################################

Make.List.regine.multicore.gb<-function(gene.list=union(a[,1],a[,2]))
{
 
  
  DAP.Seq.Network->a
  union(a[,1],a[,2])->all.genes
  union(a[,1],a[,1])->all.TFs
  unique(a[,1])->all.TFs
  
  Find.core.gb<-function(gene){
    
    return(a[which(a[,2]==gene), 1])
    

  }
  
  mclapply(as.character(gene.list), Find.core.gb) -> my.list
  

  
 names(my.list)<-toupper(gene.list)

  
  return(my.list)
  
  
  
}

###########################################################################################################
### Run this it takes a Loooong while (memoization)

Make.List.regine.multicore.gb()->DAP.Seq.Network.Fast.list

###########################################################################################################
# this prepares a dataset counting the number of connections per TF.
# also takes a while but I also have to run it once only.
###########################################################################################################

TF.number.of.connection.counter.gb<-function()
  {
  
  DAP.Seq.Network->a
  union(a[,1],a[,2])->all.genes
  union(a[,1],a[,1])->all.TFs
  matrix(nrow=length(all.TFs),ncol= 2)->out
  
  Counter.multicor.gb<-function(gene)
    {
  
  
    dim(a[grep(gene,a[,1],ignore.case=T),])[1]
  
  }
  
  mclapply(as.character(all.TFs),  Counter.multicor.gb)->list1
  for(i in 1:length(list1))
  {
    
    list1[[i]]->out[i,2]
    all.TFs[i]->out[i,1]
    
    
  }
  
  
  
   return(out)
}

TF.number.of.connection.counter.gb()->DAP.Seq.Network.TF.connections.count
row.names(DAP.Seq.Network.TF.connections.count)<-DAP.Seq.Network.TF.connections.count[,1]
dim(DAP.Seq.Network.TF.connections.count)
DAP.Seq.Network.TF.connections.count

################
# Now the actual Regine script in 2 pieces.
# one (Regine.plot) generates a plot to see the potential deviation with randomness (as compared to random genes)
# the other one find the TF that are potentially regulators 
################



Regine.plot.gb<-function(gene.list, rep.per.range.value=100, kernel.range=2:9, name.your.list="New.list")
{

  #install.packages("ggplot2")
  require(ggplot2)
 
  as.character(as.vector(gene.list))->gene.list
  #print('yes 1')
  DAP.Seq.Network->a
  as.matrix(a)->a
  TFs<-union(a[,1],a[,1])
  All<-union(a[,1],a[,2])
  
  
  
  Mat<-matrix(ncol=3)
  
  for(i in kernel.range)
  {
    
   for(j in 1:rep.per.range.value)
   {
     
     #print(i)
     #print(j)
     
    sample(gene.list,i)->g
    sample(All,i)->h
    
    
    Find.me.direct.regulators.DAP.seq.based.multicore.V2.gb(g)->reg
    Find.me.direct.regulators.DAP.seq.based.multicore.V2.gb(h)->noise
    
    
    length(reg)->lreg
    length(noise)->lnoise
    
    
    #print(lreg)
    #print(lnoise)
    
    rbind(Mat,c(lreg,"l",paste(0,i,sep='')))->Mat
    rbind(Mat,c(lnoise,"b",paste(0,i,sep='')))->Mat
    
       
    
    
   }
    }
  
 
  Mat[2:length(Mat[,1]),]->Mat
  as.data.frame(Mat)->Mat
  Mat[,1]<-as.numeric(as.vector(Mat[,1]))
  Mat[,2]<-as.factor(Mat[,2])
  Mat[,3]<-as.numeric(as.character((Mat[,3])))
  colnames(Mat)<-c('numb_of_Reg','lb','ks')
  
  
  
  paste('Regine.',name.your.list,'.pdf',sep='')->File1
  
  
  pdf(File1, width=8, height=6)
  #par(mfrow=c(1,2))
  
  boxplot(Mat[,1]~Mat[,2]*Mat[,3],col=rep(c("purple", "lightblue")),xlab='Kernel size',ylab='# of regulators', main='b = background; l = your list')
  Mat->Mat2
  log(Mat2[,1])->Mat2[,1]
  print(qplot(numb_of_Reg, ks, data = Mat2, color = factor(lb), na.rm = T , geom = c("point", "smooth"))+ coord_flip())
  
 
  dev.off()
  
  return(Mat) 
  
  #print(summary(boxplot(as.numeric(as.vector(Mat[,1]))~as.factor(Mat[,2])*as.factor(Mat[,3]))))
 
  
  
}



########################




Find.me.direct.regulators.KERNEL.STAT.gb <- function(gene.list,kernel.size=5,rep=5, name.your.list="new_list")
{
  
  paste('DAP.seq.live.report.',name.your.list,'.size',length(gene.list),'.kernel',kernel.size,'.rep',rep,'.txt',sep='')->File1
  paste('DAP.seq.results.',name.your.list,'.size',length(gene.list),'.kernel',kernel.size,'.rep',rep,'.txt',sep='')->File2
  
  
  ####
  #install.packages("devtools", dep = T)
  library(devtools)
  #install_github("rpremraj/mailR", force=T)
  
  #library(mailR)
  #####
  
  
  as.character(as.vector(gene.list))->gene.list
  #print('yes 1')
  DAP.Seq.Network->a
  as.matrix(a)->a
  TFs<-union(a[,1],a[,1])
  All<-union(a[,1],a[,2])
  
  
  ##### Background
  #intersect(AGI.non.ambig,TFs)->TFs
  #intersect(AGI.non.ambig,All)->All
  #####
  
  #print(length(TFs))
  #print(length(All))
  
  b<-matrix(nrow=length(TFs),ncol=rep)
  n<-matrix(nrow=length(TFs),ncol=rep)
  
  row.names(b)<-TFs
  row.names(n)<-TFs
  #print(b)
  
  cat(paste('Gene lists size =',length(gene.list)),file=File1,fill=T, append=F)
  cat(paste('Kernel =',kernel.size),file=File1,fill=T, append=TRUE)
  cat(paste('potential combinations =',as.numeric(choose(length(gene.list),kernel.size))->comb),file=File1,fill=T, append=TRUE)
  #print(comb)
  cat(paste('Repetitions of sampling =',rep,"which corresponds to",rep/comb*100,'% of the possibilities'),file=File1,fill=T, append=TRUE)
  cat('This script extracts common regulators from the Omaley et al [Ecker lab] dataset.',file=File1,fill=T, append=TRUE)
  cat('   ',file=File1,fill=T, append=TRUE)
  
  
  
  for(i in 1:rep)
  {
    
    cat(paste('cycle #',i,"on",rep),file=File1,fill=T, append=TRUE)
    sample(gene.list,kernel.size)->g
    #setdiff(All,gene.list)->rest.of.genome
    #sample(rest.of.genome,kernel.size)->h #comment stat de Clem
    sample(All,kernel.size)->h
    #Find.me.direct.regulators.DAP.seq.based.gb(g)->reg
    #Find.me.direct.regulators.DAP.seq.based.gb(h)->noise
    
    Find.me.direct.regulators.DAP.seq.based.multicore.V2.gb(g)->reg
    Find.me.direct.regulators.DAP.seq.based.multicore.V2.gb(h)->noise
    
    
    intersect(reg,noise)->inter
    setdiff(reg,inter)->reg.clean
    
    
    cat(paste('# Regulators =',length(reg),"; # of background regulators (this is noise) =", length(noise), '; # of regulators DeNoised',length(reg.clean)),file=File1,fill=T, append=TRUE)
    cat('List of TARGET:',file=File1,fill=T, append=TRUE)
    cat(g,file=File1,fill=T, append=TRUE)
    cat('are regulated by:',file=File1,fill=T, append=TRUE)
    cat(reg,file=File1,fill=T, append=TRUE)
    cat('   ',file=File1,fill=T, append=TRUE)
    
    if(length(reg)>0)
    {
      b[reg,i]<-1
    }
    else{}
    
    
    if(length(noise)>0)
    {
      n[noise,i]<-1
    }
    else{}
    
  
    
  }
  
  rowSums(b,na.rm = T)->d
  sort(d[d>0],decreasing=T)->regulators.out
  
  
  
  rowSums(n,na.rm = T)->o
  sort(o[o>0],decreasing=T)->noise.out	
  
  union(names(regulators.out),names(noise.out))->out.TF.rslt
  
  
  matrix(nrow=length(out.TF.rslt),ncol=4)->r
  row.names(r)<-out.TF.rslt
  
  r[names(regulators.out),1]<-regulators.out
  r[names(noise.out),2]<-noise.out
  r[is.na(r)]<-0
  r[,1]-r[,2]->r[,3]
  
  DAP.Seq.Network.TF.connections.count[row.names(r),2]->r[,4]
  
  
  r[order(as.numeric(as.matrix(r)[,3]),decreasing = T),]->r
  cbind(r,TF.descriptions[row.names(r),])->r
  
  colnames(r)<-c('Real','Noise','Real-Noise',"TF.out.connections","AGI","Description")

  
  
  
  write.table(r,file=File2,quote=F)

  cat("I'm done. You will find your results in DAP.seq.Results.txt",file=File1,fill=T, append=TRUE)
  

  
  
  
  
  
  
  
  
  
  
  return(r)
  
}


Find.me.direct.regulators.DAP.seq.based.multicore.V2.gb <- function(gene.list){
  
  

  return(Reduce(intersect,DAP.Seq.Network.Fast.list[toupper(gene.list)]))

  
  
}


###### you are now ready to go to the Regine.RUN file.
###### you are now ready to go to the Regine.RUN file.
###### you are now ready to go to the Regine.RUN file.
###### you are now ready to go to the Regine.RUN file.




