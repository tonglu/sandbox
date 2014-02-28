library(seqinr)
seq<-read.fasta(file = "sequence.fasta")
seq<-seq[[1]]
transition.matrix <- read.table("C:/Users/Lohas/Desktop/hw5-cpg-transition-matrix (1).txt", header=T, quote="\"")

#encode sequence
seq[seq=="a"]<-1
seq[seq=="c"]<-2
seq[seq=="g"]<-3
seq[seq=="t"]<-4
seq<-as.integer(seq)
n=length(seq)

#initialize
#initial distribution of p(z1)
p.z1<-rep(c(sum(seq==1)/(2*length(seq)),sum(seq==2)/(2*length(seq)),sum(seq==3)/(2*length(seq)),sum(seq==4)/(2*length(seq))),2)
x1<-seq[1]
p.x1.z1<-rep(c(0,0,1,0),2)
w1<-rep(0,2)
w1<-p.z1[which(p.x1.z1==1)] #w1
#record z1's position 
posi1<-which(p.x1.z1==1)
position<-matrix(0,2,n)
position[,1]<-posi1

#record the Zi that Wi(Zi)!=0
nonzero<-matrix(0,2,n) 
nonzero[,1]<-posi1    
for (k in 2:n){
  nonzero[1,k]<-seq[k]
  nonzero[2,k]<-seq[k]+4
}

#recursion
n=length(seq)
log.w<-matrix(0,2,n)
log.w[,1]<-log(w1)
for (k in 2:n){                           
    if (log(transition.matrix[nonzero[1,k-1],nonzero[1,k]])+log.w[1,k-1]>log(transition.matrix[nonzero[2,k-1],nonzero[1,k]])+log.w[2,k-1]){
      position[1,k]<-nonzero[1,k-1]
      log.w[1,k]<-log(transition.matrix[nonzero[1,k-1],nonzero[1,k]])+log.w[1,k-1]
    }
    else {
         position[1,k]<-nonzero[2,k-1]
         log.w[1,k]<-log(transition.matrix[nonzero[2,k-1],nonzero[1,k]])+log.w[2,k-1]
    }
    
    if (log(transition.matrix[nonzero[1,k-1],nonzero[2,k]])+log.w[1,k-1]>log(transition.matrix[nonzero[2,k-1],nonzero[2,k]])+log.w[2,k-1]){
      position[2,k]<-nonzero[1,k-1]
      log.w[2,k]<-log(transition.matrix[nonzero[1,k-1],nonzero[2,k]])+log.w[1,k-1]
    }
    else {
      position[2,k]<-nonzero[2,k-1]
      log.w[2,k]<-log(transition.matrix[nonzero[2,k-1],nonzero[2,k]])+log.w[2,k-1]
    }  
}

#back-tracking
log.w[,n] 
#log.w[1,n]<log.w[2,n], so back-tracking indicates we should choose second row
q<-matrix(0,1,n)
q<-position[2,]
q[q>=5]<-0
q[q>0]<-1
plot(q,type='l',ylab='Threshold',xlab="Base Number")



