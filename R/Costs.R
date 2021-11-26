source('SimulationsFunctions.R')
source('Parameters.R')



########## "simple" time-dependent cost functions



cgamma=2.5*c(2.5,2.5)/15 # cost of failing to treat
cbeta=2.5*c(1,1)/15 # cost of wrong treatment
M=N+H/8*cbeta[1] # death cost


tablecostsimple<-matrix(0,ncol=4,nrow=4) #cost matrix, to be multiplied by $r$
tablecostsimple[1,]<-c(0,cbeta,0)
tablecostsimple[2,]<-c(cgamma[1],0,cbeta[2],0)
tablecostsimple[3,]<-c(cgamma[2],cbeta[1],0,0)
tablecostsimple[4,]<-c(0,0,0,M)


costsimple<-function(m,s) # simple cost function in non-observed state space 
{

	t=s[1]
	if (t=="non") tt=1 else if (t=="a") tt=2 else if (t=="b") tt=3 else tt=4
	return (as.numeric(s[2])*tablecostsimple[m,tt])
}

costsimpletheta<-function(GammaIndex,s) # simple cost function in non-observed state space for element GammaIndex of Gamma with decision s = (treatment, next visit)
{
	if (GammaIndex<ngamma)
	{
		Theta=Gamma[,GammaIndex]
		Cint=sapply(Omega[,1]+1,costsimple,s=s) # note that mode starts at 0, so need to indent
		Res=Theta%*%Cint
		return(Res)
	} else return(M)
}

Finalcosttheta<-function(GammaIndex) # cost at final visit
{
	if (GammaIndex==ngamma) return (M) # death cost
	Theta=Gamma[1:(nOmega-1),GammaIndex]
	V=Omega[,2]-x0
	Res=0*Theta%*%V
	return(Res)
}


########## "marker-dependent" cost function

cbetai=1.5*c(1,1)/15 # cost of wrong treatment
cdivide=6 # how much to divide to compensate for marker value compared to visit cost
MI=N/2+H/8*cbetai[1] # death cost

tablevalue<-matrix(0,ncol=4,nrow=4) #cost matrix, to be multiplied by $r$
tablevalue[1,]<-c(0,cbetai,0)
tablevalue[2,]<-c(0,0,0,0)
tablevalue[3,]<-c(0,0,0,0)
tablevalue[4,]<-c(0,0,0,MI)




ftemp<-function(prob,z0r)
{
	return(prob%*%z0r)
}

FinalcostthetaI<-function(GammaIndex) # cost at final visit
{
	if (GammaIndex==ngamma) return (MI) 
	Theta=Gamma[1:(nOmega-1),GammaIndex]
	V=Omega[,2]-x0
	Res=Theta%*%V
	return(0)
}

costvaluetheta<-function(GammaIndex,s) # marker-dependent cost function in non-observed state space for element GammaIndex of Gamma with decision s = (treatment, next visit)
{

	if (GammaIndex==ngamma) return (MI)
	
	t=s[1]
	if (t=="non") tt=1 else if (t=="a") tt=2 else if (t=="b") tt=3 else tt=4
	indexd=(which(Time==as.numeric(s[2]))-1)*3+tt  # index that corresponds to decision s in transition matrices 
	
	

	zn0r=(Omega[,2]-x0)*as.numeric(s[2])/cdivide # value of |z_n+1 - z_0|*r
	cr=as.numeric(s[2])*(tablevalue[Omega[,1]+1,tt]) # penalty on decision
		
	Theta=Gamma[,GammaIndex] # actual filter corresponding to element GammaIndex
	T1=Theta%*%cr # penalty integrated over filter (as mode is not observed) 
	
	T2=apply(P[1:(nOmega-1),1:(nOmega-1),indexd],1,ftemp,z0r=zn0r) +P[1:(nOmega-1),nOmega,indexd]*MI # for a given x, integrating |z_n+1 - z_0|*r over P(z_n+1|x)
	T2int=Theta%*%T2 # and then integrating over filter
	

	
	return(as.numeric(T1+T2int))

}




