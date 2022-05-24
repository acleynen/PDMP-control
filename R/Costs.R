

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


Finalcosttheta<-function(GammaIndex) # cost at final visit
{
	if (GammaIndex==ngamma) return (M) # death cost
	Theta=Gamma[,GammaIndex]
	return(FinalcostthetaTraj(Theta))
}

FinalcostthetaTraj<-function(Filt)
{
	if (sum(Filt)==0) return(M) # death cost
	V=Omega[,2]-x0
	Res=Filt%*%V
	return(Res)
}

costsimpletheta<-function(GammaIndex,s) # simple cost function in non-observed state space for element GammaIndex of Gamma with decision s = (treatment, next visit)
{
	if (GammaIndex==ngamma) return (M) # death cost
	Theta=Gamma[,GammaIndex]
	return(costsimplethetaTraj(Theta,s))
}


costsimplethetaTraj<-function(Filt,s) # la fonction de coût simple pour un filtre Filt
{
	if (sum(Filt)==0) return(M)
	Cint=sapply(Omega[,1]+1,costsimple,s=s) # Mode starts at 0, necessary index shift
	Res=Filt%*%Cint
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
	return(0)
}

FinalcostthetaITraj<-function(Filt)
{
	if (sum(Filt)==0) return(MI)  # death cost
	return(0)
}


costvaluetheta<-function(GammaIndex,s) # marker-dependent cost function in non-observed state space for element GammaIndex of Gamma with decision s = (treatment, next visit)
{

	if (GammaIndex==ngamma) return (MI)
	Theta=Gamma[,GammaIndex] # actual filter corresponding to element GammaIndex
	return(costvaluethetaTraj(Theta,s))

}


costvaluethetaTraj<-function(Filt,s)
{

	if (sum(Filt)==0) return(MI)
	
	t=s[1]
	if (t=="non") tt=1 else if (t=="a") tt=2 else if (t=="b") tt=3 else tt=4
	indexd=(which(Time==as.numeric(s[2]))-1)*3+tt # index that corresponds to decision s in transition matrices 
		

	zn0r=(Omega[,2]-x0)*as.numeric(s[2])/cdivide # value of |z_n+1 - z_0|*r
	cr=as.numeric(s[2])*(tablevalue[Omega[,1]+1,tt]) # penalty on decision
		 
	T1=Filt%*%cr # penalty integrated over filter (as mode is not observed) 
	
	T2=apply(P[1:(nOmega-1),1:(nOmega-1),indexd],1,ftemp,z0r=zn0r) +P[1:(nOmega-1),nOmega,indexd]*MI # for a given x, integrating |z_n+1 - z_0|*r over P(z_n+1|x)
	T2int=Filt%*%T2 # and then integrating over filter	
	return(as.numeric(T1+T2int))

}

########## "Threshold and treatment length" cost function


cbetaiCT=4.5*c(1,1)/15 # cost of wrong treatment
cdivideCT=6 # how much to divide to compensate for marker value compared to visit cost
ThreshCT=15 # will start paying only if marker value is above this threshold
#MI=N/2+H/8*cbetai[1] # death cost
MICT=200

tablevalueCT<-matrix(0,ncol=4,nrow=4) #cost matrix, to be multiplied by $r$
tablevalueCT[1,]<-c(0,cbetaiCT,0)
tablevalueCT[2,]<-c(0,0,0,0)
tablevalueCT[3,]<-c(0,0,0,0)
tablevalueCT[4,]<-c(0,0,0,MICT)







FinalcostthetaICT<-function(GammaIndex) # cost at final visit
{
	if (GammaIndex==ngamma) return (MICT) 
	Theta=Gamma[1:(nOmega-1),GammaIndex]
	return(FinalcostthetaITrajCT(Theta))
}

FinalcostthetaITrajCT<-function(Filt)
{
	if (sum(Filt)==0) return(MICT)  
	V=Omega[,2]-x0
	Res=Filt%*%V
	return(as.numeric(0))
}


costvaluethetaCT<-function(GammaIndex,s) # marker-dependent cost function in non-observed state space for element GammaIndex of Gamma with decision s = (treatment, next visit)
{

	if (GammaIndex==ngamma) return (MICT)
	Theta=Gamma[,GammaIndex]
	return(costvaluethetaTrajCT(Theta,s))

}


costvaluethetaTrajCT<-function(Filt,s)
{

	if (sum(Filt)==0) return(MICT)
	
	t=s[1]
	if (t=="non") tt=1 else if (t=="a") tt=2 else if (t=="b") tt=3 else tt=4
	indexd=(which(Time==as.numeric(s[2]))-1)*3+tt  # index that corresponds to decision s in transition matrices 
	
	
	if (tt==1)
	{
		zn0r=(Omega[,2]-x0)*(Omega[,2]>ThreshCT)*as.numeric(s[2])/cdivideCT # value of |z_n+1 - z_0|*r if z_n+1>Threshold
		pro=	P[1:(nOmega-1),1:(nOmega-1),indexd]
	}	
	if (tt==2)
	{
		znTr=(Omega[,2]-x0)*as.numeric(s[2])/cdivideCT ## value of |z_n+1 - z_0|*r
		zn0r=znTr[(nm0+nm1+1):(nm0+nm1+nm2)]
		pro=	P[1:(nOmega-1),(nm0+nm1+1):(nm0+nm1+nm2),indexd]
	}	
	if (tt==3)
	{
		znTr=(Omega[,2]-x0)*as.numeric(s[2])/cdivideCT # value of |z_n+1 - z_0|*r
		zn0r=znTr[(nm0+1):(nm0+nm1)]
		pro=	P[1:(nOmega-1),(nm0+1):(nm0+nm1),indexd]
	}
	cr=as.numeric(s[2])*(tablevalueCT[Omega[,1]+1,tt]) # penalty on decision
		 
	T1=Filt%*%cr # penalty integrated over filter (as x is not observed) 
	
	T2=apply(pro,1,ftemp,z0r=zn0r) +P[1:(nOmega-1),nOmega,indexd]*MICT # for a given x, integrating |z_n+1 - z_0|*r over P(z_n+1|x)
	T2int=Filt%*%T2
	
	
	return(as.numeric(T1+T2int))

}

