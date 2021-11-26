MachineMin=.Machine$double.xmin

################ Trajectories

TrajExp<-function(x,t,tprime,m,j) # x= marker value at time t, tprime=time for which next marker value is wanted, m=current mode (at time t) and j=treatment 
{
	if (m==0)
		return (x0)
	if (m==1)
		if (j=="a")
			return (x*exp(-vprime[1]*(tprime-t))) else if (j=="b") return (x*exp(v1*(tprime-t))) else return (x*exp(v1not*(tprime-t)))
	if (m==2)
		if (j=="b")
			return (x*exp(-vprime[2]*(tprime-t))) else if (j=="a") return (x*exp(v2*(tprime-t)))	else return (x*exp(v2not*(tprime-t)))	
	if (m==100)
		return (D)
}

################ Jump time

Horloge<-function(x,t,m,j) # returns next jump time, and next mode. 
{

	if (m==0)
	{
		g1=rbinom(1,1,Pi[1])
		g2=rbinom(1,1,Pi[2])
		if (g1==1)
			t1=(t^(alpha[1]+1)-log(1-runif(1))*(alpha[1]+1)/(lambda[1]^alpha[1]))^(1/(alpha[1]+1))-t else t1=rnorm(1,2000,400)
		if (g2==1)
			t2=(t^(alpha[2]+1)-log(1-runif(1))*(alpha[2]+1)/(lambda[2]^alpha[2]))^(1/(alpha[2]+1))-t else t2=rnorm(1,2000,400)
		if (j=="non")
		{
			T=min(t1,t2)
			mod=which.min(c(t1,t2))
			return(c(T,mod))
		} 
		if (j=="a")
			return(c(t2,2))
		if (j=="b")
			return(c(t1,1))	
	}	
	
	if (m==1)
	{
		if (j=="a")
		{
			tfloor=1/vprime[m]*log(x/x0)
			t2=-1/(beta[2]*vprime[m]) * log( 1 + beta[2]*vprime[m]/((b[2]*x)^beta[2])*log(1-runif(1)))
			T=min(tfloor,t2)
			if (T==tfloor)
				return(c(tfloor,0)) else return(c(t2,2))			
		}	else if (j=="b")
		{
			tdeath=1/v1*log(D/x)
			return(c(tdeath,100))
		}	else
		{
			tdeath=1/v1not*log(D/x)
			return(c(tdeath,100))
		}			
	}
 
	if (m==2)
	{
		if (j=="b")
		{
			tfloor=1/vprime[m]*log(x/x0)
			t1=-1/(beta[1]*vprime[m]) * log( 1 + beta[1]*vprime[m]/((b[1]*x)^beta[1])*log(1-runif(1)))
			T=min(tfloor,t1)
			if (T==tfloor)
				return(c(tfloor,0)) else return(c(t1,1))			
		}	else if (j=="a")
		{
			tdeath=1/v2*log(D/x)
			return(c(tdeath,100))
		}	else
		{
			tdeath=1/v2not*log(D/x)
			return(c(tdeath,100))
		}		
	} else return(c(H,100))	
}



############### Simulation of a controlled trajectory over a time step 

Step<-function(num,m,x,t,w,r,j) # first argument unused: for paralel simulations. m=current mode, x=current marker value, t=current time since last jump, w=time since beginning, r=next visit, j=treatment 
{
	j1=Horloge(x,t,m,j)
	tjump=j1[1]
	mjump=j1[2]
	if (r>tjump)
	{
		posint=TrajExp(x,t,t+tjump,m,j)
		if (((m==1) && ((mjump==0) || (mjump==2)))  || ((m==2) && ((mjump==0) || (mjump==1)))) # checks for double jump during time step
		{
			hint=Horloge(posint,0,mjump,j)
			if ((r-tjump)>hint[1])
			{
				posint=TrajExp(posint,0,hint[1],mjump,j)
				tjump=tjump+hint[1]
				mjump=hint[2]
			}
		}
		xnew=TrajExp(posint,0,r-tjump,mjump,j)
		tnew=r-tjump
		mnew=mjump	
	} else
	{
		xnew=TrajExp(x,t,t+r,m,j)
		tnew=t+r
		mnew=m
	}
	return(c(mnew,xnew,tnew,w+r))
}


############### Simulation of a controlled trajectory over a time step without simulation of jump time if not necessary. In expectation same thing as previous function, but useful to compare different strategies on a "same trajectory". 

NewStep<-function(num,m1,x,t,w,r,j0,j1,tjump, mjump,TsH) # first argument unused: for paralel simulations. m1=current mode, x=current marker value, t=current time since last jump, w=time since beginning, r=next visit, j0=treatment at previous time, j1=current treatment, tjump=next jump time if no change, mjump=next mode if no change, TsH=time since last jump time was simulated 
{
	if (j0!=j1) # next jump time computed due to change in treatment
	{
		jumpa=Horloge(x,t,m1,j1)
		tjump=jumpa[1]
		mjump=jumpa[2]
		TsH=0
	}
	if (TsH+r>tjump)
	{
		posint=TrajExp(x,t,tjump,m1,j1)
		mnew=mjump
		hint=Horloge(posint,0,mjump,j1) # next jump time computed as jump just occured
		TsH=TsH+r-tjump
		tjump=hint[1]
		mjump=hint[2]		
		
		if (TsH>tjump) # checks for double jump
		{
			posint=TrajExp(posint,0,tjump,mjump,j1)
			hint=Horloge(posint,0,mjump,j1) # next jump time computed as jump just occured
			mnew=mjump
			tjump=hint[1]
			mjump=hint[2]		
		}
		xnew=TrajExp(posint,0,TsH,mnew,j1)
		tnew=TsH
	} else
	{
		xnew=TrajExp(x,t,t+r,m1,j1)
		tnew=t+r
		mnew=m1
		TsH=TsH+r
	}
	return(c(mnew,xnew,tnew,w+r,tjump,mjump,TsH))
}

############### Nearest neighbor in Omega

NNeighbor<-function(X)
{
	m=X[1] ;x=X[2]; t=X[3]
	if (m==100)
		return(c(100,0,0,nOmega))
		
	temp=Omega[Omega[,1]==m,]	
	if (m==0)
	{
		dist=(t-temp[,3])^2
		l=which.min(dist)
		return(c(temp[l,]))
	} else
	{
		M<-temp[,2:3]
		vec=c(x,t)
		dist<-apply(M,1, function(arg) sum((arg-vec)^2))
		l=which.min(dist)
		return(temp[l,])
	}
}


################## Gaussian noise

evalGauss<-function(epsilon)
{
	return( 1/sqrt(2*pi*sigma2)*exp(-epsilon^2/2/sigma2))
}

################## Link function

Flink<-function(x,fun) 
{
	if (fun=="id")
		return(x)
	if (fun=="plateau")
	{
		if (x<=(x0+2))
			return(x0)
		if (x <= (x0+2.2))
			return (11*x-22-10*x0)
		return(x)				
	}	
}


################## Noise matrix

NoiseMat<-function(epsilon, flink) # computes f(epsilon + F(w^i)-F(w^j))
{
	NMat<-matrix(nrow=nOmega-1,ncol=nOmega-1)
	for (i in 1:(nOmega-1))
		NMat[i,]<-sapply(xval,function(x) evalGauss(epsilon+Flink(xval[i],fun=flink)-sapply(x, Flink, fun=flink)))
	NMat[NMat==0]<-MachineMin	
	return(NMat)
}

################## Computes Pbar*theta 

Ptheta<-function(theta)
{
	Ptheta<-array(dim=c(nOmega-1,nd))
	for (dind in 1:nd)
			Ptheta[,dind]<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%theta # vector (in j) Sum_i P(w^j|w^i,d).theta(w^i), (cemetary not treated here) 
	return(Ptheta)
}

################## Computes Psi

EvalPsi<-function(NM,Ptheta) # returns a matrix: to each line corresponds a filter for a y' associated to a simulated noise associated to a point of Omega. 
{	
	Psi<-array(data=0,dim=c(nOmega-1,nOmega-1,nd))
	for (dind in 1:nd)
		{
			SumNum=Ptheta[,dind] # vector (einn j) Sum_i P(w^j|w^i,d).theta(w^i)
			MatNum=t(t(NM)*SumNum) # matrix (k,j) f(epsilon+F(w^k)-F(w^j))*Sum_i P(w^j|w^i,d).theta(w^i)
			VecDen=as.vector(NM%*%SumNum) # vector (k) Sum_l f(epsilon+F(w^k)-F(w^l))*Sum_i P(w^l|w^i,d).theta(w^i)
			VecDen[VecDen==0]<-1
			Psi[,,dind]=MatNum/VecDen
		}
	return(Psi)
}

############### Nearest neighbor in Gamma

NNeighborY<-function(y)
{
		dist<-(y-GammaY)^2
		l=which.min(dist)
		return(c(GammaY[l],l))
}

NNeighborT<-function(Theta)
{
		if (sum(Theta)==0) return(c(rep(0,nOmega-1),ngamma))
		dist<-apply(Gamma,2, function(arg) sum((arg-Theta)^2))
		l=which.min(dist)
		return (c(Gamma[,l],l))
}


GammaDist<-function(Filt)
{
		F0=Filt[1:nm0]
		F1=Filt[(nm0+1):(nm0+nm1)]
		F2=Filt[(nm0+nm1+1):(nm0+nm1+nm2)]
		f0=sum(F0); f1=sum(F1); f2=sum(F2)
		return(c(Filt,f0,f1,f2))
}


distblock<-function(Filt,Theta)
{
		T0=Theta[1:nm0]
		T1=Theta[(nm0+1):(nm0+nm1)]
		T2=Theta[(nm0+nm1+1):(nm0+nm1+nm2)]
		p0=sum(T0); p1=sum(T1); p2=sum(T2)
		return(abs(p0-Filt[nOmega])+abs(p1-Filt[nOmega+1])+abs(p2-Filt[nOmega+2])+sqrt(sum((Filt[1:(nOmega-1)]-Theta)^2)))	
}

NNeighborBlockT<-function(Theta)
{
		if (sum(Theta)==0) return(c(rep(0,nOmega-1),ngamma))

		dist<-apply(Gdist,2, distblock, Theta=Theta)
		l=which.min(dist)
		return (c(Gamma[,l],l))
}

############### to fill R matrix

RFill<-function(VF,d)
{
	RTemp[VF[1],VF[2],d]<-RTemp[VF[1],VF[2],d]+VF[3]
	return(RTemp)
}



################ Real cost of a trajectory

RCostReal<-function(x,m,d,v,cost) # x marker value, m= real mode, d= treatment, v=next visit time, cost="simple" or "int"
{
	L=length(x)
	if (cost=="simple")
	{
		if (x[L]==D) RC=M else RC=CVs+x[L]-x0
		for (l in 1:(L-1))
		{
			s=c(d[l],v[l])
			RC=RC+CVs+costsimple(m[l]+1,s)
		}
	} else
	{
		if (x[L]==D) RC=MI else RC=CVi+x[L]-x0
		for (l in 1:(L-1))
		{
			futur=(x[l+1]-x0)*as.numeric(v[l])
			t=d[l]
			if (t=="non") tt=1 else if (t=="a") tt=2 else if (t=="b") tt=3 else tt=4		
			pen=tablevalue[m[l]+1,tt]*as.numeric(v[l])
			RC=RC+CVi+pen+futur/cdivide
		}
	}
	return(RC)
}

################ Filtered cost of a trajectory

RCostFilt<-function(indtheta,d,v,cost) # indtheta= vector of elemnt index in grid Gamma, d= treatment, v=next visit time, cost="simple" or "int"
{
	L=length(indtheta)
	if (cost=="simple")
	{
		if (indtheta[L]==ngamma) RC=M else RC=CVs+Finalcosttheta(indtheta[L])
		for (l in 1:(L-1))
		{
			s=c(d[l],v[l])
			RC=RC+CVs+costsimpletheta(indtheta[l],s)
		}
	} else
	{
		if (indtheta[L]==ngamma) RC=M else RC=CVi+FinalcostthetaI(indtheta[L])
		for (l in 1:(L-1))
		{
			s=c(d[l],v[l])
			RC=RC+CVi+costvaluetheta(indtheta[l],s)
		}
	}
	return(RC)
}

################ Computing time spent with wrong treatment

TimeError<-function(m,d,v) #  m= real mode, d= treatment, v=next visit
{
	L=length(m)
	time=0; timewrong=0
	for (l in 1:(L-1))
	{
		s=c(d[l],v[l])
		if (costsimple(m[l]+1,s)!=0)
			time=time+as.numeric(v[l])
		if ((m[l]==1) || (m[l]==2))	
			if (costsimple(m[l]+1,s)!=0)
				timewrong=timewrong+as.numeric(v[l])
		
	}
	return(c(time,timewrong))
}


