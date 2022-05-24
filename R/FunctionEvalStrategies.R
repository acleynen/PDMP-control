# evaluation of our strategy on nsim simulations when using strategy "step" ("15" for fixed 15 days, "60" for fixed 60 days or "all", for allowing choice), using distance "distance" ("L2" or "Lm") returning confidence intervals or not (std), with precision of round digits, using time-dependent (simple) or marker-dependent (int) cost, and fixing the seed.

OurStrat<-function(nsim,step,distance,std=TRUE,round=2,cost="simple",seed=212121)
{

	RealCosts=c();RealCosti=c();FilterCosts=c();FilterCosti=c();RFCosts=c();RFCosti=c();Lengths=c();Lengthi=c();RealLengths=c();RealLengthi=c();Nds=0;Ndi=0;TimeErr=c();TimeWrong=c();
	# RealCosts and RealCosti : average cost of trajectories evaluated on the true process $X$ 
	# FilterCosts and FilterCosti : average cost evaluated using the projected filters
	# RFCosts and RFCosti : average cost evaluated using the estimated filters (before projection)
	#	Lengths and Lengthi : average length of a trajectory (for step=15 and 60, this might vary on average due to trajectories leading to death. for "all", this length varies also due to possible choice of next visit)
	# Nds and Ndi : number of death observed
	# TimeErr = average time patient spends with the wrong decision
	# TimeWrong = average time patient spends with the wrong treatment (when he needs one)
	
	if (cost=="simple")
	{
		for (sim in 1:nsim)
		{
			set.seed(seed+sim)
			xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
			Strat=list(dec=rep("non",N),vis=rep(0,N))
			yn=0
			Psin=as.matrix(c(1,rep(0,nOmega-2)))
			decision<-Decsimple[1,1]
			Strat$dec[1]<-Decisions[decision,1]
			Strat$vis[1]<-Decisions[decision,2]		
			RealTime=0
			n=1
			IndexGamma<-c(1)
			while(RealTime<H)
			{
				xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
				xn$traj[n+1]=xnew[2]
				xn$time[n+1]=xnew[3]
				xn$mode[n+1]=xnew[1]
				xn$realtime[n+1]=xnew[4]
				RealTime<-RealTime+as.numeric(Strat$vis[n])
	
				if(xnew[1]==100)
				{
					ynew=0 ; Psin<-cbind(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H; IndexGamma=c(IndexGamma,ngamma);
				} else
				{	
					vmult=which(Time==Strat$vis[n])
					dmult=which(Treatment==Strat$dec[n])
					dind=(vmult-1)*3+dmult
					ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
					if (trunc) 
						if (ynew<2) ynew=0
					PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
					eG<-evalGauss(ynew-xval)
					Num=PTemp*eG
					Den=eG%*%PTemp
					Psitemp=Num/as.vector(Den)
					if (distance=="L2")	NT=NNeighborT(Psitemp)[nOmega] else		NT=NNeighborBlockT(Psitemp)[nOmega] 
					Psin<-cbind(Psin,Psitemp)
					indGamma=NT
					IndexGamma<-c(IndexGamma,indGamma)
					indvisit=RealTime/delta
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
					yn=c(yn,ynew)
				}
				n<-n+1
			}

			ll=which(xn$realtime==H)
			if (sum(xn$mode==100)>0) {ll=which(xn$mode==100); Nds=Nds+1}
			rl=xn$realtime[ll]

			RealCosts=c(RealCosts,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
			FilterCosts=c(FilterCosts,RCostFilt(IndexGamma,Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
			RFCosts=c(RFCosts,RCostFiltTraj(Psin,Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
			Lengths=c(Lengths,ll)
			RealLengths=c(RealLengths,rl)	
			TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
			TimeErr=c(TimeErr,TT[1])
			TimeWrong=c(TimeWrong,TT[2])
		}
		if (std)
			return(c(round(Vprogsimple[1,1],round), round(mean(RealCosts),round), paste('(', round(1.96*sd(RealCosts)/sqrt(nsim),round),')',sep=""), round(mean(FilterCosts),round), paste('(',round(1.96*sd(FilterCosts)/sqrt(nsim),round),')',sep=""), round(mean(RFCosts),round), paste('(', round(1.96*sd(RFCosts)/sqrt(nsim),round),')',sep=""),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Lengths),2),round(mean(RealLengths),2),Nds))
		return(c(round(Vprogsimple[1,1],round), round(mean(RealCosts),round),round(mean(RFCosts),round),round(mean(FilterCosts),round),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Lengths),2),round(mean(RealLengths),2),Nds))		
	}
	
	for (sim in 1:nsim)
	{
		set.seed(seed+sim)
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		Strat=list(dec=rep("non",N),vis=rep(0,N))
		yn=0
		Psin=as.matrix(c(1,rep(0,nOmega-2)))
		decision<-Decint[1,1]
		Strat$dec[1]<-Decisions[decision,1]
		Strat$vis[1]<-Decisions[decision,2]		
		RealTime=0
		n=1
		IndexGamma<-c(1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psin<-cbind(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H; IndexGamma=c(IndexGamma,ngamma);
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc)
					if (ynew<2) ynew=0				
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				if (distance=="L2")	NT=NNeighborT(Psitemp)[nOmega] else		NT=NNeighborBlockT(Psitemp)[nOmega] 
				Psin<-cbind(Psin,Psitemp)
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				indvisit=RealTime/delta
				if (RealTime<H)
				{
					decision<-Decint[indvisit+1,indGamma]
					Strat$dec[n+1]<-Decisions[decision,1]
					Strat$vis[n+1]<-Decisions[decision,2]
				} else 
				{
					Strat$dec[n+1]<-nd+1	
				}	
				yn=c(yn,ynew)
			}
			n<-n+1
		}

		ll=which(xn$realtime==H)
		if (sum(xn$mode==100)>0) 
			{ll=which(xn$mode==100); Ndi=Ndi+1}
		rl=xn$realtime[ll]

		RealCosti=c(RealCosti,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		FilterCosti=c(FilterCosti,RCostFilt(IndexGamma,Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		RFCosti=c(RFCosti,RCostFiltTraj(Psin,Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		Lengthi=c(Lengthi,ll)
		RealLengthi=c(RealLengthi,rl)
		TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
		TimeErr=c(TimeErr,TT[1])
		TimeWrong=c(TimeWrong,TT[2])
	}
	if (std)
		return(c(round(Vprogint[1,1],round),round(mean(RealCosti),round),paste('(',round(1.96*sd(RealCosti)/sqrt(nsim),round),')',sep=""),round(mean(FilterCosti),round),paste('(',round(1.96*sd(FilterCosti)/sqrt(nsim),round),')',sep=""),round(mean(RFCosti),round),paste('(',round(1.96*sd(RFCosti)/sqrt(nsim),round),')',sep="") ,round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Lengthi),2), round(mean(RealLengthi),2),Ndi))
	return(c( round(Vprogint[1,1],2),round(mean(RealCosti),round),round(mean(FilterCosti),round),round(mean(RFCosti),round),round(mean(TimeErr),2),round(mean(TimeWrong),2), round(mean(Lengthi),2), round(mean(RealLengthi),2),Ndi))
}





##################  "See all" strategy

# Evaluation of cost of the strategies where decisions are taken knowing true process $X$. Hence optimal decisions are to not treat unless the mode is non-negative and treat with appropriate treatment when mode has jumped, all the while returning every fixed "step" days.

SeeStrat<-function(nsim,step,std=TRUE,round=2,seed=212121)
{

	RealCosts=c();RealCosti=c();Length=c();TimeErr=c();TimeWrong=c();

	for (sim in 1:nsim)
	{
		set.seed(seed+sim)
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ;  xn$realtime[n+2]=H; xn$traj[n+2]=D;  RealTime<-H
			} else if (RealTime==H)
			{
				Strat$dec[n+1]<-nd+1	
			} else
			{	
				if (xnew[1]==0)
				{
					Strat$dec[n+1]<-Decisions[1,1]
				}
				if (xnew[1]==1)
				{
					Strat$dec[n+1]<-Decisions[2,1]
				}			
				if (xnew[1]==2)
				{
					Strat$dec[n+1]<-Decisions[3,1]
				}			
				Strat$vis[n+1]<-step
			}
			n<-n+1
		}

		ll=which(xn$realtime==H)
		if (sum(xn$mode==100)>0) ll=which(xn$mode==100)


		RealCosts=c(RealCosts,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
		RealCosti=c(RealCosti,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		Length=c(Length,ll)
		TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
		TimeErr=c(TimeErr,TT[1])
		TimeWrong=c(TimeWrong,TT[2])			
	}
	if (std)
		return(c(round(mean(RealCosts),round), paste("(",round(1.96*sd(RealCosts)/sqrt(nsim),round),')',sep=""),round(mean(RealCosti),round), paste('(',round(1.96*sd(RealCosti)/sqrt(nsim),round),')',sep=""),round(mean(TimeErr),2),round(mean(TimeWrong),2), round(mean(Length),2)))
	return(c(round(mean(RealCosts),round),round(mean(RealCosti),round),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Length),2)))

}





############ Kalman-like strategy

# Evaluation of cost of the strategies where decisions are taken estimating the filter and using the mode with maximal mass as ground truth. Hence optimal decisions are to not treat unless the maximal mode is non-negative and treat with appropriate treatment when maximal mode is positive. This can be done by fixing the return dates to 15 or 60 days.


ModeStrat<-function(nsim,step,std=TRUE,round=2,seed=212121)
{

	RealCosts=c();RealCosti=c();RFCosts=c();RFCosti=c();Length=c();TimeErr=c();TimeWrong=c();

	for (sim in 1:nsim)
	{
		set.seed(seed+sim)
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
	Psin=as.matrix(c(1,rep(0,nOmega-2)))
		if (step==60) Strat=list(dec=rep("non",N),vis=c(60,rep(0,N-1))) else Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psin<-cbind(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D;  RealTime<-H
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc) { if (ynew<2) ynew=0 }
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				Psin<-cbind(Psin,Psitemp)
				Mprob<-c(sum(Psitemp[1:nm0]),sum(Psitemp[(nm0+1):(nm0+nm1)]),sum(Psitemp[(nm0+nm1+1):(nm0+nm1+nm2)]))
				Mmax=which.max(Mprob) 
				if (RealTime<H)
				{
					Strat$dec[n+1]<-Decisions[Mmax,1]
					Strat$vis[n+1]<-step
				} else 
				{
					Strat$dec[n+1]<-nd+1	
				}	
				yn=c(yn,ynew)
			}
			n<-n+1
		}

		ll=which(xn$realtime==H)
		if (sum(xn$mode==100)>0) ll=which(xn$mode==100)


		RealCosts=c(RealCosts,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
		RFCosts=c(RFCosts,RCostFiltTraj(Psin,Strat$dec[1:ll],Strat$vis[1:ll],"simple"))	
		RealCosti=c(RealCosti,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		RFCosti=c(RFCosti,RCostFiltTraj(Psin,Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		Length=c(Length,ll)
		TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
		TimeErr=c(TimeErr,TT[1])
		TimeWrong=c(TimeWrong,TT[2])		
	}
	if (std)
		return( c(NA,round(mean(RealCosts),round), paste('(', round(1.96*sd(RealCosts)/sqrt(nsim),round),')',sep=""), NA, NA, round(mean(RFCosts),round), paste('(',round(1.96*sd(RFCosts)/sqrt(nsim),round),')',sep=""), NA, round(mean(RealCosti),round),paste('(',round(1.96*sd(RealCosti)/sqrt(nsim),round),')',sep=""),NA,NA,round(mean(RFCosti),round),paste('(',round(1.96*sd(RFCosti)/sqrt(nsim),round),')',sep=""),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Length),2)))
	return(c(NA, round(mean(RealCosts),round),NA,round(mean(RFCosts),round), NA, round(mean(RealCosti),round),NA,round(mean(RFCosti),round),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Length),2)))

}




#########   Classical strategy

s1=5 # threshold for relapse
srem=1 # threshold for remission
Treatm=c("a","b")

# classic strategies use fixed thresholds to detect relapse and remission. Treatments are then adapted by always trying treatment $a$ first, and changing treatment if marker value increases. Return dates are then adapted: 15 days when relapse is detected, 60 days when remission is assumed.

RealStrat<-function(nsim,std=TRUE,round=2,seed=212121)
{

	RealCosts=c();RealCosti=c();Length=c();TimeErr=c();TimeWrong=c();

	for (sim in 1:nsim)
	{
		set.seed(seed+sim)
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		yn=0
		Psin=c(1,rep(0,N-1))
		Strat=list(dec=rep("non",N),vis=c(delta,rep(0,N-1)))
		RealTime=0
		n=1
		IndexGamma<-c(1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psinewproj=ngamma; xn$realtime[n+2]=H; xn$traj[n+2]=D;  RealTime<-H
			} else
			{	
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc)
					if (ynew<2) ynew=0				
				if (RealTime==H) 
				{
					Strat$dec[n+1]<-nd+1	
				} else
				{
					if  (Strat$dec[n]=="non")
					{
						if (ynew<s1)
						{
							Strat$dec[n+1]<-"non"
							if ((H-RealTime)<60) Strat$vis[n+1]<-(H-RealTime) else 		Strat$vis[n+1]<-60				
						}	else
						{
							Strat$dec[n+1]<-"b"
							Strat$vis[n+1]<-15
							}
					} else 	
					{
						if (ynew>yn[n])
						{
							Strat$dec[n+1]<-Treatm[!is.element(Treatm,Strat$dec[n])]			
						} else if (ynew<srem)
						{
							Strat$dec[n+1]<-"non"
						} else
						{
							Strat$dec[n+1]<-Strat$dec[n]
						}
						Strat$vis[n+1]<-15
					} 
				}
				yn=c(yn,ynew)
			}
			n<-n+1
		}

		ll=which(xn$realtime==H)
		if (sum(xn$mode==100)>0) ll=which(xn$mode==100)


		RealCosts=c(RealCosts,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
		RealCosti=c(RealCosti,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		Length=c(Length,ll)
		TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
		TimeErr=c(TimeErr,TT[1])
		TimeWrong=c(TimeWrong,TT[2])				
	}
	if (std)
		return(c(NA,round(mean(RealCosts),round), paste("(",round(1.96*sd(RealCosts)/sqrt(nsim),round),')',sep=""),NA,NA,NA, NA, NA, round(mean(RealCosti),round), paste('(',round(1.96*sd(RealCosti)/sqrt(nsim),round),')',sep=""),NA,NA,NA,NA,round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Length),2)))
	return(c(NA,round(mean(RealCosts),round),NA,NA,NA,round(mean(RealCosti),round),NA,NA,round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Length),2)))

}

#################################"  Using projected filter as current filter


OurStratchap<-function(nsim,step,distance,std=TRUE,round=2,cost="simple",seed=212121)
{

	RealCosts=c();RealCosti=c();FilterCosts=c();FilterCosti=c();RFCosts=c();RFCosti=c();Lengths=c();Lengthi=c();RealLengths=c();RealLengthi=c();Nds=0;Ndi=0;TimeErr=c();TimeWrong=c();
	# RealCosts and RealCosti : average cost of trajectories evaluated on the true process $X$ 
	# FilterCosts and FilterCosti : average cost evaluated using the projected filters
	# RFCosts and RFCosti : average cost evaluated using the estimated filters (before projection)
	#	Lengths and Lengthi : average length of a trajectory (for step=15 and 60, this might vary on average due to trajectories leading to death. for "all", this length varies also due to possible choice of next visit)
	# Nds and Ndi : number of death observed
	# TimeErr = average time patient spends with the wrong decision
	# TimeWrong = average time patient spends with the wrong treatment (when he needs one)
	
	if (cost=="simple")
	{
		for (sim in 1:nsim)
		{
			set.seed(seed+sim)
			xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
			Strat=list(dec=rep("non",N),vis=rep(0,N))
			yn=0
			Psin=as.matrix(c(1,rep(0,nOmega-2)))
			decision<-Decsimple[1,1]
			Strat$dec[1]<-Decisions[decision,1]
			Strat$vis[1]<-Decisions[decision,2]		
			RealTime=0
			n=1
			IndexGamma<-c(1)
			while(RealTime<H)
			{
				xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
				xn$traj[n+1]=xnew[2]
				xn$time[n+1]=xnew[3]
				xn$mode[n+1]=xnew[1]
				xn$realtime[n+1]=xnew[4]
				RealTime<-RealTime+as.numeric(Strat$vis[n])
	
				if(xnew[1]==100)
				{
					ynew=0 ; Psin<-cbind(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H; IndexGamma=c(IndexGamma,ngamma);
				} else
				{	
					vmult=which(Time==Strat$vis[n])
					dmult=which(Treatment==Strat$dec[n])
					dind=(vmult-1)*3+dmult
					ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
					if (trunc) 
						if (ynew<2) ynew=0
					PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
					eG<-evalGauss(ynew-xval)
					Num=PTemp*eG
					Den=eG%*%PTemp
					Psitemp=Num/as.vector(Den)
					if (distance=="L2")	NT=NNeighborT(Psitemp)[nOmega] else		NT=NNeighborBlockT(Psitemp)[nOmega] 
					Psin<-cbind(Psin,Gamma[,NT])
					indGamma=NT
					IndexGamma<-c(IndexGamma,indGamma)
					indvisit=RealTime/delta
					if (RealTime<H)
					{
						decision<-Decsimple[indvisit+1,indGamma]
						Strat$dec[n+1]<-Decisions[decision,1]
						Strat$vis[n+1]<-Decisions[decision,2]
					} else 
					{
						Strat$dec[n+1]<-nd+1	
					}	
					yn=c(yn,ynew)
				}
				n<-n+1
			}

			ll=which(xn$realtime==H)
			if (sum(xn$mode==100)>0) {ll=which(xn$mode==100); Nds=Nds+1}
			rl=xn$realtime[ll]

			RealCosts=c(RealCosts,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
			FilterCosts=c(FilterCosts,RCostFilt(IndexGamma,Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
			RFCosts=c(RFCosts,RCostFiltTraj(Psin,Strat$dec[1:ll],Strat$vis[1:ll],"simple"))
			Lengths=c(Lengths,ll)
			RealLengths=c(RealLengths,rl)	
			TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
			TimeErr=c(TimeErr,TT[1])
			TimeWrong=c(TimeWrong,TT[2])
		}
		if (std)
			return(c(round(Vprogsimple[1,1],round), round(mean(RealCosts),round), paste('(', round(1.96*sd(RealCosts)/sqrt(nsim),round),')',sep=""), round(mean(FilterCosts),round), paste('(',round(1.96*sd(FilterCosts)/sqrt(nsim),round),')',sep=""), round(mean(RFCosts),round), paste('(', round(1.96*sd(RFCosts)/sqrt(nsim),round),')',sep=""),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Lengths),2),round(mean(RealLengths),2),Nds))
		return(c(round(Vprogsimple[1,1],round), round(mean(RealCosts),round),round(mean(RFCosts),round),round(mean(FilterCosts),round),round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Lengths),2),round(mean(RealLengths),2),Nds))		
	}
	
	for (sim in 1:nsim)
	{
		set.seed(seed+sim)
		xn=list(mode=rep(0,N),traj=rep(x0,N),time=rep(0,N),realtime=rep(0,N))
		Strat=list(dec=rep("non",N),vis=rep(0,N))
		yn=0
		Psin=as.matrix(c(1,rep(0,nOmega-2)))
		decision<-Decint[1,1]
		Strat$dec[1]<-Decisions[decision,1]
		Strat$vis[1]<-Decisions[decision,2]		
		RealTime=0
		n=1
		IndexGamma<-c(1)
		while(RealTime<H)
		{
			xnew=Step(1,xn$mode[n],xn$traj[n],xn$time[n],RealTime,as.numeric(Strat$vis[n]),Strat$dec[n])
			xn$traj[n+1]=xnew[2]
			xn$time[n+1]=xnew[3]
			xn$mode[n+1]=xnew[1]
			xn$realtime[n+1]=xnew[4]
			RealTime<-RealTime+as.numeric(Strat$vis[n])
	
			if(xnew[1]==100)
			{
				ynew=0 ; Psin<-cbind(Psin,rep(0,nOmega-1)); xn$realtime[n+2]=H; xn$traj[n+2]=D; RealTime<-H; IndexGamma=c(IndexGamma,ngamma);
			} else
			{	
				vmult=which(Time==Strat$vis[n])
				dmult=which(Treatment==Strat$dec[n])
				dind=(vmult-1)*3+dmult
				ynew=max(Flink(xnew[2],flink)+rnorm(1,mean=0,sd=sqrt(sigma2)),0)
				if (trunc)
					if (ynew<2) ynew=0				
				PTemp<-t(P[1:(nOmega-1),1:(nOmega-1),dind])%*%Psin[,n]
				eG<-evalGauss(ynew-xval)
				Num=PTemp*eG
				Den=eG%*%PTemp
				Psitemp=Num/as.vector(Den)
				if (distance=="L2")	NT=NNeighborT(Psitemp)[nOmega] else		NT=NNeighborBlockT(Psitemp)[nOmega] 
				Psin<-cbind(Psin,Gamma[,NT])
				indGamma=NT
				IndexGamma<-c(IndexGamma,indGamma)
				indvisit=RealTime/delta
				if (RealTime<H)
				{
					decision<-Decint[indvisit+1,indGamma]
					Strat$dec[n+1]<-Decisions[decision,1]
					Strat$vis[n+1]<-Decisions[decision,2]
				} else 
				{
					Strat$dec[n+1]<-nd+1	
				}	
				yn=c(yn,ynew)
			}
			n<-n+1
		}

		ll=which(xn$realtime==H)
		if (sum(xn$mode==100)>0) 
			{ll=which(xn$mode==100); Ndi=Ndi+1}
		rl=xn$realtime[ll]

		RealCosti=c(RealCosti,RCostReal(xn$traj[1:ll],xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		FilterCosti=c(FilterCosti,RCostFilt(IndexGamma,Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		RFCosti=c(RFCosti,RCostFiltTraj(Psin,Strat$dec[1:ll],Strat$vis[1:ll],"int"))
		Lengthi=c(Lengthi,ll)
		RealLengthi=c(RealLengthi,rl)
		TT=TimeError(xn$mode[1:ll],Strat$dec[1:ll],Strat$vis[1:ll])
		TimeErr=c(TimeErr,TT[1])
		TimeWrong=c(TimeWrong,TT[2])
	}
	if (std)
		return(c(round(Vprogint[1,1],round),round(mean(RealCosti),round),paste('(',round(1.96*sd(RealCosti)/sqrt(nsim),round),')',sep=""),round(mean(FilterCosti),round),paste('(',round(1.96*sd(FilterCosti)/sqrt(nsim),round),')',sep=""),round(mean(RFCosti),round),paste('(',round(1.96*sd(RFCosti)/sqrt(nsim),round),')',sep="") ,round(mean(TimeErr),2),round(mean(TimeWrong),2),round(mean(Lengthi),2), round(mean(RealLengthi),2),Ndi))
	return(c( round(Vprogint[1,1],2),round(mean(RealCosti),round),round(mean(FilterCosti),round),round(mean(RFCosti),round),round(mean(TimeErr),2),round(mean(TimeWrong),2), round(mean(Lengthi),2), round(mean(RealLengthi),2),Ndi))
}


