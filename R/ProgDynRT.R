ProgDRT<-function()
{

	######################################################  All (see ProgDyn for details)###################################################################


	########## simple cost #################

	Vprogsimple=matrix(ncol=ngamma,nrow=N+1) 
	Decsimple=matrix(ncol=ngamma,nrow=N) 
	Vprogsimple[N+1,]<-sapply(1:ngamma,Finalcosttheta)


	stagecost<-matrix(ncol=ngamma,nrow=nrow(Decisions))
	for (i in 1:nd)
		stagecost[i,]<-sapply(1:ngamma,costsimpletheta,s=Decisions[i,])
	
		
	for (k in (N-1):0)
	{
		if((N-k)<tmax)
		{
			ttemp=N-k
			p=sum(Time/delta<=ttemp)
			ndtemp=p*length(Treatment)
			Matk<-matrix(ncol=ngamma,nrow=ndtemp)
			for (tint in 1:p)
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVs+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogsimple[k+tmult+1,],0)	
				}
			Decsimple[k+1,]=apply(Matk,2,which.min)
			Decsimple[k+1,ngamma]<-10
			Vprogsimple[k+1,]=apply(Matk,2,min)
		} else
		{
			Matk<-matrix(100000,ncol=ngamma,nrow=nd)
			for (tint in 1:length(Time))
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVs+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogsimple[k+tmult+1,],0)	
				}
			Decsimple[k+1,]=apply(Matk,2,which.min)
			Decsimple[k+1,ngamma]<-10
			Vprogsimple[k+1,]=apply(Matk,2,min)
		}	
	}


	########## marker-dependent #################

	Vprogint=matrix(ncol=nT+1,nrow=N+1) 
	Decint=matrix(ncol=nT+1,nrow=N)
	Vprogint[N+1,]<-sapply(1:ngamma,FinalcostthetaI)


	stagecost<-matrix(ncol=nT+1,nrow=nrow(Decisions))
	for (d in 1:nd)
		stagecost[d,]<-sapply(1:ngamma,costvaluetheta,s=Decisions[d,])
	
		
	for (k in (N-1):0)
	{
		if((N-k)<tmax)
		{
			ttemp=N-k
			p=sum(Time/delta<=ttemp)
			ndtemp=p*length(Treatment)
			Matk<-matrix(ncol=nT+1,nrow=ndtemp)
			for (tint in 1:p)
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVi+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogint[k+tmult+1,],0)	
				}		
			Decint[k+1,]=apply(Matk,2,which.min)
			Decint[k+1,ngamma]<-10
			Vprogint[k+1,]=apply(Matk,2,min)
		} else
		{
			Matk<-matrix(100000,ncol=ngamma,nrow=nd)
			for (tint in 1:length(Time))
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVi+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogint[k+tmult+1,],0)	
				}
			Decint[k+1,]=apply(Matk,2,which.min)
			Decint[k+1,ngamma]<-10
			Vprogint[k+1,]=apply(Matk,2,min)
		}	
	}

	###################################################"

	save(Decsimple,Vprogsimple,Decint,Vprogint,file=paste('RData/ProgDynRT_', nT,"-", distance,'_timeall_CV', CVs,'.RData', sep=""))



	################################################################# Time 15 ########################################################


	########## simple cost #################

	Vprogsimple=matrix(ncol=ngamma,nrow=N+1) 
	Decsimple=matrix(ncol=ngamma,nrow=N)
	Vprogsimple[N+1,]<-sapply(1:ngamma,Finalcosttheta)


	stagecost<-matrix(ncol=ngamma,nrow=nrow(Decisions))
	for (i in 1:nd)
		stagecost[i,]<-sapply(1:ngamma,costsimpletheta,s=Decisions[i,])
	
		
	for (k in (N-1):0)
	{
		Matk<-matrix(100000,ncol=ngamma,nrow=nd)
		for (tint in 1)
			for (d in 1:3)
			{
				tmult=(Time/delta)[tint]
				inddec=(tint-1)*3+d
				Matk[inddec,]<-CVs+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogsimple[k+tmult+1,],0)	
			}
		Decsimple[k+1,]=apply(Matk,2,which.min)
		Decsimple[k+1,ngamma]<-10
		Vprogsimple[k+1,]=apply(Matk,2,min)
		
	}


	########## marker-dependent #################

	Vprogint=matrix(ncol=nT+1,nrow=N+1) 
	Decint=matrix(ncol=nT+1,nrow=N)
	Vprogint[N+1,]<-sapply(1:ngamma,FinalcostthetaI)


	stagecost<-matrix(ncol=ngamma,nrow=nrow(Decisions))
	for (d in 1:nd)
		stagecost[d,]<-sapply(1:ngamma,costvaluetheta,s=Decisions[d,])
	
		
	for (k in (N-1):0)
	{
			Matk<-matrix(100000,ncol=ngamma,nrow=nd)
			for (tint in 1)
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVi+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogint[k+tmult+1,],0)	
				}
			Decint[k+1,]=apply(Matk,2,which.min)
			Decint[k+1,ngamma]<-10
			Vprogint[k+1,]=apply(Matk,2,min)

	}

	###################################################"

	save(Decsimple,Vprogsimple,Decint,Vprogint,file=paste('RData/ProgDynRT_', nT,"-", distance,'_time15_CV', CVs,'.RData', sep=""))



	####################################################  Time 60 #############################################################################"


	##########  simple #################

	Vprogsimple=matrix(ncol=ngamma,nrow=N+1) 
	Decsimple=matrix(ncol=ngamma,nrow=N)
	Vprogsimple[N+1,]<-sapply(1:ngamma,Finalcosttheta)


	stagecost<-matrix(ncol=ngamma,nrow=nrow(Decisions))
	for (i in 1:nd)
		stagecost[i,]<-sapply(1:ngamma,costsimpletheta,s=Decisions[i,])
	
		
	for (k in (N-1):0)
	{
		if((N-k)<tmax)
		{
			ttemp=N-k
			p=sum(Time/delta<=ttemp)
			ndtemp=p*length(Treatment)
			Matk<-matrix(ncol=ngamma,nrow=ndtemp)
			for (tint in 1:p)
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVs+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogsimple[k+tmult+1,],0)	
				}			
			Decsimple[k+1,]=apply(Matk,2,which.min)
			Decsimple[k+1,ngamma]<-10
			Vprogsimple[k+1,]=apply(Matk,2,min)
		} else
		{
			Matk<-matrix(100000,ncol=ngamma,nrow=nd)
			for (tint in length(Time))
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVs+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogsimple[k+tmult+1,],0)	
				}
			Decsimple[k+1,]=apply(Matk,2,which.min)
			Decsimple[k+1,ngamma]<-10
			Vprogsimple[k+1,]=apply(Matk,2,min)
		}	
	}


	########## marker-dependent #################

	Vprogint=matrix(ncol=nT+1,nrow=N+1) 
	Decint=matrix(ncol=nT+1,nrow=N)
	Vprogint[N+1,]<-sapply(1:ngamma,FinalcostthetaI)


	stagecost<-matrix(ncol=nT+1,nrow=nrow(Decisions))
	for (d in 1:nd)
		stagecost[d,]<-sapply(1:ngamma,costvaluetheta,s=Decisions[d,])
	
		
	for (k in (N-1):0)
	{
		if((N-k)<tmax)
		{
			ttemp=N-k
			p=sum(Time/delta<=ttemp)
			ndtemp=p*length(Treatment)
			Matk<-matrix(ncol=nT+1,nrow=ndtemp)
			for (tint in 1:p)
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVi+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogint[k+tmult+1,],0)	
				}		
			Decint[k+1,]=apply(Matk,2,which.min)
			Decint[k+1,ngamma]<-10
			Vprogint[k+1,]=apply(Matk,2,min)
		} else
		{
			Matk<-matrix(100000,ncol=ngamma,nrow=nd)
			for (tint in length(Time))
				for (d in 1:3)
				{
					tmult=(Time/delta)[tint]
					inddec=(tint-1)*3+d
					Matk[inddec,]<-CVi+stagecost[inddec,]+c(RT[,,inddec]%*%Vprogint[k+tmult+1,],0)	
				}
			Decint[k+1,]=apply(Matk,2,which.min)
			Decint[k+1,ngamma]<-10
			Vprogint[k+1,]=apply(Matk,2,min)
		}	
	}

	###################################################"


	save(Decsimple,Vprogsimple,Decint,Vprogint,file=paste('RData/ProgDynRT_', nT,"-", distance,'_time60_CV', CVs,'.RData', sep=""))

}

