


#Code requires plotting functions at the end of this script file, run those first
# heatmap FUN 
SpatiotempPlot=function(spacetime, Xax=1:dim(spacetime)[2], Yax=1:dim(spacetime)[1], XaxN="X", YaxN="Y", figtitle="Title", 
					Zlim=c(0, max(spacetime, na.rm=TRUE)), cont=NULL, cexAx=1, contPlot=spacetime, cexCont=1.5*cexAx, Conts=NULL, contSpan=1, palette=1){
	require(fields); spacetime[is.na(spacetime)] = Zlim[1] - 1;
	COL = rev(rainbow(1000, start = 0, end = 0.7))
	if (palette > 1) { require(pals); COL = list(parula(1000), head(tail(parula(1000), -50, -50)))[[palette - 1]]; }
	if (length(Zlim) == 1) {
		Zlim = quantile(spacetime, c(Zlim, 1 - Zlim), na.rm = T)
		Rng = range(spacetime)
		if (Zlim[1] < Rng[1]) Zlim[1] = Rng[1]
		if (Zlim[2] > Rng[2]) Zlim[2] = Rng[2]
	}
	spacetime[which(is.na(spacetime), arr.ind = TRUE)] = max(Zlim) + 1
	image.plot(x = Xax, y = Yax, z = t(spacetime), zlim = Zlim, 
			 xlab = XaxN, ylab = YaxN, cex.axis = cexAx, cex.lab = cexAx, 
			 legend.cex = cexAx, main = figtitle, col = COL); box();
	if (!is.null(cont)) {
		if (abs(log(max(contPlot, na.rm = T), 10)) > 2) options(scipen = 0)
		if (contSpan > 1) {
			smoo1 = t(apply(contPlot, 1, function(x) supsmu(1:ncol(contPlot), x, span = contSpan/ncol(contPlot))$y))
			smoo2 = t(apply(smoo1, 1, function(x) supsmu(1:ncol(smoo1), x, span = contSpan/ncol(smoo1))$y))
			contPlot = smoo2
		}
		if (is.null(Conts)) contour(x = Xax, y = Yax, z = t(contPlot), add = T, col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
		if (!is.null(Conts)) contour(x = Xax, y = Yax, z = t(contPlot), levels = Conts, add = T, col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
		if (abs(log(max(contPlot, na.rm = T), 10)) > 2) options(scipen = 0)
	}
}
gampred=function(dat,varSeq=NULL,varResp=NULL,k=5:30,fam=gaussian,wt=NULL,seqRng=NULL,GAM=NULL, plotgive=TRUE,plotfun=plot,seFact=1.96,spanSignif=0.1,seFill=0.2,seCol=NULL,Lcol=NULL,predType="link",...){
	#assumes that either (a) GAM is built on variables named in dat or (b) gam is unbuilt
	if("matrix"%in%class(dat)) dat=data.frame(dat);
	require(mgcv); xvr=varSeq; yvr=varResp;
	if(!is.null(GAM)){ k="k given"
	if(is.numeric(varSeq)) varResp=colnames(dat)[varSeq]; if(is.numeric(varResp)) varResp=colnames(dat)[varResp];
	if(is.null(varSeq)) varSeq=all.vars(formula(GAM))[2]; if(is.null(varResp)) varResp=all.vars(formula(GAM))[1];
	} else {
		if(is.numeric(varSeq)) colnames(dat)[varSeq]="X" else colnames(dat)[colnames(dat)==varSeq]="X"; varSeq="X";
		if(is.numeric(varResp)) colnames(dat)[varResp]="Y" else colnames(dat)[colnames(dat)==varResp]="Y"; varResp="Y";
		if(length(k)){ Bs=1[-1]; for(ki in k) Bs=c(Bs,tryCatch(BIC(gam(Y~s(X,k=ki),data=dat)),error=function(e) NA)); k=k[c(tail(1+which(diff(Bs)<(-2)),1),1)[1]]; }
		if(length(wt)==1) wt=dat[,wt]; GAM=gam(Y~s(X,k=k), data=dat,weights=wt,family=fam);
	}
	
	if(is.null(seqRng)) seqRng=quantile(dat[,varSeq],c(0.01,0.99),na.rm=TRUE);
	N=1e2; datNew=data.frame(dat[rep(c(which(!is.na(dat))[1],1)[1],N),]);
	varNum=unlist(lapply(dat,class))=="numeric"; datNew[,varNum]=colMeans(t(t(dat[,varNum])),na.rm=TRUE); 
	xseq=seq(seqRng[1],seqRng[2],len=N); datNew[,varSeq]=xseq;
	
	PRED=predict.gam(GAM,newdata=datNew,se.fit=TRUE,type=predType); Bnd=cbind(c(xseq,rev(xseq)),c(PRED$fit+seFact*PRED$se.fit,rev(PRED$fit-seFact*PRED$se.fit)));
	rsq=tryCatch(cor(GAM$fitted.values,GAM$y)^2,error=function(e) NA); pval=summary(GAM)$s.table[,"p-value"][1]; dSignif=1+0*xseq;
	if(!is.na(spanSignif)){
		require(tsgam); #remotes::install_github("gavinsimpson/tsgam") #tsgam not available on CRAN
		Deriv=fderiv(GAM,newdata=datNew); muDeriv=Deriv$derivatives[[varSeq]]$deriv; se2Deriv=1.96*Deriv$derivatives[[varSeq]]$se.deriv;
		dSignif[(muDeriv+se2Deriv)<0]=4; dSignif[(muDeriv-se2Deriv)>0]=2;
		dSrle=rle(dSignif); dSignif[rep(dSrle$lengths,dSrle$lengths)<(spanSignif*N)]=1;
	}
	OUT=list(mu=cbind(xseq,PRED$fit),  polyse=Bnd,   pts=cbind(dat[,varSeq],dat[,varResp]),  dSignif=dSignif,  rsq=rsq,  pval=pval,  k=k); if(!plotgive) return(OUT);
	
	xpt=list(...); chkin=(!c("main","ylab","xlab","col")%in%names(xpt)); rgbCol=c(0,0,0); if(!chkin[4]) rgbCol=col2rgb(xpt$col)[,1]/255;
	if(chkin[2]) xpt$ylab=yvr; if(chkin[3]) xpt$xlab=xvr; if(chkin[1]) xpt$main=paste0(c("R2=",round(rsq,3),"__Pval=",round(pval,3),"__k=",k),collapse="");
	if(!is.null(seCol)) seCol=col2rgb(seCol)/255 else seCol=rgbCol;
	do.call(plotfun, c(list(OUT[["pts"]]), xpt)); polygon(OUT[["polyse"]], border=0, col=do.call(rgb, c(as.list(seCol),alpha=seFill))); 
	if(!is.na(spanSignif)){ segs=c(0,cumsum(diff(dSignif)!=0)); for(si in unique(segs)) lines(OUT[["mu"]][segs==si,],lwd=4,col=OUT[["dSignif"]][segs==si][1]); do.call(points, c(list(OUT[["pts"]]), xpt));
	} else {
		if(is.null(Lcol)) Lcol=do.call(rgb, as.list(rgbCol))
		lines(OUT[["mu"]],lwd=4,col=Lcol)
	}
}




library(deSolve); library(abind); library(nloptr); library(mgcv); 
splinefun2=function(dat,k=20,cats=1:3){
	require(mgcv); datDF=data.frame(dat); datnm=names(datDF);
	xgam=gam(as.formula(paste(datnm[1], "~", "s(", datnm[cats[2]], ",", datnm[cats[3]], ",k=", k, ")")), data = datDF)
	function(a,b){ newd=data.frame(a,b); names(newd)=names(xgam$var.summary);
	predict.gam(xgam, newdata = newd, type = "link"); }
}
sig=function(x,thresh,smoo=15) 1-1/(1+exp(smoo*(x-thresh)))
assdo2=function(parms,varLims=cbind(k=c(0.1,0.5),piA=c(0,1)),Trun=1e2,give=0,ntests=40,ic=NULL,XN=colnames(varLims)[2],YN=colnames(varLims)[1],plotgive=1){
	NP=ncol(varLims); N=ntests^NP; varnms=colnames(varLims); vars=round(apply(varLims,2,function(x) seq(x[1],x[2],len=ntests)),4);
	Pmat=do.call("rbind",rep(list(parms),2*N)); Pmat[,varnms]=apply(as.matrix(expand.grid(vars[,1],vars[,2])),2,rep,2)[,1:NP];
	ic0=c(1,0.05,0.75,0.75,0.5,0.5,  0.1,0.75,0.05,0.05,0.3,0.01); NV=length(ic0)/2; imitig=c(TRUE,rep(FALSE,NV-1));
	if(length(ic)==1){ Pmat[(N+1):(2*N),"bt"]=ic; ic=rep(ic0[1:NV],2); } else if(is.null(ic)) ic=ic0;
	
	if(ntests<5) Pmat[,"s"]=apply(Pmat,1,attnOptim) else {
		nt2=5; PmatBrief=Pmat[1:(nt2^2),]; PmatBrief[,varnms]=as.matrix(expand.grid(seq(varLims[1,1],varLims[2,1],len=nt2),seq(varLims[1,2],varLims[2,2],len=nt2))); 
		if(max(Pmat[,"p"])>0){ PmatBrief[,"s"]=apply(PmatBrief,1,attnOptim); sspln=splinefun2(PmatBrief[,c("s",varnms)],k=15); Pmat[,"s"]=sspln(Pmat[,varnms[1]],Pmat[,varnms[2]]); }; }
	
	res0=lsoda(c(rep(head(ic,NV),N),rep(tail(ic,NV),N)),1:Trun,wind4,Pmat)[,-1]; res=round(apply(res0[,imitig], 2, function(x) tail(na.omit(x),1)),3); 
	if(give==-1) return(res0); if(give==-2) return(list(res0[,1:6],Pmat[1,]));
	HIC=matrix(head(res,N),ntests); LIC=matrix(tail(res,N),ntests); dS=HIC-LIC; rownames(dS)=vars[,1]; colnames(dS)=vars[,2];
	dS1=dS2=dS4=dS; dS1[1:N]=head(res,N); dS2[1:N]=tail(res,N); if(give==1) return(dS1); if(give==2) return(dS2); if(give==3) return(abind(dS1,dS2,along=3));
	outP=2+0*dS; outP[dS>0.1]=0; outP[LIC>0.9*max(HIC)]=1; out=pmax(abind(outP,HIC,LIC,dS,along=3),0)
	SpatiotempPlot(out[,,plotgive],Yax=vars[,1],Xax=vars[,2],XaxN=XN,YaxN=YN,figtitle=paste0(c(XN," (X) vs. ",YN," (Y)"),collapse=""),cexAx=1.3); out;
}
#if bt is being varied, assdoTclose assumes it is the second column of varLims and second element of ntests
assdoTclose=function(parms,varLims=cbind(k=c(0.1,0.5),bt=c(0,10)),Trun=1e2,ntests=40,ic=NULL,closeProp=0.05,plotgive=1,Conts=c(20,40,80,100),Zlim=c(1.8,4.6),YN="Proportion of mitigation dependent on buy-in",XN="Time delay in starting mitigation",Fill=FALSE,sfix=NA){
	ntests=c(head(ntests,1),tail(ntests,1))
	varnms=colnames(varLims); vars=apply(rbind(varLims,ntests),2,function(x) round(seq(x[1],x[2],len=x[3]),4));
	if(is.matrix(vars)) vars=split(t(vars),1:ncol(vars)); N=prod(unlist(lapply(vars,length)));
	Pmat=do.call("rbind",rep(list(parms),N)); Pmat[,varnms]=as.matrix(expand.grid(vars[[1]],vars[[2]]))[,1:ncol(varLims)];
	ic0=c(1,0.05,0.75,0.025,1,0.01,  0.75,0.05,0.75,0.75,0.3,0.01); if(is.null(ic)) ic=ic0; NV=length(ic0)/2; imitig=c(TRUE,rep(FALSE,NV-1));
	
	#Here, calibrate pop size s so resource reliance and public attention have equal influence
	#calibration takes time and s changes smoothly, so if there's lots of trials, calculate s at a few parm sets and interpolate using splines
	if(N<50) Pmat[,"s"]=apply(Pmat,1,attnOptim) else {
		nt2=5; PmatBrief=Pmat[1:(nt2^2),]; PmatBrief[,varnms]=as.matrix(expand.grid(seq(varLims[1,1],varLims[2,1],len=nt2),seq(varLims[1,2],varLims[2,2],len=nt2))); 
		if(max(Pmat[,"p"])>0){ PmatBrief[,"s"]=apply(PmatBrief,1,attnOptim); sspln=splinefun2(PmatBrief[,c("s",varnms)],k=15); Pmat[,"s"]=sspln(Pmat[,varnms[1]],Pmat[,varnms[2]]); }; }
	if(!is.na(sfix)) Pmat[,"s"]=sfix
	
	res0=lsoda(rep(ic[1:NV],N),1:Trun,wind4,Pmat)[,-1][,imitig];
	if(plotgive=="bfn") return(matrix(res0[Trun,],ntests[1],ntests[2]))
	if(!"bt"%in%varnms){
		Conts=c(10,20,40,80); res0=tail(res0,-parms["bt"]); parms["bt"]=0; base0=lsoda(rep(ic[NV+1:NV],N),1:Trun,wind4,Pmat)[,-1][,imitig]; 
	} else base0=lsoda(rep(ic[NV+1:NV],ntests[1]),1:Trun,wind4,Pmat[1:ntests[1],])[,-1][,imitig]
	base=round(apply(base0, 2, function(x) tail(na.omit(x),1)),3);
	resTclose=apply((res0[-1,]%*%diag(rep(1/base,N/length(base))))>(1-closeProp), 2, function(x) DelDegr+c(which(x[-(1:DelDegr)])[1],NA)[1])
	if(!"bt"%in%varnms){ resTclose[base<mean(base)]=0; resTclose[resTclose>0.66*Trun]=NA; }
	out0=out=matrix(log(1+resTclose),ntests[1],ntests[2]); out[is.na(out)]=-1; out=pmin(out,quantile(out,0.95)); 
	if("bt"%in%varnms) out[out==min(out,na.rm=TRUE)]=-1 else { out0[out0==0]=NA; out[out==0]=log(2); }
	out[apply(out,1,sd)==0,]=-1; if(plotgive==0) return(out);
	if(max(out)>Zlim[2]) Zlim[2]=1.01*max(out); if(Fill) out[out>0 & out<Zlim[1]]=Zlim[1];
	SpatiotempPlot(out,Yax=vars[[1]],Xax=vars[[2]]-DelDegr*(varnms[2]=="bt"),XaxN=XN,YaxN=YN,cexAx=1.3,contPlot=exp(out0)-1,cont=1,Conts=Conts,contSpan=1,cexCont=1.5,Zlim=Zlim); out;
}
attnOptim=function(parms,IC=c(N=1,A=0.05,U=.75,m=0.025,S=1,I=0.01),fitWindow=20+DelDegr+RampDegr){
	Loglikfun=function(x){ parms["s"]=x; TSi=lsoda(IC,0:fitWindow,wind4,parms)[,1+c(6,3)]; if(DelDegr>0) TSi=tail(TSi,-DelDegr); diff(c(x,1)*colMeans(TSi))^2; }
	fit=nloptr(1,Loglikfun,lb=0.1,ub=50,opts=list("algorithm"="NLOPT_GN_DIRECT", "print_level"=0,  "xtol_rel"=0.1, "maxeval"=10)); fit$solution;
}


#ICs are c(R,A,U,m,S,I)
Dtype=2; DR=0; DelDegr=0; RampDegr=0;
WindVis=function(parms,BT=c(0,15),Trun=1e2,cols=c("black","indianred","blue","magenta","orange"),colsP=c("darkgreen","red"),RW=0,...){
	tseries=assdo2(parms,varLims=cbind(bt=BT,omg=c(1,1)),ic=rep(c(1,0.05,0.95,0.025,.99,0.01),2),give=-1,Trun=Trun,ntests=2)[,c(1:4,6:10,12)]; tseries[,c(5,10)]=tseries[,c(5,10)]/max(tseries[,c(5,10)]);
	cols[c(2*(parms["piA"]==0),5*(parms["p"]==0))]="white"; tseries[,cols=="white"]=NA;
	matplot((1:nrow(tseries))-DelDegr,tseries,col=rep(cols,2),lty=rep(c(1,3),each=5),lwd=c(3+RW,rep(3,4)),type="l",ylim=c(0,1.1))
}
wind4=function(t,start,Parms){ #for parallel replicates, start=rep(c(R,A,U,m),ntrials)
	IC=matrix(pmax(start,1e-5),ncol=6,byrow=TRUE); R=IC[,1]; A=IC[,2]; U=IC[,3]; m=IC[,4]; S=IC[,5]; I=IC[,6];
	if(nrow(IC)==1) Parms=t(Parms); f=Parms[,"f"]; 
	if("FDYN"%in%colnames(Parms)) FDYN=Parms[1,"FDYN"];
	if(!is.na(FDYN)) f=fdyn(t,fb=f,tmax=FDYN);
	Parms=cbind(Parms,b2=1); Parms[,c("b","k","b2")[Dtype]]=Parms[,c("b","k","b2")[Dtype]]*(t>=Parms[,"bt"]);
	dR = Parms[,"rR"]*R*(1-R) - (t>DelDegr)*Parms[,"mu"]*R*(1-Parms[,"b"]*(1-f+f*Parms[,"b2"]*m)) #+ 0.01
	dA = Parms[,"rA"]*A*(1-A) - Parms[,"alpha"]*R*A
	dU = Parms[,"xi"]*(R - U)
	dm = m*Parms[,"k"]*(1-m)*(Parms[,"piR"]*(U+Parms[,"p"]*(I*Parms[,"s"] - U)) - Parms[,"piA"]*A - 1)
	
	learn=Parms[,"beta"]*S*I*(pmax(omega-R,0))
	dS = - learn; 	dI = learn - I*Parms[,"delta"];
	return(list(as.numeric(t(cbind(dR,dA,dU,dm,dS,dI)))))
}
wind4null=function(t,start,Parms){ #for parallel replicates, start=rep(c(R,A,U,m),ntrials)
	IC=matrix(pmax(start,1e-5),ncol=6,byrow=TRUE); R=IC[,1]; A=IC[,2]; U=IC[,3]; m=IC[,4]; S=IC[,5]; I=IC[,6];
	if(nrow(IC)==1) Parms=t(Parms); f=Parms[,"f"]; 
	Parms=cbind(Parms,b2=1); Parms[,c("b","k","b2")[Dtype]]=Parms[,c("b","k","b2")[Dtype]]*(t>=Parms[,"bt"]);
	dR = Parms[,"rR"]*R*(1-R) - (t>DelDegr)*Parms[,"mu"]*R*(1-Parms[,"b"]*(1-f+f*Parms[,"b2"]*m)) #+ 0.01
	dA = Parms[,"rA"]*A*(1-A) - Parms[,"alpha"]*R*A
	dU = Parms[,"xi"]*(R - U)
	dm = m*Parms[,"k"]*(1-m)
	
	learn=Parms[,"beta"]*S*I*(pmax(omega-R,0))
	dS = - learn; 	dI = learn - I*Parms[,"delta"];
	return(list(as.numeric(t(cbind(dR,dA,dU,dm,dS,dI)))))
}


#Note some differences in variable names between paper and code:
# Natural-state service (N) = R
# Reliance (R) = U
# Unaware (U) = S
# Aware (D) = I

#And in parameters:
# kappa = k
# r_N = rR
# pi_N = piR
# f in paper (=0.15) is 1-f here (=0.85)




#Set parameters for each conservation scenario
FDYN=NA; ZLM=c(2.7,4.52); COL=c("black","red3","blue","magenta","orange");
omega=0.77; parmsS1=c(rR=.42,rA=.47,mu=.31,b=0.73,bt=0,f=0.85,alpha=1,  k=0.3,piR=3.5,piA=0,xi=0.22,omg=1,  p=0.3,beta=9,delta=1.4,s=1)
parmsS2=parmsS1; parmsS2[c("p","mu")]=c(0.7,0.22);
parmsS3=parmsS1; parmsS3[c("f","piA","p","mu")]=c(0.7,4,0,.375); c(Fmsy=parmsS3["rR"]/4,F.Fmsy=parmsS3["mu"]/(parmsS3["rR"]/4))



################ FIGURE 2
DelDegr=5; par(mfrow=c(1,3))
WindVis(parmsS1,BT=5+c(0),Trun=55,cols=COL,las=1,RW=3);
WindVis(parmsS1,BT=5+c(12),Trun=55,cols=COL,las=1,RW=3);
WindVis(parmsS1,BT=5+c(18),Trun=55,cols=COL,las=1,RW=3);




################ FIGURE 3
#Note: approximate runtime is 10-15 minutes. Set NTR=25 to speed up.
par(mfrow=c(1,3))
btm=60; piRlm=c(3,7.5); NTR=50; DelDegr=0;
adapt=assdoTclose(Trun=5e2,varLims=cbind(piR=piRlm,bt=c(0,btm+DelDegr)),ntests=NTR,Zlim=ZLM,plotgive=2,parms=parmsS1,YN="Natural-state resource payoff",XN="Time delay in starting mitigation",Fill=TRUE)
opinion=assdoTclose(Trun=5e2,varLims=cbind(piR=piRlm,bt=c(0,btm+DelDegr)),ntests=NTR,Zlim=ZLM,plotgive=2,parms=parmsS2,YN="Natural-state resource payoff",XN="Time delay in starting mitigation")
#Note: "pit" in restoration times @pi>6 is because it "helps" to wait till things get bad enough for outcry
altres=assdoTclose(Trun=2e2,varLims=cbind(piR=piRlm,bt=c(0,btm+DelDegr)),ntests=NTR,Zlim=ZLM,plotgive=2,parms=parmsS3,YN="Natural-state resource payoff",XN="Time delay in starting mitigation",Fill=TRUE)





# FIGURE 3a: Simplify heatmaps into plots
trlplot=function(sim,piShow=c(4.25,7)){
	idx=sapply(piShow,function(x) which.min((seq(piRlm[1],piRlm[2],len=nrow(sim))-x)^2))
	oot=t(sim[idx,,drop=FALSE]); oot[oot==-1]=NA; oot[oot==max(sim)]=NA
	plt=cbind(seq(0,btm,len=ncol(sim)), oot)
	matplot(plt[,1],plt[,-1],type="l",lty=1,lwd=3,col=c(1,4),las=1,ylim=ZLMi,xlim=c(0,xmx))
}
par(mfrow=c(1,3)); DF=4.35; ZLMi=c(ZLM[1],5); xmx=30; trlplot(adapt,DF); trlplot(opinion,DF); trlplot(altres,DF)



################ FIGURE 3b
#Window lengths for each level of resource value
plotdat=function(d,rumn=4,nas=FALSE){
	BTseq=seq(0-DelDegr,delMax,len=ncol(d))
	WLs=apply(d,1,function(x) BTseq[which(x==-1)[1]])
	if(nas) WLs[is.na(WLs)]=999; WLs[WLs==DelDegr]=NA;
	if(rumn>0) WLs=as.numeric(filter(WLs,rep(1/rumn,rumn)))
	cbind(seq(piRlm[1],piRlm[2],len=nrow(d)),WLs)
}
#heatmap-based visual:
delMax=25; RM=3; par(mfrow=c(3,1))
xd=plotdat(adapt,0); SpatiotempPlot(t(t(xd[,2])),Yax=xd[,1],Zlim=c(0,delMax),pal=2)
xd=plotdat(opinion,RM); SpatiotempPlot(t(t(xd[,2])),Yax=xd[,1],Zlim=c(0,delMax),pal=2)
xd=plotdat(altres,RM); SpatiotempPlot(t(t(xd[,2])),Yax=xd[,1],Zlim=c(0,delMax),pal=2)
par(mfrow=c(1,1))




########### FIGURE S1
bfnplot=function(ResLvls,X=seq(0,8,len=nrow(ResLvls))) matplot(X,ResLvls,type="l",lwd=3,col=c(1,2),ylab="",las=1)#Steady-state natural ecosystem service level",xlab="Value of natural ecosystem service")
DelDegr=150; FDYN=NA; COL=c("black","red3","blue","magenta","orange");
ass1=assdoTclose(Trun=1e3,varLims=cbind(piR=c(0,8),bt=c(0,200+DelDegr)),ntests=c(30,2),plotgive="bfn",parms=parmsS1)
ass2=assdoTclose(Trun=1e3,varLims=cbind(piR=c(0,8),bt=c(0,200+DelDegr)),ntests=c(30,2),plotgive="bfn",parms=parmsS2)
ass3=assdoTclose(Trun=3e2,varLims=cbind(piR=c(0,8),bt=c(0,200+DelDegr)),ntests=c(30,2),plotgive="bfn",parms=parmsS3)
par(mfrow=c(1,3)); bfnplot(ass1); bfnplot(ass2); bfnplot(ass3); par(mfrow=c(1,1)); 



########### FIGURE S4
wind4=wind4null; NTR=25; DelDegr=5;
adaptNULL=assdoTclose(Trun=1e2,varLims=cbind(piR=piRlm,bt=c(0,btm+DelDegr)),ntests=NTR,Zlim=ZLM,plotgive=2,parms=parmsS1,Fill=TRUE)
opinionNULL=assdoTclose(Trun=1e2,varLims=cbind(piR=piRlm,bt=c(0,btm+DelDegr)),ntests=NTR,Zlim=ZLM,plotgive=2,parms=parmsS2)
altresNULL=assdoTclose(Trun=1e2,varLims=cbind(piR=piRlm,bt=c(0,40+DelDegr)),ntests=NTR,Zlim=ZLM,plotgive=2,parms=parmsS3,Fill=TRUE)

par(mfrow=c(1,3)); DF=4.33; ZLMi=c(ZLM[1],5); xmx=30; trlplot(adaptNULL,DF); trlplot(opinionNULL,DF); trlplot(altresNULL,DF)







################ FIGURE S3a: role of concentrating government effort up-front
# given a fixed investment and time range, calculate window size along a constant-cost (isocost) curve.
FDYN=NA; fdyn=function(t,fb,tmax) cbind(1,fb)[,1+(t<=tmax)]
isocostWinds=function(parmsi,tendLim=75*c(0.8,1.2),costLim=c(10,30),ntests=c(tend=15,fs=25,bt=15),BTmx=50){
	tendSeq=round(seq(tendLim[1],tendLim[2],len=ntests[1]),1); out=matrix(0,ntests[1],ntests[2]); 
	colnames(out)=round(seq(costLim[2],costLim[1],len=ntests[2]),1); rownames(out)=tendSeq;
	winds=out; for(i in 1:ntests[1]){ 
		assign("FDYN",tendSeq[i],envir=.GlobalEnv); fis=1-costLim/FDYN; 
		if(min(fis)<0) print("NOTE: govt effort 1-f exceeds 1");
		tri=assdoTclose(Trun=4e2,varLims=cbind(f=c(pmax(fis[2],0.01),fis[1]),bt=c(0,BTmx+DelDegr)),ntests=ntests[2:3],Zlim=ZLM,plotgive=0,parms=parmsi)
		winds[i,]=apply(tri,1,function(y) seq(0,BTmx+DelDegr,len=ntests[3])[which(y==-1)[1]])
	}; assign("FDYN",NA,envir=.GlobalEnv); winds;
}

#Hybrid case. Note: approximate runtime is 5 minutes
Mag=0.3; isocostFig3=isocostWinds(c(piR=4.5,(parmsS1+parmsS2)/2),tendLim=c(20,60),costLim=60*0.15*c(1-Mag,1+Mag),ntests=c(24,3,60),BTmx=50)
X=as.numeric(rownames(isocostFig3)); Y=rev(as.numeric(colnames(isocostFig3)));
EffortWasteCheck=matrix(apply(expand.grid(X,Y),1,function(x) x[2]>x[1]), length(X),length(Y))
cand=t(apply(isocostFig3,1,rev)); cand[EffortWasteCheck]=NA; #matplot(X,cand,type="o",lwd=2,main=Y);
plot(supsmu(X,cand[,2]),type="l",lwd=3,main=Y,las=1,ylim=c(25,40))






############################# FIGURE 1a - Public opinion dynamics 
# This analysis is built on scraping data from newspapers.com. See APPENDIX at end of script for functions that did this.
library(rio)
ia=import("IssueAttnCycles_Fig1.csv")
X=-5:15; iap=apply(matrix(ia$Hits[ia$YrSinceMax%in%X],nrow=length(X)),2,function(x) (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
gampred(cbind(rep(X,ncol(iap)),as.vector(iap)),1,2,spanSignif=NA,col="lightgray",seFill=1)
matlines(X,iap,type="o",pch=16,lwd=1,lty=1,col=1:7); 
legend("topright",text.col=8,c("Chesapeake pfiesteria, MD","DDT, USA","Chestnut Blight, PA","Lake pollution, WA","Zebra mussels, OH","Erie pollution, OH","Erie algae, OH"),seg.len=1.5,lwd=4,col=1:7,box.col=0)



############################# Environmental trends and parameterization
etrends=import("/ErieEnvTrends.csv"); colnames(etrends)
########## FIGURE 4a
Ar=etrends$maxBloom_km2; pctAr=c(30,NA,NA,30, 100*(1-Ar[-(1:4)]/2867)); X=etrends$Year; 
gampred(cbind(X,pctAr),1,2,spanSignif=NA,ylim=c(20,100),cex=1.2,cex.axis=1.15,las=1,seFill=0.07); points(X,pctAr,pch=16,col=8);

########## FIGURE 4b
gampred(cbind(100*etrends$maxBloom_km2/2867,etrends$medNdays),1,2,k=3,col=1,pch=16,seFill=0.1,spanSignif=NA,ylim=c(5,50))
par(new=TRUE)
gampred(cbind(100*etrends$maxBloom_km2/2867,etrends$medSevere),1,2,k=3,pch=16,col="indianred",seFill=0.1,spanSignif=NA,axes=FALSE,ylim=c(0,7))
axis(4,col="indianred",col.axis="indianred")


RunAvg=function(x,nvals=5,interpol=TRUE,past=FALSE,vec=is.null(dim(x)),hlf=(nvals-1)/2){
	if(vec) x=cbind(1,1,x); if(interpol) x[,3]=approx(x[,3],xout=1:nrow(x))$y;
	for(i in unique(x[,1])){ use=x[,1]==i;	if(sum(use)>nvals) x[use,]=cbind(x[use,1:2],c(head(x[use,3],hlf),tail(head(filter(as.numeric(x[use,3]),rep(1/nvals,nvals),sides=2-past),-hlf),-hlf),tail(x[use,3],hlf))); }; 
	if(vec) return(x[,3]); x;
}
cumPlot=function(dat,colors,x=1:nrow(dat),smoothe=15,nonsmoo=NULL,center=FALSE,...){ #smoothe<5 drives running averages of window |smoothe|
	out=dat; rownames(out)=1:nrow(out); N=ncol(dat); if(smoothe<(-1)) out=apply(dat,2,RunAvg,nvals=abs(smoothe)); 
	if(smoothe>0) for(i in 1:N){ pred=predict(gam(dat[,i]~s(x,k=smoothe))); out[rownames(out)%in%names(pred),i]=pred; }
	if(!is.null(nonsmoo)) out[nonsmoo,]=as.matrix(dat[nonsmoo,])
	out[is.na(out)]=0; out2=cbind(0,t(apply(out,1,cumsum))); out2=out2 - center*matrix(rowSums(out2)/(N+1),length(x),(N+1));
	matplot(x,out2,col=0,...)
	for(i in 1:N) polygon(c(x,rev(x)),c(out2[,i+1],rev(out2[,i])),col=colors[i])
}
#Smooth DRP with extrapolation to early years
etrends$DRPs=RunAvg(etrends$AnnualDRPload_t)
etrends$TTP=rowSums(etrends[,c("NonpointTP","UnmonitoredAdjTP","PointTP","AtmosphericTP","HuronTP","UndeterminedTP")],na.rm=TRUE)
etrends$TTPs=RunAvg(etrends$TTP)
drpms=gam(AnnualDRPload_t~s(TTP)+s(Year),data=etrends); #summary(drpms)
drpsp=predict(drpms,newdata=data.frame(Year=etrends$Year,TTP=etrends$TTPs),link="response",se.fit=TRUE)
dats=head(tail(cbind(etrends$Year,drpsp$fit,1.96*drpsp$se.fit),-2),-2)

########## FIGURE 5g
AL=1; COLS=c(rgb(.75,.75,.75,alpha=AL),rgb(0,.55,.27,alpha=AL),rgb(.8,0,0,alpha=AL))
cumPlot(etrends[,c("UndeterminedTP","NonpointTP","PointTP")],COLS,smoothe=-5,nonsmoo=1:5,x=etrends$Year,ylab="")
par(new=TRUE); plot(dats[,1:2],axes=FALSE,xlab="",ylab="",type="l",lwd=4,col="orange",ylim=c(100,900)); axis(side=2,line=3,col="orange",col.axis=6); 
polygon(c(dats[,1],rev(dats[,1])),c(dats[,2]+dats[,3],rev(dats[,2]-dats[,3])),border=0,col=rgb(1,.65,0,alpha=1/3))











############################# IV - Model parameterization
#finding r and mu values
nrm=function(x,Min=TRUE){
	minScl=0; if(Min) minScl=min(x,na.rm=TRUE);
	return((x-minScl)/(max(x,na.rm=TRUE)-minScl))
}
datFit=cbind(Y=X,N=pctAr/1e2,E=nrm(etrends$AnnualDRPload_t,FALSE),Ns=predict.gam(gam((pctAr/100)~s(X)), newdata=data.frame(X=X)))
datFit=cbind(datFit,Es=RunAvg(datFit[,"E"]))
ode=function(t,start,parmsFit,Var){
	N=pmax(start,1e-6); list(N*parmsFit["r"]*(1-N) - N*(parmsFit["mu0"]*parmsFit["r"]+parmsFit["gP"]*Var))
}
LLfun=function(parmsFiti,stats=FALSE){
	parmsFit[pfit]=parmsFiti
	Np=dati[1,"Ntm1"]; for(i in 1:nrow(dati)) Np=c(Np,lsoda(tail(Np,1),0:1,ode,parmsFit,Var=dati[i,"E"])[2,2]); Np=Np[-1];
	if(stats) return(cbind(dati,Np)); sum((Np-dati[,"N"])^2);
}
V="Ns"; lm(datFit[,V]~datFit[,"Es"]); plot(datFit[,c("Es",V)],xlim=0:1)
# Note the lm intercept: 1-mu0/r=Nbar=0.969. 
# So  mu0=0.03*r - baseline proportion of eastern Lake Erie covered by blooms without nutrient pollution


#Now fit model to get natural-state service growth rate, r_N
NT=350; parmsFit=c(r=1,mu0=0.03,gP=0.3); Ncol="Ns"; Ecol="E"; pfit=c(1,3); library(deSolve); library(nloptr);
UB=c(2,1,1); LB=c(0.1,0.01,0); parmnames=names(LB)=names(UB)=names(parmsFit);
dati=na.omit(cbind(Y=datFit[,"Y"],N=datFit[,Ncol],E=datFit[,Ecol],Ntm1=head(c(NA,datFit[,Ncol]),-1)))
fit=nloptr(opts=list("algorithm"="NLOPT_GN_DIRECT","xtol_rel"=1e-5,"maxeval"=NT,"print_level"=0), x0=parmsFit[pfit],eval_f=LLfun,lb=LB[pfit],ub=UB[pfit],stats=FALSE)

#Check model fit
s=fit$solution; x=LLfun(s,TRUE); c(rN=s[1],mu=s[1]*parmsFit["mu0"]+s[2],SSE=LLfun(s)); round(cor(x,use="pairwise.complete.obs")^2,2)["Np",2:4]
g_realscale=s[2]/max(RunAvg(etrends$AnnualDRPload_t),na.rm=TRUE) #correct for different units of runAvg SRP (normalized relative to max)
plot(datFit[,c("Y","N")],pch=16); lines(x[,1],x[,"Np"],lwd=4,col=4); lines(datFit[,c("Y","Ns")],col=8,lwd=3);

#Now back-calculate pollution rates
mu_1970s=mean(drpsp$fit[etrends$Year<=1970]*g_realscale+parmsFit["mu0"]); mu_1970s
mu_2010s=mean(drpsp$fit[etrends$Year%in%2005:2020]*g_realscale+parmsFit["mu0"]); mu_2010s


#Calculating upper bound on social learning rate - specifically, kappa*Pi_t
Pt=etrends$PointTP[etrends$Year%in%(1970:1982)]
mt=1-(Pt-tail(Pt,1))/(Pt[1]-tail(Pt,1)); mtm1=c(NA,head(mt,-1));
rpt=log(mt/mtm1)/(1-mtm1); rpt=rpt[rpt<Inf & !is.na(rpt)][-1];
plot(rpt[rpt<Inf]); 
mean(rpt); sd(rpt)/sqrt(length(rpt)) #Mean kappa times Pi_t +/- SE
#As Pi_t was likely >1, this is an upper-bound estimate














############## APPENDIX
#####Obtaining newspapers.com data for Figure 1a
#####This code produced IssueAttnCycles_Fig1.csv
#This code is meant for documentation: it no longer works as of 2023 due to website requiring a paid account for any search
#There wasn't a way to get the total number of articles, so we just took the number of articles mentioning any of the baseline terms
#Carefully selected baselines to avoid any plausibly marketing related terms, as marketing heavily shifted onto internet around 2005
# popular=function(states,term,years,baselines=c("the*","high","sky")){
# 	ns=length(states); res=matrix(0,length(years),2+2*ns); res[,1]=years;
# 	colnames(res)=c("Yr","CumFreq",paste0(rep(states,2),rep(c(".hits",".numpubs"),each=ns)))
# 	numfind=function(needle,haystack) c(as.numeric(gsub("[^\\d]+", "", haystack[grep(needle,haystack)+2], perl=TRUE)),0)[1];
# 	for(t in years){
# 		hits=readLines(paste0("https://go.newspapers.com/results.php?query=",term,"&dr_year=",t));
# 		base1=base2=base3=0; if(!is.na(baselines)){
# 			base1=readLines(paste0("https://go.newspapers.com/results.php?query=",baselines[1],"&dr_year=",t));
# 			base2=readLines(paste0("https://go.newspapers.com/results.php?query=",baselines[2],"&dr_year=",t));
# 			base3=readLines(paste0("https://go.newspapers.com/results.php?query=",baselines[3],"&dr_year=",t)); }
# 		for(j in states) res[t==years,2+which(states==j)+c(0,ns)] = c(numfind(j,hits), (numfind(j,base1)+numfind(j,base2)+numfind(j,base3))/0.07);
# 	}
# 	res[,2]=rowSums(t(t(res[,2+(1:ns)])))/rowSums(t(t(res[,-(1:(2+ns))])));
# 	if(!is.na(baselines)) plot(res[,1:2],type="o",ylab="Proportion of newspaper articles",main=term); res;
# }
#substitute + for spaces in multiple-word searches. %2a is wildcard, eg "alga%2a+bloom".
# cb=popular("Pennsylvania","chestnut+blight",1908:1925)
# zm=popular("Ohio","zebra+mussel",1984:2020); 
# polWA=popular("Washington","lake+pollution",1964:1980)
# polCB=popular("Maryland",c("Chesapeake+bacteria"),1991:2020)
# aboh1=popular("Ohio","pollution+Erie",1960:1980);
# aboh2=popular("Ohio","pollution+Erie",2005:2021);
# dt=popular("USA","DDT",1964:1990);




