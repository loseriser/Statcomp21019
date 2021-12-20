#' @title A dose regimen dataset
#' @name Sa
#' @description A dataset used to show differnet  dose regimen.
#' @examples
#' \dontrun{
#' data(Sa)
#' attach(Sa)
#' t<-seq(from=0,to=35,by=0.1)
#' plot(t,Concentration(t=t,d=rep(25,7)),type = "l",lty=2,ylab="Concentration ")
#' lines(t,Concentration(t=t,d=dose2),type = "l")
#' }
NULL

#' @title drug concentration 
#' @description A 1-compartment infusion model for the PK
#' @param t the time 
#' @param d the dose of drug
#' @param t0 Administration time
#' @param Cl Clearance of elimination
#' @param V Volume of distribution
#' @param tinf the duration of the infusion of the jth administration
#' @return  the cytokine response
#' @examples
#' \dontrun{
#' t<-seq(from=0,to=35,by=0.1)
#' plot(t,Concentration(t=t,d=rep(25,7)),type = "l",lty=2,ylab="Concentration ")
#' lines(t,Concentration(t=t,d=c(1,5,10,rep(25,4))),type = "l")
#' }
#' @export
Concentration<-function(t,d,t0=c(1, 5, 9, 13, 17, 21, 25),Cl=1.36*10,V=3.4,tinf=rep(1e-5,7)){
  k=Cl/V
  stopifnot((length(d)==length(t0))&&(length(tinf)==length(t0)))
  n<-length(d)
  cal<-function(t,t0,tinf,d,k,V){
    cc<-d/(tinf*k*V)
    ifelse(t-t0<=0,0,ifelse(t-t0>tinf,cc*(1-exp(-k*tinf))*exp(-k*(t-t0-tinf)),cc*(1-exp(-k*(t-t0)))))
  }
  SCC<-sapply(t,function(t){
    sum(50*unlist(sapply(1:n,function(i){
      cal(t,t0=t0[i],tinf=tinf[i],d=d[i],k=k,V=V)
    })))
  }
  )
}

#' @title cytokine Response 
#' @description A cytokine response function for the PD model
#' @param t the time 
#' @param d the dose of drug
#' @param t0 Administration time
#' @param EMAX Maximum cytokine release rate
#' @param EC50 Drug exposure for half-maximum release
#' @param H Hill coefficient for cytokine release
#' @param IMAX Maximum inhibition of cytokine release
#' @param IC Cytokine exposure for half-maximum inhibition
#' @param KDEG Degradation rate for cytokine
#' @param K Priming factor for cytokine release
#' @return  A list of the cytokine response , its derivatives and its AUC.
#' @examples
#' \dontrun{
#' t<-seq(from=0,to=35,by=0.02)
#' plot(t,Cytokine(t=t,d=rep(25,7))$EE,type="l",lty=2,ylab="E(t)")
#' abline(h=max(Cytokine(t=t,d=rep(25,7))$EE),lty=2)
#' lines(t,Cytokine(t=t,d=c(1,5,10,rep(25,4)))$EE,type="l")
#' abline(h=max(Cytokine(t=t,d=c(1,5,10,rep(25,4)))$EE),lty=1)
#' }
#' @export
Cytokine<-function(t,d,t0=c(1, 5, 9, 13, 17, 21, 25),EMAX=3.59e5,EC50=1e4,H=0.92*0.9,
             IMAX=0.995,IC=1.82e4,KDEG=3.6,K=2.83){
  C_1<-Concentration(t=t,d=d)
  J<-length(t0)
  n<-length(t)
  dE<-EE<-AUCE<-aa<-numeric(n)
  dd<-diff(t)
  for(i in 2:n){
    dE[i-1]<-EMAX*C_1[i-1]^H/(EC50^H+C_1[i-1]^H)*(1-((IMAX*AUCE[i-1])/(IC/K^(J-1)+AUCE[i-1])))-KDEG*EE[i-1]
    EE[i]<-EE[i-1]+(dE[i-1])*dd[i-1]
    AUCE[i]<-AUCE[i-1]+(EE[i-1])*dd[i-1]
  }
  return(list(dE=dE,EE=EE,AUCE=AUCE))
}
