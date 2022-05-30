objective_function<-function(ob,error,bsz,subind,lmd,R,t){
### This prepares the objective function

    tsz=length(t)

    bsz=bsz
    ob=ob
    sz=length(ob)

    subind=subind

### This part depends on how the error is specified.
    if(is.list(error)){
        family=error[[1]]
        if(family=="Normal"){
            mu=error[[2]]
            sigma=error[[3]]
            integrate_density<-function(t,infn0){
                return((-infn0-mu)*stats::pnorm(-infn0,mu,sigma)+sigma^2*stats::dnorm(-infn0,mu,sigma))
            }
        }else if(family=="Laplace"){
            a=error[[2]]
            b=error[[3]]
            integrate_density<-function(t,infn0){
                infn1=a+infn0+b*rmutil::plaplace(-infn0,a,b)
                infn1[-infn0>a]=(b-b*rmutil::plaplace(-infn0,a,b))[-infn0>a]
                return(infn1)
            }
        }else if(family=="Beta"){
            a=error[[2]]
            b=error[[3]]
            integrate_density<-function(t,infn0){
                return(-infn0*stats::pbeta(-infn0,a,b)-a/(a+b)*stats::pbeta(-infn0,a+1,b))
            }
        }else{
            stop(error=paste("Error family \"",family,"\" not yet implemented",sep=""))
        }
    }else{
        esz=length(error)
        integrate_density<-function(t,infn0){
            return(apply(infn0,c(1,2),function(x){sum((x+error)[x+error>0])/(esz)}))
        }
    }

    subind=subind



    getMatrices<-function(t){
        dett=t[2]-t[1]
### integrate density

        infn0=-ob%*%t(rep(1,tsz))+rep(1,sz)%*%t(t)  # n*T
        infn1=integrate_density(t,infn0)            # n*T

### support basis

        nnbasis=t(abs(t-rep(1,tsz)%*%t(t[seq(5,tsz-5,length.out=30)])))/2 # T*30
### These appear to be needlessly limited to values from t.
### They are not evenly spaced. For tsz<30, there could be repeats.
### In practice tsz=1000, so it should be OK.
### Could make the number of support points a parameter

        basr=rbind(infn1[subind,],t,rep(1,tsz),t^2,nnbasis)

        Ur.vals=rbind(t,rep(1,tsz),t^2,t(abs(t-rep(1,tsz)%*%t(t[seq(5,tsz-5,length.out=30)])))/2)
        Ur=Ur.vals%*%t(basr)*dett
        Ur3=Ur[,-((bsz+1):(bsz+3))]

###        print(basr)

### get initial values
        prd=t(abs(t%*%t(rep(1,tsz))-rep(1,tsz)%*%t(t))/2)%*%t(basr)*dett

        y.inio=rep(1,tsz)
        ydst=stats::density(ob,n=1e5)
        for(i in 1:tsz){ ## evaluate density at the nearest calculated point to each t[i]
            y.inio[i]=ydst$y[which.min(abs(ydst$x-t[i]))]
        }

        fit.inio=stats::lm(y.inio~prd-1)
        inio0=fit.inio$coefficients
        inio0[is.na(inio0)]=1e-18

        return(list("infn1"=infn1,"basr"=basr,"Ur"=Ur,"Ur3"=Ur3,"inio0"=inio0))
    }

### calculate lambda

    if(is.null(lmd)){
        matrices<-getMatrices(t)
        dett<-t[2]-t[1]


        x<-matrices$inio0
        Ur<-matrices$Ur
        Ur3<-matrices$Ur3
        basr<-matrices$basr
        infn1<-matrices$infn1

        x[bsz+1]=(t((Ur[2,bsz+3]*Ur3[1,]-Ur[1,bsz+3]*Ur3[2,])*(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2])- (Ur3[2,]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur3[3,])*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2]))%*%x[-((bsz+1):(bsz+3))]-2*Ur[2,bsz+3]*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2]))/ ((Ur[2,bsz+1]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+1])*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2])- (Ur[1,bsz+1]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+1])*(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2]))
        x[bsz+2]=((Ur[2,bsz+3]*Ur3[3,]-Ur3[2,]*Ur[3,bsz+3])%*%x[-((bsz+1):(bsz+3))]-(Ur[2,bsz+1]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+1])*x[bsz+1]-2*Ur[2,bsz+3])/(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2])
        x[bsz+3]=(-Ur3[1,]%*%x[-((bsz+1):(bsz+3))]-Ur[1,bsz+1]*x[bsz+1]-Ur[1,bsz+2]*x[bsz+2])/Ur[1,bsz+3]

        lmd=sum(abs(-rowSums((basr%*%t(infn1))%*%diag(1/(as.vector(((t(x)%*%basr)%*%t(infn1))))))))/sum(abs(2*(basr%*%t(t(x)%*%basr)*dett)))/R

    }

    getfunction<-function(t){
### Calculate function using a support grid of t
        matrices<-getMatrices(t)

        dett<-t[2]-t[1]
        inio0<-matrices$inio0
        Ur<-matrices$Ur
        Ur3<-matrices$Ur3
        basr<-matrices$basr
        infn1<-matrices$infn1

        fr<-function(y){
            x=c(y,inio0[-(1:bsz)])
            x[bsz+1]=(t((Ur[2,bsz+3]*Ur3[1,]-Ur[1,bsz+3]*Ur3[2,])*(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2])-
                        (Ur3[2,]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur3[3,])*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2]))%*%
                      x[-((bsz+1):(bsz+3))]-2*Ur[2,bsz+3]*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2]))/
                ((Ur[2,bsz+1]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+1])*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2])-
                 (Ur[1,bsz+1]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+1])*(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2]))
            x[bsz+2]=((Ur[2,bsz+3]*Ur3[3,]-Ur3[2,]*Ur[3,bsz+3])%*%x[-((bsz+1):(bsz+3))]-(Ur[2,bsz+1]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+1])*x[bsz+1]-2*Ur[2,bsz+3])/(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2])
            x[bsz+3]=(-Ur3[1,]%*%x[-((bsz+1):(bsz+3))]-Ur[1,bsz+1]*x[bsz+1]-Ur[1,bsz+2]*x[bsz+2])/Ur[1,bsz+3]

            v1=-sum(log((t(x)%*%basr)%*%t(infn1)*dett))+lmd*(((t(x)%*%basr)%*%t(t(x)%*%basr))*dett)
            return(v1)
        }

        evaluate<-function(coef1){
            ## evaluates the function on the support points
            ft1=(t(abs(t%*%t(rep(1,tsz))-rep(1,tsz)%*%t(t)))/2)%*%(t(basr)%*%coef1)*dett
            return(ft1)
        }

        get_coefficients<-function(par){

            coef1=inio0
            coef1[c(1:bsz)]=par
            coef1[bsz+1]=(t((Ur[2,bsz+3]*Ur3[1,]-Ur[1,bsz+3]*Ur3[2,])*(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2])- (Ur3[2,]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur3[3,])*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2]))%*% coef1[-((bsz+1):(bsz+3))]-2*Ur[2,bsz+3]*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2]))/ ((Ur[2,bsz+1]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+1])*(Ur[1,bsz+2]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+2])- (Ur[1,bsz+1]*Ur[2,bsz+3]-Ur[1,bsz+3]*Ur[2,bsz+1])*(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2]))
            coef1[bsz+2]=((Ur[2,bsz+3]*Ur3[3,]-Ur3[2,]*Ur[3,bsz+3])%*%coef1[-((bsz+1):(bsz+3))]-(Ur[2,bsz+1]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+1])*coef1[bsz+1]-2*Ur[2,bsz+3])/(Ur[2,bsz+2]*Ur[3,bsz+3]-Ur[2,bsz+3]*Ur[3,bsz+2])
            coef1[bsz+3]=(-Ur3[1,]%*%coef1[-((bsz+1):(bsz+3))]-Ur[1,bsz+1]*coef1[bsz+1]-Ur[1,bsz+2]*coef1[bsz+2])/Ur[1,bsz+3]
            return(coef1)
        }

        return(list("objective"=fr,"initial"=inio0[seq_len(bsz)],"evaluate"=evaluate,"get_coefficients"=get_coefficients,"basr"=basr,"lmd"=lmd))
    }
    return(getfunction)
}



optfn=function(Objective,linit,uinit,tsz){
    lind=uind=0
    t=seq(linit,uinit,length.out=tsz)

    while((lind!=1)|(uind!=tsz)){
##################optimize##############
#########################################


        ObjectiveFun=Objective(t)
        try1=stats::optim(ObjectiveFun$initial,ObjectiveFun$objective,control=list(maxit=10000000,abstol=1e-10,reltol=1e-10))

        coef1=ObjectiveFun$get_coefficients(try1$par)

##################update boundaries
        ft1=ObjectiveFun$evaluate(coef1)

        ## uind is the minimiser of the upper half of ft1
        ## lind is the minimiser of the lower half of ft1
        uind=round((tsz+max(which.min(ft1[which.max(ft1):tsz]))+which.max(ft1)-1)/2)
        lind=round((min(which.min(ft1[1:which.max(ft1)]))+1)/2)
        if((uind==tsz)&(lind==1)&(min(ft1)<0)){
            uind=tsz-1
            lind=2
        }
        u=t[uind]
        l=t[lind]

####################initials#########
#######################################
        t=seq(l,u,length.out=tsz)
    }
    return(list(t=t,u=u,l=l,basr=ObjectiveFun$basr,coef1=coef1,lmd0=ObjectiveFun$lmd))
}


get_init_limits<-function(ob,error){
### error can either be a sample of an empirical pure error distribution
### or a list the first elements is a distribution, and the remaining elements are parameters.
    lstart<-min(stats::density(ob)$x)
    ustart<-max(stats::density(ob)$x)

    if(is.list(error)){
### Named error distribution
### should make this more robust and
### insist on a named list of parameter values to allow better checking
        family=error[[1]]
        if(family=="Normal"){
            usub<-error[[2]]-stats::qnorm(0.9999)*error[[3]]  # Assumes first parameter is mean
            lsub<-error[[2]]+stats::qnorm(0.9999)*error[[3]] # second is std. dev.
        }else if(family=="Laplace"){
            usub<-error[[2]]-stats::qexp(0.9998)*error[[3]]  # Assumes first parameter is mean
            lsub<-error[[2]]+stats::qexp(0.9998)*error[[3]] # second is scale
        }else if(family=="Beta"){
### Note that this is a beta distribution on [0,1] scaled beta is not supported.
            usub<-stats::qbeta(0.0001,error[[2]],error[[3]])
            lsub<-stats::qbeta(0.9999,error[[2]],error[[3]])
        }else{
            stop(error=paste("Error family \"",family,"\" not yet implemented",sep=""))
        }
    }else{
        ## Not a great approach
        usub<-min(stats::density(error)$x)
        lsub<-max(stats::density(error)$x)
    }
    return(c(lstart-lsub,ustart-usub))
}


#' @export
pmledecon=function(ob,error,supp,n=1000,lmd=NULL,R=1e5,tsz=1000,stsz,bsz,subid=TRUE,conv=FALSE){
  ### error can either be a sample of an empirical pure error distribution
  ### or a list the first elements is a distribution, and the remaining elements are parameters.
  limits=get_init_limits(ob,error)
  linit=limits[1]
  uinit=limits[2]
  if(!missing(stsz)){
    ## stsz overrides tsz
    tsz=(uinit-linit)/stsz
  }
  t=seq(linit,uinit,length.out=tsz)

  ob=sort(ob)
  sz=length(ob)
  if(!missing(bsz)){
    if(sz>30){
      bsz=30
      sbsz=round(sz/10)
    }
    else{
      bsz=20
      sbsz=5
    }
  }
  subsp=round(seq(0,sz,length.out=bsz+1))
  each=subsp[2:(bsz+1)]-subsp[1:bsz]
  group=rep(1:bsz,times=each)
  data1=data.frame(1:sz,group)

  ### Counters for number of successful subsamples and number of unsuccessful subsamples.
  m1=0
  failed_subsamples=0


  ct=matrix(0,sbsz,tsz)
  lmdt=rep(0,sbsz)
  bs=matrix(0,tsz,sbsz)

  while(m1<sbsz){
    subind=unlist(splitstackshape::stratified(data1,"group",size=1)[,1])

    #        print("subbasis")
    #        print(subind)

    objective<-objective_function(ob,error,bsz,subind,lmd,R,t)

    result=try(optfn(objective,linit,uinit,tsz),silent=TRUE)
    if(inherits(result,"try-error")){
      ## Subsample failed.
      print(result)
      failed_subsamples=failed_subsamples+1
    }else{
      basr=result$basr
      coef1=result$coef1
      m1=m1+1
      ct[m1,]=result$t
      bs[,m1]=t(basr)%*%coef1
      lmdt[m1]=result$lmd0
      #            print(subind)
      #            print(paste("lambda=",result$lmd0,sep=""))
      #            print(bs[,m1])
    }
    if(subid){
      print(m1)
    }
  }
  l=min(ct[,1])
  u=max(ct[,tsz])
  if(missing(supp)){
    t.est=seq(l,u,length.out=n)
  }else{
    t.est=supp
    n=length(supp)
  }
  ft=matrix(0,sbsz,n)
  for(i in 1:sbsz){
    t=ct[i,]
    dett=t[2]-t[1]
    ft[i,]=(t(abs(t%*%t(rep(1,length(t.est)))-rep(1,length(t))%*%t(t.est)))/2)%*%(bs[,i])*dett
  }

  #    print(ft)
  mft=apply(ft,2,mean)

  #average_functions(ct,ft,t.est)
  if(conv==FALSE){
    convll=NULL
  }else{
    l1=NULL
    for(i in 1:sbsz){
      t=ct[i,]
      dett=t[2]-t[1]
      if(is.list(error)){
        if(error[[1]]=="Normal"){
          teinfn0=sort(ob)%*%t(rep(1,length(t)))-rep(1,sz)%*%t(t)
          #############dnorm
          teinfn1=(teinfn0-error[[2]])*stats::pnorm(teinfn0,error[[2]],error[[3]])+error[[3]]^2*stats::dnorm(teinfn0,error[[2]],error[[3]])
        }else if(error[[1]]=="Laplace"){
          teinfn0=sort(ob)%*%t(rep(1,length(t)))-rep(1,sz)%*%t(t)
          #############dlaplace
          teinfn1=error[[3]]*rmutil::plaplace(teinfn0,error[[2]],error[[3]])
          teinfn1[teinfn0>error[[2]]]=(teinfn0-error[[2]]+error[[3]]-error[[3]]*rmutil::plaplace(teinfn0,error[[2]],error[[3]]))[teinfn0>error[[2]]]
        }else if(error[[1]]=="Beta"){
          teinfn0=sort(ob)%*%t(rep(1,length(t)))-rep(1,sz)%*%t(t)
          #############dbeta
          teinfn1=teinfn0*stats::pbeta(teinfn0,error[[2]],error[[3]])-error[[2]]/(error[[2]]+error[[3]])*stats::pbeta(teinfn0,error[[2]]+1,error[[3]])
        }
      }else{
        teinfn0=-ob%*%t(rep(1,length(t)))+rep(1,sz)%*%t(t)
        teinfn1=apply(teinfn0,c(1,2),function(x){sum((x+error)[x+error>0])/(length(error))})
      }
      lik=t(bs[,i])%*%t(teinfn1)*dett
      lik[abs(lik)<1e-10]=1e-10
      l1.n=sum(log(lik))
      l1=c(l1,l1.n)
    }
    convll=mean(l1)
  }
  return(list(sup=t.est,f=mft,conll=convll,lmd.sub=lmdt))
}
