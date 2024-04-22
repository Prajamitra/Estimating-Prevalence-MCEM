
rm(list=ls())

Bootstrap_repli=500          # Number of bootstrap samples



########## Saarling's formula for computing factorial ###########

log_fact<-function(z){
val=max(z*log(z)-z,0)
val[val=="NaN"]=0
return(val)
}

log_fact=Vectorize(log_fact)

##########  Initial Values  ###############


alpha_1=0.2
alpha_2=0.1
alpha_3=0.2  
alpha_4=0.1

########### Arrays for results   ######################

LCI=c()
UCI=c()
cover=0
LCIC=c()
UCIC=c()
coverC=0

##### Real Data ##############

# WTC Data
X111=174
X110=88
X101=1658
X011=750
X100=1702
X010=270
X001=4323

#################   MCEM Estimation     ######################

X111_data=X111
X110_data=X110
X101_data=X101
X011_data=X011
X100_data=X100
X010_data=X010
X001_data=X001

X=array(0,3)

X[1]=X111+X110+X101+X100
X[2]=X111+X110+X001+X010
X[3]=X111+X101+X011+X001

X0=X111+X110+X101+X011+X100+X010+X001
X0_data=X0

X0_obs=X0


### Initialization for conditional expectation

n_hat = 13000 #X0+max(ceiling((X111/X110)*(X100/X101)*(X001/X011)*X010),1)

n_0 <- n_hat
alpha1_0 =alpha_1
alpha2_0 =alpha_2
alpha3_0 =alpha_3
alpha4_0 =alpha_4
  m1=0
  m2=0
  m3=0
  n1=0
  n2=0
  n3=0

HH=c(1,5)

while(abs(HH[length(HH)]-HH[length(HH)-1]) > 1) {

  n<-max(n_0,X0)
  alpha1<-alpha1_0
  alpha2<-alpha2_0
  alpha3<-alpha3_0
  alpha4<-alpha4_0
  alpha0 = alpha1+alpha2+alpha3+alpha4


n_obs=100  # No. of observations to be drawn

f<-function(x){
  N<-x[1] 
  a1<-x[2]
  a2<-x[3]
  a3<-x[4]
  a4<-x[5]
  a0<-a1+a2+a3+a4


if ( N > X0 && a1>0 && a2 >0 && a3>0 && a4 > 0 && a0 <1 ) {

  set.seed(1)

  g1<-rbeta(n_obs,m1+1,n1+1)
  g2<-rbeta(n_obs,m2+1,n2+1)
  g3<-rbeta(n_obs,m3+1,n3+1)

  
  P111=((1-alpha0)*g1*g2*g3)+(alpha1*g1*g3)+(alpha2*g1*g2)+(alpha3*g1*g2)+(alpha4*g1)
  P110=((1-alpha0)*g1*g2*(1-g3))+(alpha1*g1*(1-g3))
  P011=((1-alpha0)*(1-g1)*g2*g3)+(alpha2*(1-g1)*g2)
  P100=((1-alpha0)*g1*(1-g2)*(1-g3))+(alpha2*g1*(1-g2))
  P101=((1-alpha0)*g1*(1-g2)*g3)+(alpha3*g1*(1-g2))
  P010=((1-alpha0)*(1-g1)*g2*(1-g3))+(alpha3*(1-g1)*g2)
  P001=((1-alpha0)*(1-g1)*(1-g2)*g3)+(alpha1*(1-g1)*g3)
  P000=((1-alpha0)*(1-g1)*(1-g2)*(1-g3))+(alpha1*(1-g1)*(1-g3))+(alpha2*(1-g1)*(1-g2))+(alpha3*(1-g1)*(1-g2))+(alpha4*(1-g1))


  
L01= log_fact(N)+((X111*(1-alpha0)*g1*g2*g3/P111)*log((1-a0)*g1*g2*g3))+((X111*alpha1*g1*g3/P111)*log(a1*g1*g3))+(((X111*alpha2*g1*g2)/P111)*log(a2*g1*g2))+(((X111*alpha3*g1*g2)/P111)*log(a3*g1*g2))+(((X111*alpha4*g1)/P111)*log(a4*g1))
L02= (((X110*(1-alpha0)*g1*g2*(1-g3))/P110)*log((1-a0)*g1*g2*(1-g3)))+(((X110*alpha1*g1*(1-g3))/P110)*log(a1*g1*(1-g3)))
L03= (((X011*(1-alpha0)*(1-g1)*g2*g3)/P011)*log((1-a0)*(1-g1)*g2*g3))+(((X011*alpha2*(1-g1)*g2)/P011)*log(a2*(1-g1)*g2))
L04= (((X100*(1-alpha0)*g1*(1-g2)*(1-g3))/P100)*log((1-a0)*g1*(1-g2)*(1-g3)))+(((X100*alpha2*g1*(1-g2))/P100)*log(a2*g1*(1-g2)))
L05= (((X101*(1-alpha0)*g1*(1-g2)*g3)/P101)*log((1-a0)*g1*(1-g2)*g3))+(((X101*alpha3*g1*(1-g2))/P101)*log(a3*g1*(1-g2)))
L06= (((X010*(1-alpha0)*(1-g1)*g2*(1-g3))/P010)*log((1-a0)*(1-g1)*g2*(1-g3)))+(((X010*alpha3*(1-g1)*g2)/P010)*log(a3*(1-g1)*g2))
L07= (((X001*(1-alpha0)*(1-g1)*(1-g2)*g3)/P001)*log((1-a0)*(1-g1)*(1-g2)*g3))+(((X001*alpha1*(1-g1)*g3)/P001)*log(a1*(1-g1)*g3))
L08= ((((n-X0)*(1-alpha0)*(1-g1)*(1-g2)*(1-g3))/P000)*log((1-a0)*(1-g1)*(1-g2)*(1-g3)))+((((n-X0)*alpha1*(1-g1)*(1-g3))/P000)*log(a1*(1-g1)*(1-g3)))+((((n-X0)*alpha2*(1-g1)*(1-g2))/P000)*log(a2*(1-g1)*(1-g2)))+((((n-X0)*alpha3*(1-g1)*(1-g2))/P000)*log(a3*(1-g1)*(1-g2)))+((((N-X0)-(n-X0)*(1-alpha4*(1-g1)/P000)))*log(a4*(1-g1)))


L1=L01+L02+L03+L04+L05+L06+L07+L08 +  log(dbeta(g1,m1+1,n1+1)) + log(dbeta(g2,m2+1,n2+1)) + log(dbeta(g3,m3+1,n3+1))


## Monte varlo Expection of the denominator of log-likelihood based on y's

prob_y111=cbind(((1-alpha0)*g1*g2*g3),(alpha1*g1*g3),(alpha2*g1*g2),(alpha3*g1*g2),(alpha4*g1))/P111

y111_vec=array(NA, dim=c(n_obs,5))

for(i in 1:n_obs){
y111_vec[i,]=rmultinom(1,X111,prob_y111[i,])
}

y1111<<-y111_vec[,1]
y1112<<-y111_vec[,2]
y1113<<-y111_vec[,3]
y1114<<-y111_vec[,4]
y1110<<-y1111+y1112+y1113+y1114
y1101<<-rbinom(n_obs,X110,((1-alpha0)*g1*g2*(1-g3))/P110)   # if y1101 becomes 0, take it 0.00001
y0111<<-rbinom(n_obs,X011,((1-alpha0)*(1-g1)*g2*g3)/P011)
y1001<<-rbinom(n_obs,X100,((1-alpha0)*g1*(1-g2)*(1-g3))/P100)
y1011<<-rbinom(n_obs,X101,((1-alpha0)*g1*(1-g2)*g3)/P101)
y0101<<-rbinom(n_obs,X010,((1-alpha0)*(1-g1)*g2*(1-g3))/P010)
y0011<<-rbinom(n_obs,X001,((1-alpha0)*(1-g1)*(1-g2)*g3)/P001)

prob_y000=cbind(((1-alpha0)*(1-g1)*(1-g2)*(1-g3)),(alpha1*(1-g1)*(1-g3)),(alpha2*(1-g1)*(1-g2)),(alpha3*(1-g1)*(1-g2)),(alpha4*(1-g1)))/P000


y000_vec=array(NA, dim=c(n_obs,5))

for(i in 1:n_obs){
y000_vec[i,]=rmultinom(1,(n-X0),prob_y000[i,])
}

y0001<<-y000_vec[,1]
y0002<<-y000_vec[,2]
y0003<<-y000_vec[,3]
y0004<<-y000_vec[,4]
y0000<<-y0001+y0002+y0003+y0004

N00=N-X0-y0000
N00[N00<=0]=0

L2=log_fact(y1111)+log_fact(y1112)+log_fact(y1113)+log_fact(y1114)+log_fact(X111-y1110)+
log_fact(y1101)+log_fact(X110-y1101)+
log_fact(y0111)+log_fact(X011-y0111)+
log_fact(y1001)+log_fact(X100-y1001)+
log_fact(y1011)+log_fact(X101-y1011)+
log_fact(y0101)+log_fact(X010-y0101)+
log_fact(y0011)+log_fact(X001-y0011)+
log_fact(y0001)+log_fact(y0002)+log_fact(y0003)+log_fact(y0004)+log_fact(N00)
  
L_avg<- mean(L1-L2)  # this average is due to integration over P1, P2, P3


} else {
L_avg<- NaN
}

return(-L_avg)
} # end of function

R=optim(c(15000,alpha_1,alpha_2,alpha_3,alpha_4),f,method="BFGS",hessian = TRUE)

HH=c(HH, R$par[1])


n_0<-round(R$par[1])
  alpha1_0<-R$par[2]
  alpha2_0<-R$par[3]
  alpha3_0<-R$par[4]
  alpha4_0<-R$par[5]

m1=X[1]
n1=(n_0-X[1])

m2=y1111+y1113+y1114+y1101+X011+X010
n2=X100+X101+y0011+y0001+y0003+y0004

m3=y1111+y1112+y0111+y1011+X001
n3=X110+y1001+y0101+y0001+y0002

} # end of while






############# Bootstrap Analysis ##############


Sim_N_BS=array(NA, dim=Bootstrap_repli)
vec_MAE=c()
Sim_alpha1_BS=array(NA, dim=Bootstrap_repli)
Sim_alpha2_BS=array(NA, dim=Bootstrap_repli)
Sim_alpha3_BS=array(NA, dim=Bootstrap_repli)
Sim_alpha4_BS=array(NA, dim=Bootstrap_repli)
Sim_alpha_BS=array(NA, dim=Bootstrap_repli)



Sim_P1_BS=array(NA, dim=Bootstrap_repli)
Sim_P2_BS=array(NA, dim=Bootstrap_repli)
Sim_P3_BS=array(NA, dim=Bootstrap_repli)



for(s in 1:Bootstrap_repli){

### Bootstrap Data Generation ###########

BS_data=rmultinom(1,Sim_N,prob=c(X111_data/Sim_N,X110_data/Sim_N,X101_data/Sim_N,X011_data/Sim_N,X100_data/Sim_N,X010_data/Sim_N,X001_data/Sim_N,(Sim_N-X0_data)/Sim_N))
X111=BS_data[1,1]
X110=BS_data[2,1]
X101=BS_data[3,1]
X011=BS_data[4,1]
X100=BS_data[5,1]
X010=BS_data[6,1]
X001=BS_data[7,1]
X000=BS_data[8,1]

X=array(0,3)

X[1]=X111+X110+X101+X100
X[2]=X111+X110+X001+X010
X[3]=X111+X101+X011+X001

X0=X111+X110+X101+X011+X100+X010+X001

# Initialization
n_hat= round(rnorm(1,Sim_N,100)) #X0+max(ceiling((X111/X110)*(X100/X101)*(X001/X011)*X010),1)

n_0 <-n_hat  # to safeguard from n-X0 being 0 if generated X000=0
alpha1_0 =Sim_alpha1
alpha2_0 =Sim_alpha2
alpha3_0 =Sim_alpha3
alpha4_0 =Sim_alpha4
  m1=0
  m2=0
  m3=0
  n1=0
  n2=0
  n3=0

HH1=c(1,5)

while(abs(HH1[length(HH1)]-HH1[length(HH1)-1]) > 1) {

  n<-n_0
  alpha1<-alpha1_0
  alpha2<-alpha2_0
  alpha3<-alpha3_0
  alpha4<-alpha4_0
  alpha0 = alpha1+alpha2+alpha3+alpha4

n_obs=100  # No. of observations to be drawn

f_B<-function(x){
  N<-x[1] 
  a1<-x[2]
  a2<-x[3]
  a3<-x[4]
  a4<-x[5]
  a0<-a1+a2+a3+a4


if ( N > X0 && a1>0 && a2 >0 && a3>0 && a4 > 0 && a0 <1 ) {

  set.seed(1)

  g1<-rbeta(n_obs,m1+1,n1+1)
  g2<-rbeta(n_obs,m2+1,n2+1)
  g3<-rbeta(n_obs,m3+1,n3+1)


  
  P111=((1-alpha0)*g1*g2*g3)+(alpha1*g1*g3)+(alpha2*g1*g2)+(alpha3*g1*g2)+(alpha4*g1)
  P110=((1-alpha0)*g1*g2*(1-g3))+(alpha1*g1*(1-g3))
  P011=((1-alpha0)*(1-g1)*g2*g3)+(alpha2*(1-g1)*g2)
  P100=((1-alpha0)*g1*(1-g2)*(1-g3))+(alpha2*g1*(1-g2))
  P101=((1-alpha0)*g1*(1-g2)*g3)+(alpha3*g1*(1-g2))
  P010=((1-alpha0)*(1-g1)*g2*(1-g3))+(alpha3*(1-g1)*g2)
  P001=((1-alpha0)*(1-g1)*(1-g2)*g3)+(alpha1*(1-g1)*g3)
  P000=((1-alpha0)*(1-g1)*(1-g2)*(1-g3))+(alpha1*(1-g1)*(1-g3))+(alpha2*(1-g1)*(1-g2))+(alpha3*(1-g1)*(1-g2))+(alpha4*(1-g1))
  
L01= log_fact(N)+((X111*(1-alpha0)*g1*g2*g3/P111)*log((1-a0)*g1*g2*g3))+((X111*alpha1*g1*g3/P111)*log(a1*g1*g3))+(((X111*alpha2*g1*g2)/P111)*log(a2*g1*g2))+(((X111*alpha3*g1*g2)/P111)*log(a3*g1*g2))+(((X111*alpha4*g1)/P111)*log(a4*g1))
L02= (((X110*(1-alpha0)*g1*g2*(1-g3))/P110)*log((1-a0)*g1*g2*(1-g3)))+(((X110*alpha1*g1*(1-g3))/P110)*log(a1*g1*(1-g3)))
L03= (((X011*(1-alpha0)*(1-g1)*g2*g3)/P011)*log((1-a0)*(1-g1)*g2*g3))+(((X011*alpha2*(1-g1)*g2)/P011)*log(a2*(1-g1)*g2))
L04= (((X100*(1-alpha0)*g1*(1-g2)*(1-g3))/P100)*log((1-a0)*g1*(1-g2)*(1-g3)))+(((X100*alpha2*g1*(1-g2))/P100)*log(a2*g1*(1-g2)))
L05= (((X101*(1-alpha0)*g1*(1-g2)*g3)/P101)*log((1-a0)*g1*(1-g2)*g3))+(((X101*alpha3*g1*(1-g2))/P101)*log(a3*g1*(1-g2)))
L06= (((X010*(1-alpha0)*(1-g1)*g2*(1-g3))/P010)*log((1-a0)*(1-g1)*g2*(1-g3)))+(((X010*alpha3*(1-g1)*g2)/P010)*log(a3*(1-g1)*g2))
L07= (((X001*(1-alpha0)*(1-g1)*(1-g2)*g3)/P001)*log((1-a0)*(1-g1)*(1-g2)*g3))+(((X001*alpha1*(1-g1)*g3)/P001)*log(a1*(1-g1)*g3))
L08= ((((n-X0)*(1-alpha0)*(1-g1)*(1-g2)*(1-g3))/P000)*log((1-a0)*(1-g1)*(1-g2)*(1-g3)))+((((n-X0)*alpha1*(1-g1)*(1-g3))/P000)*log(a1*(1-g1)*(1-g3)))+((((n-X0)*alpha2*(1-g1)*(1-g2))/P000)*log(a2*(1-g1)*(1-g2)))+((((n-X0)*alpha3*(1-g1)*(1-g2))/P000)*log(a3*(1-g1)*(1-g2)))+((((N-X0)-(n-X0)*(1-alpha4*(1-g1)/P000)))*log(a4*(1-g1)))


L1=L01+L02+L03+L04+L05+L06+L07+L08 +  log(dbeta(g1,m1+1,n1+1)) + log(dbeta(g2,m2+1,n2+1)) + log(dbeta(g3,m3+1,n3+1))





## Monte carlo Expection of the denominator of log-likelihood based on y's

prob_y111=cbind(((1-alpha0)*g1*g2*g3),(alpha1*g1*g3),(alpha2*g1*g2),(alpha3*g1*g2),(alpha4*g1))/P111


y111_vec=array(NA, dim=c(n_obs,5))
for(i in 1:n_obs){
y111_vec[i,]=rmultinom(1,X111,prob_y111[i,])
}

y1111<<-y111_vec[,1]
y1112<<-y111_vec[,2]
y1113<<-y111_vec[,3]
y1114<<-y111_vec[,4]
y1110<<-y1111+y1112+y1113+y1114


y1101<<-rbinom(n_obs,X110,((1-alpha0)*g1*g2*(1-g3))/P110)   # if y1101 becomes 0, take it 0.00001
y0111<<-rbinom(n_obs,X011,((1-alpha0)*(1-g1)*g2*g3)/P011)
y1001<<-rbinom(n_obs,X100,((1-alpha0)*g1*(1-g2)*(1-g3))/P100)
y1011<<-rbinom(n_obs,X101,((1-alpha0)*g1*(1-g2)*g3)/P101)
y0101<<-rbinom(n_obs,X010,((1-alpha0)*(1-g1)*g2*(1-g3))/P010)
y0011<<-rbinom(n_obs,X001,((1-alpha0)*(1-g1)*(1-g2)*g3)/P001)

prob_y000<<-cbind(((1-alpha0)*(1-g1)*(1-g2)*(1-g3)),(alpha1*(1-g1)*(1-g3)),(alpha2*(1-g1)*(1-g2)),(alpha3*(1-g1)*(1-g2)),(alpha4*(1-g1)))/P000

y000_vec=array(NA, dim=c(n_obs,5))

for(i in 1:n_obs){
y000_vec[i,]=rmultinom(1,(n-X0),prob_y000[i,])
}

y0001<<-y000_vec[,1]
y0002<<-y000_vec[,2]
y0003<<-y000_vec[,3]
y0004<<-y000_vec[,4]
y0000<<-y0001+y0002+y0003+y0004

N00=N-X0-y0000
N00[N00<=0]=1e-6

L2=log_fact(y1111)+log_fact(y1112)+log_fact(y1113)+log_fact(y1114)+log_fact(X111-y1110)+
log_fact(y1101)+log_fact(X110-y1101)+
log_fact(y0111)+log_fact(X011-y0111)+
log_fact(y1001)+log_fact(X100-y1001)+
log_fact(y1011)+log_fact(X101-y1011)+
log_fact(y0101)+log_fact(X010-y0101)+
log_fact(y0011)+log_fact(X001-y0011)+
log_fact(y0001)+log_fact(y0002)+log_fact(y0003)+log_fact(y0004)+log_fact(N00)
  
L_avg<- mean(L1-L2)  # this average is due to integration over P1, P2, P3

} else {

L_avg<- NaN
}

return(-L_avg)

} # end of function

EE=try(optim(c(n_hat,Sim_alpha1,Sim_alpha2,Sim_alpha3,Sim_alpha4),f_B, method="BFGS",hessian = TRUE), silent=TRUE)
if ('try-error' %in% class(EE)){s=s
} else {
R1=EE
}

HH1=c(HH1,R1$par[1])

n_0<-round(R1$par[1])
  alpha1_0<-R1$par[2]
  alpha2_0<-R1$par[3]
  alpha3_0<-R1$par[4]
  alpha4_0<-R1$par[5]

m1=X[1]
n1=(n_0-X[1])


m2=y1111+y1113+y1114+y1101+X011+X010
n2=X100+X101+y0011+y0001+y0003+y0004

m3=y1111+y1112+y0111+y1011+X001
n3=X110+y1001+y0101+y0001+y0002

} # end of while in Bootstrap

Sim_N_BS[s]=R1$par[1]

vec_MAE=c(vec_MAE,abs(Sim_N_BS[s]-Sim_N)) #estimate of mean absolute error based on bootstrap estimates


Sim_alpha1_BS[s]=R1$par[2]
Sim_alpha2_BS[s]=R1$par[3]
Sim_alpha3_BS[s]=R1$par[4]
Sim_alpha4_BS[s]=R1$par[5]
Sim_alpha_BS[s]=R1$par[2]+R1$par[3]+R1$par[4]+R1$par[5]


Sim_P1_BS[s]=m1/(n1+m1)
Sim_P2_BS[s]=mean(m2/(n2+m2))
Sim_P3_BS[s]=mean(m3/(n3+m3))




set.seed(NULL)

print(s)

} # end of s loop for bootstrap

######## Bootstrap Ends #######################

sd_N_hat_BS=sd(Sim_N_BS)

Medi=median(Sim_N_BS)
MAE_N_hat_BS=sum(abs(Sim_N_BS-Medi))/length(Sim_N_BS)

C_BS= exp(1.96*sqrt(log(1+(sd_N_hat_BS/(Sim_N-X0_obs))^2)))

LCI_Chao=(X0_obs+(Sim_N-X0_obs)/C_BS)

UCI_Chao=(X0_obs+(Sim_N-X0_obs)*C_BS)



################## Results  ###############


#########  Checking convergence  ########


R$convergence  # 0 means convergence achived


##### Estimate of N from real data ###############

Sim_N=round(R$par[1])  

Sim_N # Point estimate


LCI_Chao # Lower confidence level
UCI_Chao # Upper confidence level


sd_N_hat_BS # SE of estimator of N

MAE_N_hat_BS # MAE of estimator of N


##### Estimate of dependence parameters from real data ###############

Sim_alpha1=R$par[2]
Sim_alpha2=R$par[3]
Sim_alpha3=R$par[4]
Sim_alpha4=R$par[5]

Sim_alpha1 
Sim_alpha2
Sim_alpha3
Sim_alpha








