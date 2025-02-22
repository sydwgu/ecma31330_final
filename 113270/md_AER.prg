/* Minimum distance estimation                                      */

new;

ma          =1;								/* Choose if you want MA(1) for income                         	*/
taste       =1;								/* Choose if you want taste shocks          				   	*/
varying_ins =1;								/* Choose if you want time-varying insurance coeff				*/

library optmum;
#include optmum.ext;
optset;

open f1=c:\papers\partial\revise2\coh7992;
x      =readr(f1,rowsf(f1));
f1     =close(f1);

m      =x[.,1];                       /*Vector of covariances*/
var    =x[.,2:cols(x)-1];             /*Weighting matrix - Variance matrix for optimal mindist*/

T           =x[1,cols(x)];            /*Lenght of the panel (in FD)*/
initial_year=x[3,cols(x)];
final_year  =x[4,cols(x)];

T_pz=1+varying_ins;
T_pe=1+varying_ins;
T_e =T-2;
T_z =T-4;
T_u=9;

extra=ma+taste;           /* Number of extra parameters*/


@ Starting parameters and model type@

    b0                                       			=zeros(T_e+T_z+T_pe+T_pz+extra+T_u,1);
    b0[extra-1]                              			=-0.1;                 /*teta*/
    b0[extra]                                			=0.01;                 /*var csi*/
    b0[1+extra:T_z+extra]                    			=0.03*ones(T_z,1);     /*variance of the perm shocks*/
    b0[1+T_z+extra:T_z+T_e+extra]            			=0.03*ones(T_e,1);     /*variance of the trans shocks*/
    b0[1+T_z+T_e+extra:T_z+T_e+T_pz+extra]   			=1*ones(T_pz,1);       /*coeff of partial insurance perm shock*/
    b0[1+T_z+T_e+T_pz+extra:T_z+T_e+T_pz+T_pe+extra]	=0*ones(T_pe,1);       /*coeff of partial insurance of trans shock*/
    b0[1+T_z+T_e+T_pz+T_pe+extra:rows(b0)]    			=0.06*ones(T_u,1);     /*variance of meas. error consumption*/
	
__title = "Error-Components Model";
__row=0;
__miss=1;

{b,fun,grad,ok} = optmum(&cov,b0);                          @ Optimisation procedure @
{dfb}           = gradfd(&fcn,b);                           @ Gradient of covariances at minimum - df/db @

/*Choose weighting matrix*/

ww=invpd(diagrv(zeros(rows(var),cols(var)),diag(var)));
covth=invpd(dfb'*ww*dfb)*dfb'*ww*var*ww*dfb*invpd(dfb'*ww*dfb); @ Covariances of parameters under EWMD @

fmn  =fcn(b);                                               @ Value of covariances at minimum @
df   =rows(m) - rank(dfb);                                  @ Degrees of freedom @
seth =sqrt(diag(covth));                                    @ Standard Errors @

A=dfb*invpd(dfb'*ww*dfb)*dfb'*ww;
W=eye(rows(A))-A;
mf=(m-fmn)'*pinv(W*var*W')*(m-fmn);
prb  =cdfchic(mf,df);                                       @ P-value of the above test @

insur_est=2~2;

results=(b~seth)|(mf~prb)|insur_est;
y=saved(results,"c:\\papers\\partial\\revise2\\MD_res",0);

output file=c:\papers\partial\revise2\batchres.log on;

print "teta    "   b[extra-1]~seth[extra-1];
print "varcsi  "   b[extra]~seth[extra];
print "zt      "   b[1+extra:T_z+extra]~seth[1+extra:T_z+extra]; 
print "et      "   b[1+T_z+extra:T_z+T_e+extra]~seth[1+T_z+extra:T_z+T_e+extra];
print "phit    "   b[1+T_z+T_e+extra:T_z+T_e+T_pz+extra]~seth[1+T_z+T_e+extra:T_z+T_e+T_pz+extra]; 
print "psit    "   b[1+T_z+T_e+T_pz+extra:T_e+T_z+T_pe+T_pz+extra]~seth[1+T_z+T_e+T_pz+extra:T_e+T_z+T_pe+T_pz+extra];
print "varv    "   b[1+T_z+T_e+T_pz+T_pe+extra:T_e+T_z+T_pe+T_pz+extra+T_u]~seth[1+T_z+T_e+T_pz+T_pe+extra:T_e+T_z+T_pe+T_pz+extra+T_u];
print "Chi-sq. "   mf~df~prb;

if varying_ins==1;
	tttphi=((b[T_z+T_e+T_pz+extra]-b[1+T_z+T_e+extra])^2)/(covth[T_z+T_e+T_pz+extra,T_z+T_e+T_pz+extra]
			+covth[1+T_z+T_e+extra,1+T_z+T_e+extra]-2*covth[T_z+T_e+T_pz+extra,1+T_z+T_e+extra]);						/*test constant phi*/
	cdfchic(tttphi,1);
	tttpsi=((b[T_e+T_z+T_pe+T_pz+extra]-b[1+T_z+T_e+T_pz+extra])^2)/(covth[T_e+T_z+T_pe+T_pz+extra,T_e+T_z+T_pe+T_pz+extra]
			+covth[1+T_z+T_e+T_pz+extra,1+T_z+T_e+T_pz+extra]-2*covth[1+T_z+T_e+T_pz+extra,T_e+T_z+T_pe+T_pz+extra]);	/*test constant psi*/
	cdfchic(tttpsi,1);
endif;

output off;

@system;@
end;

/*NOTE: BELOW YOU WILL FIND THE MINIMUM DISTANCE PROC(EDURE) IN CASE YOU WANT TO MODIFY IT
  or ACCESS IT*/

@ Minimum Distance Procedure @
proc cov(b);
local fmm,f,peso;
fmm=fcn(b);
peso=diagrv(zeros(rows(var),cols(var)),diag(var));
f=(m-fmm)'*invpd(peso)*(m-fmm);     @ m is vector of actual covariances - fmm is implied covariances @
retp(f);
endp;

@ Gradient Procedure  - defines implied covariance structure from hypothesised model @

proc fcn(b);

local teta,psit,phit,zt,et,lev,levy,levc,levyc,levaux1,levaux2,lev1,
j,j1,i,difa1,difa2,fm,levt,levp,L,k,r,varv,
dify,difcd,difc,difcme,difyc,dif,varcsi,vareta,missing_v;

k=extra;

teta   =b[k-1];	
vareta =0;
varcsi =b[k];
zt     =b[1+k:T_z+k]; 
et     =b[1+T_z+k:T_z+T_e+k];
phit   =b[1+T_z+T_e+k:T_z+T_e+T_pz+k]; 
psit   =b[1+T_z+T_e+T_pz+k:T_e+T_z+T_pe+T_pz+k];
varv   =b[1+T_z+T_e+T_pz+T_pe+k:T_e+T_z+T_pe+T_pz+k+T_u];

dify  =zeros(T,T);
difcd =zeros(T,T); /* Consumption */
difc  =zeros(T,T); /* Consumption */
difcme=zeros(T,T); /* Measurement error of consumption */
difyc =zeros(T,T);
dif   =zeros(2*T,2*T);

/* This is the variance of Income */

dify[1,1]=zt[1]+et[1]+(1-teta)^2*et[1]+teta^2*et[1];
dify[2,2]=zt[1]+et[2]+(1-teta)^2*et[1]+teta^2*et[1];
dify[3,3]=zt[1]+et[3]+(1-teta)^2*et[2]+teta^2*et[1];
j=4;
do while j<=T-3;
    dify[j,j]=zt[j-2]+et[j]+(1-teta)^2*et[j-1]+teta^2*et[j-2];
j=j+1;
endo;          
dify[T-2,T-2]=zt[T-4]+et[T-2]+(1-teta)^2*et[T-3]+teta^2*et[T-4];
dify[T-1,T-1]=zt[T-4]+et[T-2]+(1-teta)^2*et[T-2]+teta^2*et[T-3];
dify[T,T]    =zt[T-4]+et[T-2]+(1-teta)^2*et[T-2]+teta^2*et[T-2];

dify[1,2]=-(1-teta)^2*et[1];
j=3;
do while j<=T-1;
    dify[j-1,j]=-(1-teta)*et[j-1]+teta*(1-teta)*et[j-2];
j=j+1;
endo;
dify[T-1,T-2]=-(1-teta)^2*et[T-2];

j=3;
do while j<=T;
    dify[j-2,j]=-teta*et[j-2];
j=j+1;
endo;

i=2;
do while i<=T;
j=i;
do while j<=T;
    dify[j,i-1]=dify[i-1,j];
j=j+1;
endo;
i=i+1;
endo;

/* This is the variance of Consumption */

difcd[1,1]=phit[1]^2*zt[1]+psit[1]^2*et[1]+varcsi;
difcd[2,2]=phit[1]^2*zt[1]+psit[1]^2*et[2]+varcsi;
difcd[3,3]=phit[1]^2*zt[1]+psit[1]^2*et[3]+varcsi;
j=4;
do while j<=6;
    difcd[j,j]=phit[1]^2*zt[j-2]+psit[1]^2*et[j]+varcsi;
j=j+1;
endo;

j=7;
do while j<=T-3;
    difcd[j,j]=phit[1+varying_ins]^2*zt[j-2]+psit[1+varying_ins]^2*et[j]+varcsi;
j=j+1;
endo;

difcd[T-2,T-2]=phit[1+varying_ins]^2*zt[T-4]+psit[1+varying_ins]^2*et[T-2]+varcsi;
difcd[T-1,T-1]=phit[1+varying_ins]^2*zt[T-4]+psit[1+varying_ins]^2*et[T-2]+varcsi;
difcd[T,T]    =phit[1+varying_ins]^2*zt[T-4]+psit[1+varying_ins]^2*et[T-2]+varcsi;

missing_v=(varv[1]+varv[2]+varv[3]+varv[4]+varv[5]+varv[6]+varv[7]+varv[8]+varv[9])/9;
difcme[1,1]=2*varv[1];
j=2;
do while j<=7;
 difcme[j,j]=varv[j]+varv[j-1];
j=j+1;
endo;
difcme[8,8]=varv[7]+missing_v;
j=9;
do while j<=11;
 difcme[j,j]=2*missing_v;
j=j+1;
endo;
difcme[12,12]=varv[8]+missing_v;
difcme[13,13]=varv[9]+varv[8];
difcme[14,14]=2*varv[9];


j=1;
do while j<=7;
 difcme[j,j+1]=-varv[j];
j=j+1;
endo;
j=8;
do while j<=11;
 difcme[j,j+1]=-missing_v;
j=j+1;
endo;
j=12;
do while j<=13;
 difcme[j,j+1]=-varv[j-4];
j=j+1;
endo;

difc=difcme+difcd;

i=2;
do while i<=T;
j=i;
do while j<=T;
    difc[j,i-1]=difc[i-1,j];
j=j+1;
endo;
i=i+1;
endo;

/* This is the Covariance of Income and Consumption */

difyc[1,1]=phit[1]*zt[1]+psit[1]*et[1];
difyc[2,2]=phit[1]*zt[1]+psit[1]*et[2];
difyc[3,3]=phit[1]*zt[1]+psit[1]*et[3];

j=4;
do while j<=6;
    difyc[j,j]=phit[1]*zt[j-2]+psit[1]*et[j];
j=j+1;
endo;

j=7;
do while j<=T-3;
    difyc[j,j]=phit[1+varying_ins]*zt[j-2]+psit[1+varying_ins]*et[j];
j=j+1;
endo;

difyc[T-2,T-2]=phit[1+varying_ins]*zt[T-4]+psit[1+varying_ins]*et[T-2];
difyc[T-1,T-1]=phit[1+varying_ins]*zt[T-4]+psit[1+varying_ins]*et[T-2];
difyc[T,T]    =phit[1+varying_ins]*zt[T-4]+psit[1+varying_ins]*et[T-2];

j=2;
do while j<=7;
    difyc[j-1,j]=-(1-teta)*psit[1]*et[j-1];
j=j+1;
endo;

j=8;
do while j<=T-1;
    difyc[j-1,j]=-(1-teta)*psit[1+varying_ins]*et[j-1];
j=j+1;
endo;

difyc[T-1,T]=-(1-teta)*psit[1+varying_ins]*et[T-2];

j=3;
do while j<=8;
    difyc[j-2,j]=-teta*psit[1]*et[j-2];
j=j+1;
endo;

j=9;
do while j<=T-1;
    difyc[j-2,j]=-teta*psit[1+varying_ins]*et[j-2];
j=j+1;
endo;

/* Final matrix */

dif[1:T,1:T]            =difc;
dif[T+1:2*(T),1:T]      =difyc;
dif[1:T,T+1:2*(T)]      =difyc';
dif[T+1:2*(T),T+1:2*(T)]=dify;

difa1=dif[1:8,.]|dif[12:2*T,.];			/*Accounts for the fact that PSID did not collect food data in some years*/
difa2=difa1[.,1:8]~difa1[.,12:2*T];

fm=vech(difa2);
retp(fm);
endp;
end;
