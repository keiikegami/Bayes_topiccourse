#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import<maximize>

main(){

decl nT,np,ns,vy,mX,mdata;
decl i,nburn;
decl vb,dphi,mb,vphi, dlamsq, mtauinv;
decl vm0,mC0,da0,db0, da1,db1;
decl vmu,mOm,mS;

decl vtau,dlam;

mdata = loadmat("simdata.csv");

vy = mdata[][0];
mX = mdata[][1:];

nT = sizer(mX);
np = sizec(mX);

// prior
vm0 = ones(np,1);
mC0 = unit(np);
da0 = 1;
db0 = 1;
// lambdaのshape parameterは大きくないとだめ
// OK : 12,15, 20, 25, 50		（ちっちゃい方が良さげ）
// NO : 10
da1 = 30;
db1 = 1;

// initialize
vb = vm0;
dphi = da0/db0;
vtau = rangamma(np,1,ones(np,1),0.5*ones(np,1));
dlamsq = 1;

// Simulation size
ns = 1000;
nburn = 0.1*ns;

// storage
mb = zeros(ns,np);
vphi = zeros(ns,1);

// use these in MCMC
	   
// START MCMC
// tau and lambda are not our interest so these are not stored in this MCMC.
// lamsqのaが結構重要っぽい
// tauに一個でもNaNが含まれると全部NaNになる
// tauから始めるのがよくないのか？
for(i=-nburn;i<ns;i++){
mtauinv = diag(1 ./ vtau);
vmu = invert(mX'*mX + mtauinv)*mX'*vy;
mS = choleski(invert(dphi*(mX'*mX + mtauinv)));
vb = vmu + mS * rann(np, 1);
dphi = rangamma(1,1,0.5*(np+nT) + da0, 0.5*(vy - mX*vb)'*(vy - mX*vb) + 0.5*vb'*mtauinv*vb + db0);
dlamsq = rangamma(1,1,da1,vtau'*ones(np,1)/2 + db1);
vtau = rangig(np,1,da0-1/2,fabs(sqrt(dphi)*vb),dlamsq*ones(np,1));
//vtau = rangig(np,1,0.5,dlamsq*ones(np,1),dphi*vb.*vb);

if(i == -nburn+4){
print(vtau);
print(dlamsq);
print(vmu);
print(mS);
print(vb);
print(dphi);
print(dphi*vb.*vb);
print(dlamsq*ones(np,1));
}

if(i>=0){
mb[i][] = vb';
vphi[i] = dphi;
}

}

print(vphi);

// Summarize
	for(i=0;i<np;i++){
	DrawDensity(i,mb[][i]',sprint("Coef. No.",i+1),1,1);  DrawAdjust(ADJ_COLOR,2,10); 	DrawTitle(i, "");												
	}
	SaveDrawWindow(sprint("beta.eps"));
  	CloseDrawWindow();

	DrawDensity(0,1 ./ sqrt(vphi)',"$\\sigma  = \\phi ^{-1/2}$",1,1);  DrawAdjust(ADJ_COLOR,2,10); 	DrawTitle(i, "");												
	SaveDrawWindow(sprint("sigma2.eps"));
  	CloseDrawWindow();

// Convergence diagnosis by sample path

	DrawT(0, vphi',1,1,1);		DrawAdjust(ADJ_COLOR,2,12);		DrawTitle(0, "Path of $\\phi$");												
	SaveDrawWindow(sprint("path_phi.eps"));
  	CloseDrawWindow();


}