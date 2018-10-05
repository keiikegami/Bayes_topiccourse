#include<oxstd.h>
#include<oxprob.h>
#include<oxfloat.h>
#include<oxdraw.h>
#import<maximize>

main(){
decl vx, va, vb;
va = (1,2,3);
vb = (0.1,0.2,0.3);

vx = rangig(3,1,0.5,va,vb);
print(vx);
}