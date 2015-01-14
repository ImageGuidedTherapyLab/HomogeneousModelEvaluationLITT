function [aa]= mathematicaSSS_symbolicToolbox(u0,ua,kCond,wPerf,Powr,rVar,mua,mus,r1,r2,anfact,cblood);

mutr  = mua+mus*(1.0-anfact);
mueff = sqrt(3.0*mua*mutr);

BB = sqrt(cblood*wPerf/kCond);
CC = wPerf*cblood;
DD = kCond*mueff^2;

% Testing the symbolic toolbox
syms bb cc dd rr MuEff R1 R2 UA U0 POWER

BBB = sym(BB);
CCC = sym(CC);
DDD = sym(DD);
rrr = sym(rVar);
mumueff = sym(mueff);
RR1 = sym(r1);
RR2 = sym(r2);
uaua = sym(ua);
u0u0 = sym(u0);
PowPow = sym(Powr);

% ff = sym('exp(a+b)')
% syms a b
% c = sym(CC)
% d = sym(DD)
% vpa(subs(ff,[a, b],[c, d] ) )

ff_temperature = sym(' (exp(-bb*rr-MuEff*(R1+R2+rr))*(-exp(bb*(R1+2*R2)+MuEff*(R2+rr))*MuEff^2*POWER*(-1+bb*R2)+exp(MuEff*(R1+R2)+bb*(2*R2+rr))*MuEff^2*POWER*(-1+bb*R2)+exp(MuEff*(R1+R2)+bb*(2*R1+rr))*MuEff^2*POWER*(1+bb*R2)-exp(MuEff*(R2+rr)+bb*(R1+2*rr))*MuEff^2*POWER*(1+bb*R2)-exp(bb*(2*R1+R2)+MuEff*(R1+rr))*MuEff^2*POWER*(1+MuEff*R2)+exp(MuEff*(R1+rr)+bb*(R2+2*rr))*MuEff^2*POWER*(1+MuEff*R2)+4*(cc-dd)*exp(bb*(R1+2*R2)+MuEff*(R1+R2+rr))*pi*R1*(-1+bb*R2)*U0+4*(cc-dd)*exp(MuEff*(R1+R2+rr)+bb*(R1+2*rr))*pi*R1*(1+bb*R2)*U0+4*(cc-dd)*exp(MuEff*(R1+R2+rr)+bb*(2*R2+rr))*pi*(-1+bb*R2)*rr*UA+4*(cc-dd)*exp(bb*(2*R1+rr)+MuEff*(R1+R2+rr))*pi*(1+bb*R2)*rr*UA))/(4*(cc-dd)*pi*(exp(2*bb*R2)*(-1+bb*R2)+exp(2*bb*R1)*(1+bb*R2))*rr) ');

aa = vpa( subs( ff_temperature, [bb, cc, dd, rr, MuEff, R1, R2, UA, U0, POWER], [BBB, CCC, DDD, rrr, mumueff, RR1, RR2, uaua, u0u0, PowPow]) );


% ff = sym('exp(a)')
% b = sym('3000')
% vpa(subs(ff,'a',b))
% 
% ff = sym('exp(a+b)')
% syms a b
% c = sym('1')
% d = sym('3000')
% vpa(subs(ff,[a, b],[c, d] ) )

end