u0 = -16;
ua = 0;
k = 0.527;
w = 6;
P = 1.185;
mueff=146.942;
r1=0.00075;
r2=1;
cblood=3840;

r = linspace(r1+0.00001, 0.1, 1000);
t_sample = zeros(1, length(r));
for j=1:length(r)
    
    t_sample(j)=mathematicaSSS6(u0,ua,k,w,P,r(j),mueff,r1,r2,cblood);
end