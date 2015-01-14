% t_math2 and 3 are very nearly the same. t_math2 is computed more quickly.
close all
clear

u0 = -16;
ua = 0;
k_cond = 0.5270;
w_perf = 6;
Power = 12;
mua = 1473.54;
mus = 8000;
r1 = 0.00075;
r2 = 1;
anfact = 0.88;
cblood = 3840;
mutr  = mua+mus*(1-anfact);
mueff = sqrt(3.0*mua*mutr);
r_var = linspace(0.00075, 0.1, 10000);

u0_new = u0 + 310.15;
ua_new = ua + 310.15;
% u0_new = u0 + 37;
% ua_new = ua + 37;

t_trace1 = zeros(1000,1);
t_math1 = t_trace1;
t_math2 = t_trace1;
t_math3 = t_trace1;
tic
for ii = 1:1000;
    
    t_trace1 (ii) = old1D(u0,ua,k_cond,w_perf,Power,r_var(ii),mua,mus,r1,r2,anfact);
    %t_math1 (ii)  = mathematicaSSS5(u0,ua,k_cond,w_perf,Power,r_var(ii),mua,mus,r1,r2,anfact,cblood);
    %t_math2 (ii) = change_AnalyticSSS_mod5(u0,ua,k_cond,w_perf,Power,r_var(ii),mua,mus,r1,r2,anfact,cblood);
    
    t_math3 (ii) = mathematicaSSS6(u0,ua,k_cond,w_perf,Power,r_var(ii),mueff,r1,r2,cblood);
    
    %t_math1 (ii)  = mathematicaSSS2(r_var(ii),mua,mus,anfact);
    %t_math2 (ii)  = mathematicaSSS3(r_var(ii),mua,mus,anfact);
    t_math2 (ii)  = mathematicaSSS_symbolicToolbox(u0,ua,k_cond,w_perf,Power,r_var(ii),mua,mus,r1,r2,anfact,cblood);

end
toc
t_math_diff1 =t_math1-t_math2;
t_math_diff2 =t_trace1-t_math1;
t_math_diff3 =t_math3-t_math2;

figure(1); plot(r_var(1:1000), t_trace1);
figure(2); plot(r_var(1:1000), t_math1);
figure(3); plot(r_var(1:1000), t_math2);
figure(4); plot(r_var(1:1000), t_math3);
figure(5); plot(r_var(1:1000), t_math_diff1);
figure(6); plot(r_var(1:1000), t_math_diff2);
figure(7); plot(r_var(1:1000), t_math_diff3);