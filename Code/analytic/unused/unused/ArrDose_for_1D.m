%Dose function uses Arrhenius dose to calculate the thermal dose, w.  A w>=1
%means the voxel is judged as dead. This function takes temperature maps with uncertainty as
%input and outputs Arrhenius dose with uncertainty at the 2 sigma level.

% 'w' is the output, (x-dim,y-dim,z-dim). 'T' is the temperature image, (x-dim,y-dim,z-dim,time). 'w.mean_time' finds the differential dose accrued
% in one time step, 'dt'.  'Iso', (x-dim,y-dim,z-dim) finds the dead voxels.
function [w,Iso]=ArrDose_for_1D(T,dt);
T=double(T);


A=3.1*10^98;   %(1/s)  Constants and convert T from degrees Celsius to Kelvin
E=6.28*10^5;    %(J/mol)  Reference: Fahrenholtz et al. May/June 2013 in IJH
R=8.314;        %(J/mol/K)
%dt=5;           %(s)    Perhaps should be an input because T-maps could have different time steps
T1=T+273.15;    %Convert deg C into K

a=size(T);            %Initialize w1 and Iso
w=zeros(a(1),a(2),a(3),a(4));
Iso=zeros(a(1),a(2),a(3));

w=A*exp(-(E/(R*T1)));  %This is the discretized Arrhenius dose, no time yet
w=single(sum(w,4)*dt); %This is the summation part of discretized dose

%w=squeeze(w);

Iso( w>1 ) = 1;  %Make binary live/dead maps

Iso=logical(Iso);

end

% w1=A.*exp(-(E./(R.*T1))).*dt;  %This is the discretized Arrhenius dose with the uncertainties
% w_plus=A.*exp(-(E./(R.*T_plus))).*dt;
% w_minus=A.*exp(-(E./(R.*T_minus))).*dt;
% w.mean=sum(w1,4); %This is the summation part of discretized dose
% w.plus=sum(w_plus,4);
% w.minus=sum(w_minus,4);
% 
% 
% 
% for j=1:a(1)        %Make binary live/dead map for Iso
%     for jj=1:a(2)
%         for jjj=1:a(3)
%             if w.mean(j,jj,jjj) >= 1
%                 Iso.mean(j,jj,jjj)=1;
%             end
%             if w.plus(j,jj,jjj) >=1
%                 Iso.plus(j,jj,jjj)=1;
%             end
%             if w.minus(j,jj,jjj) >=1
%                 Iso.minus(j,jj,jjj) =1;
%             end
%         end
%     end
% end
% 
% end