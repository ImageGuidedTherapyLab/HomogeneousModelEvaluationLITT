%This script is meant to be run in SI units.
%As a function of power, P, this makes a prediction of the temperature
%within the domain, dom. 'dom' is a structure representing the simulated domain with
%three spacial dimensions in meters and the number of points in each
%dimension (a total of 6 fields). The number of points in 'dom' should be
%odd in order for the centroid to be on the origin and a point to be on the
%origin.

%'source' is a structure with the details of the source. The 'n' field is
%the number of sources. The 'length' is the diffusing tip length, commonly
%0.01 m.

function [tmap]=Bioheat1D_SSS_reza_comparison (P,dom,source,w,k,g,mua,mus,probe_u,robin_co,cblood);

%List of space and time details
r1=0.0015/2; % (m) R1 is the distance from the isotropic laser source point and the edge of the fiber
r2=1; % (m) R2 is the maximum edge of the domain;
%Npowers=size(P,1)+1; %Returns how many timesteps will be calculated
% P(2:Npowers,:)=P(1:(Npowers-1),:); %Inserts the P=0 timestep
% P(1,1)=1;
% P(1,2)=0;
P = 10;
P=P/source.n;  %Scales the power to the number of source points

%List of constants
% mua=500; % 1/m
% mus=14000; %1/m
% g=0.88; % Unity
% k=0.527; % W/(m * K)
% w=6; % kg / (m^3 * s)
u0 = probe_u - 37; % K
ua = 0; % K

%
%Points structure
points.x=linspace(-0.015, 0.015, 101);
points.y=points.x;
%points.z=linspace(-dom.z/2,dom.z/2,dom.pointz);
points.z=linspace(-0.015,0.015,7);

dom.pointx = length(points.x);
dom.pointy = length(points.y);
dom.pointz = length(points.z);

%Initialize t_sample and r

t_sample=zeros(dom.pointx,dom.pointy,dom.pointz,source.n);
r=zeros(source.n,1);
%r = zeros(dom.pointx,dom.pointy,dom.pointz,source.n);
%Spatial locations of the sources; My convention is that the long axis of
%the laser is parallel to the y-axis.
%laser=linspace((-source.length/2),(source.length/2),source.n);

%Giant for loop vector for each source, calculate the t_sample; i, ii, iii
%are spatial; j and jj are non-spatial
mutr  = mua+mus*(1-g);
mueff = sqrt(3.0*mua*mutr);

for i=1:dom.pointx   %Spatial loop for i, ii, iii
    
    for ii=1:dom.pointy
        
        for iii=1:dom.pointz % This takes advantage of symmetry
            
            for j=1:source.n    %Loop for the separate isotropic sources
                %r(i,ii,iii,j)=sqrt( (source.laser(j,1)-points.x(i))^2 + (source.laser(j,2)-points.y(ii))^2+(source.laser(j,3)-points.z(iii))^2);
                r(j)=sqrt( (source.laser(j,1)-points.x(i))^2 + (source.laser(j,2)-points.y(ii))^2+(source.laser(j,3)-points.z(iii))^2);  %Distance for each isotropic source
                if (r(j) <= r1) && (source.n ==1)
                    t_sample(i,ii,iii,j) = u0;
                elseif (r(j) <= r1) && (source.n == 10)
                    t_sample(i,ii,iii,j) = u0+55;   % This is tuned to the number of source points... i think
                elseif (r(j) <= r1) && (source.n > 1)
                    t_sample(i,ii,iii,j) = u0;
                else
                    %[t_sample(i,ii,iii,j)]=sammodel1D(u0,ua,k,w,P,r(j),mua,mus,R1,R2,g); % Calculate the temperature from one source
                    %t_sample(i,ii,iii,j)=mathematicaSSS(u0,ua,k,w,P,r(j),mua,mus,r1,r2,g,cblood);
                    %t_sample(i,ii,iii,j)=change_AnalyticSSS_mod3(u0,ua,k,w,P,r(j),mua,mus,r1,r2,g,cblood);
                    %t_sample=mathematicaSSS3(r(j),mua,mus,g);
                    t_sample(i,ii,iii,j)=mathematicaSSS6(u0,ua,k,w,P,r(j),mueff,r1,r2,cblood);
                end
            end
        end
    end
end

clear i ii iii j

aa = size(t_sample);

if length(aa) == 3
%     t_summing = zeros(aa(1),aa(2),aa(3));
%     t_summing(:,:,1) = t_sample(:,:,1);
%     for ii = 2:aa(3)
%         t_summing(:,:,ii) = 2.*t_sample(:,:,ii);
%     end
%     tmap = (1/dom.z_subslice).*sum(t_summing,3);
else
%     t_summing = zeros(aa(1),aa(2),aa(3),aa(4));
%     t_summing(:,:,1,:) = t_sample(:,:,1,:);
%     for ii = 2:aa(3)
%         t_summing(:,:,ii,:) = 2.*t_sample(:,:,ii,:);
%     end
%     t_avg = squeeze((1/dom.z_subslice).*sum(t_summing,3));
%     tmap = sum(t_avg,3);
tmap = squeeze( sum(t_sample,4) );
end
    
% if j == 1
%     t_avg = (1/dom.z_subslice).*(t_sample(:,:,1)+2.*t_sample(:,:,2)+2.*t_sample(:,:,3));
% end
% tmap = sum(t_sample,4);  %sum(t_sample(j,jj,1));

end