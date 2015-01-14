%This function accepts the unique tmaps and returns the tmaps that
%correspond to MRTI.

function [tmap]=Build_tmap_history(tmap_unique,delta_P);

tmap=zeros(size(tmap_unique,1),size(tmap_unique,2),size(tmap_unique,3),size(delta_P,1));

j=1;
for i=1:length(delta_P)

    if delta_P(i) ~= 0
        j=j+1;
    end
    
    tmap(:,:,:,i)=tmap_unique(:,:,:,j);
    
end
    
end