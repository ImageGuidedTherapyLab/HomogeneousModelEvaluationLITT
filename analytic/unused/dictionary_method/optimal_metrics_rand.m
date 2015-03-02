function [ total ] = optimal_metrics_rand (total,ii);

% Record optimal L information
total{ii,7} = zeros(3,4);
[total{ii,7}(1,1) , index] = min (total{ii,2}(:,3));  % record optimal L1
total{ii,7}(1,2:3) = total{ii,2}(index,1:2);   % record values that produces optimal L1
total{ii,7}(1,4) = index; % record index that produces optimal L1
[total{ii,7}(2,1) , index] = min (total{ii,2}(:,4));  % record optimal L2
total{ii,7}(2,2:3) = total{ii,2}(index,1:2);   % record value that produces optimal L2
total{ii,7}(2,4) = index; % record index that produces optimal L2
[total{ii,7}(3,1) , index] = min (total{ii,2}(:,5));  % record optimal L_inf
total{ii,7}(3,2:3) = total{ii,2}(index,1:2);   % record value that produces optimal L_inf
total{ii,7}(3,4) = index; % record index that produces optimal L_inf

% Record optimal 57 C isotherm DSC information
total{ii,8} = zeros(1,4);
[total{ii,8}(1) , index] = max (total{ii,3}(:,7));  % record optimal Dice
total{ii,8}(2:3) = total{ii,2}(index,1:2);   % record value that produces optimal Dice
total{ii,8}(4) = index; % record index that produces optimal Dice

% Record optimal 57 C Hausdorff distance information
total{ii,9} = zeros(1,4);
[total{ii,9}(1) , index] = min (total{ii,4}(:,7)); % record optimal Hausdorff Distance
total{ii,9}(2:3) = total{ii,2}(index,1:2); % record value that produces optimal Hausdorff distance
total{ii,9}(4) = index;

% Record optimal temperature MI and 57 C isotherm MI
total{ii,10} = zeros(2,4);
[total{ii,10}(1,1) , index] = max (total{ii,2}(:,5));  % MI for temperature
total{ii,10}(1,2:3) = total{ii,2}(index,1:2); % record value that produces optimal MI for temperature
total{ii,10}(1,4) = index;
[total{ii,10}(2,1) , index] = max (total{ii,5}(:,7));  % MI for 57 C isotherm
total{ii,10}(2,2:3) = total{ii,2}(index,1:2); % record value that produces optimal MI for 57 C isotherm
total{ii,10}(2,4) = index;

% Record optimal number of false pixels for 57 C isotherm
total{ii,11} = zeros(1,4);
[total{ii,11}(1,1) , index] = min (total{ii,6}(:,7,3));  % False pixel number for 57 C isotherm
total{ii,11}(1,2:3) = total{ii,2}(index,1:2); % record value that produces number of false pixels
total{ii,11}(1,4) = index;

end