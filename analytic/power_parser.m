%This parses the power further

function [P,unique_P,delta_P, DF_Intervals]=power_parser(power_log);

%Find where the power changes at all;
delta_P = (diff( power_log(:,5) )~=0) + ( diff(power_log(:,6) )~=0);  %Find which elements change from columns 5 and 6; then add the changes into one column
delta_P(1+1:end+1,:) = delta_P (1:end,:);  %add row back in because diff function eliminates the first one.
delta_P(1,:) = delta_P (2,:); 
%At this point, delta_P lists all the times that columns 5 and 6 change

delta_P = power_log (:,5) .* delta_P;  %Only keeps power changes while laser is on (power setting can change while the laser is off, and I don't care about those changes)

P = find (delta_P); %Column 1 of P records the times that the power changes

P(:,2) = power_log (P(:,1),6); %Use the times from column 1, P to record the corresponding powers


P(:,2) = P(:,2) * 15/100; %Convert % power to W power
unique_P = unique (P(:,2));

Power_intervals = zeros ( length( unique_P ) , 1 );
Power_intervals(1+1:end+1) = unique_P (1:end);  %add row back in because diff function eliminates the first one.
Power_intervals(1) = 0; 

DF_Intervals = cell ( 2 , 2 );  %Initiate cell array
DF_Intervals { 1 , 1 } = 'Time'; % Time and power headers
DF_Intervals { 1 , 2 } = 'Power';
DF_Intervals { 2 , 1 } = find (delta_P); % Time and power data
DF_Intervals { 2 , 2 } = Power_intervals;

end