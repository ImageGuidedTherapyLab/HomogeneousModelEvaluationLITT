%This parses the power further so that there are no repeats, redundant
%intervals.

function [Power_interval_str]=power_parser_write_DF_array(power_log);

%Find where the power changes at all;
delta_P = (diff( power_log(:,5) )~=0) + ( diff(power_log(:,6) )~=0);  %Find which elements change from columns 5 and 6; then add the changes into one column
delta_P(1+1:end+1,:) = delta_P (1:end,:);  %add row back in because diff function eliminates the first one.
delta_P(1,:) = delta_P (2,:);
%At this point, delta_P lists all the times that columns 5 and 6 change

% Check to see if there is anything changing
if isempty ( find ( delta_P )) == 1;
    
    Power_intervals (1,1) = power_log ( end, 4 );
    Power_intervals (1,1) = floor ( Power_intervals (1,1) / 2 ); % cut the number of times in two
    Power_intervals (1,2) = power_log ( end, 6 );
    
    
    str_power = [ '0' '0' ];
    
    Power_interval_str = '[[';
    Power_interval_str = strcat( Power_interval_str, str_power( 1, 1) );
    Power_interval_str = strcat( Power_interval_str, '],[' );
    Power_interval_str = strcat( Power_interval_str, str_power( 1, 2) );
    Power_interval_str = strcat( Power_interval_str, ']]' );
    
    return
end

keep = find ( delta_P) ; % Collects the indices of where the power changes
on_off = zeros ( length( find ( delta_P ) )  , 1);

%This loop captures the on/off state of the power_log
for ii = 1 : length ( keep )
    
    on_off (ii) = power_log( keep(ii) , 5 );
    
end

clear ii

P = find (delta_P); %Column 1 of P records the times that the power changes

P(:,2) = power_log (P(:,1),6); %Use the times from column 1, P to record the corresponding powers


P(:,2) = P(:,2) * 15/100; %Convert % power to W power
P ( :,3 ) = on_off; %Add a column of on/off for next calculation

k_P (:,1)=P(:,1);
k_P (:,2)=P(:,2).*P(:,3);

times = k_P( : , 1 ); %Split in two for easy indexing.
powers = k_P( :,2);

% Add a 0 power at the beginning.
powers = cat (1,0,powers);
% Add the ending time to the time series and divide by 2 to round
times(end+1)=power_log(end,4);
times = floor (times/2);

Power_intervals (:,1) = times;  %Write the Power_intervals
Power_intervals (:,2) = powers;

%Cut repitious power/time information where there isn't a change.
Power_intervals_diff = diff (Power_intervals);
cut_indices = find ( not ( Power_intervals_diff ( :,2) ));
Power_intervals ( cut_indices , : )= [];

%Make the Power_intervals a cell array of strings
str_power = cellstr( num2str( Power_intervals )); % Make the cell array
str_power = strtrim( str_power ); % Cut the left and right white space
str_power = regexp( str_power , '\s+' , 'split'); % Use the remaining middle white space to make two columns.

Power_interval_str = '[[';

for ii = 1 : ( size( str_power, 1 ) - 1)
    
    Power_interval_str = strcat( Power_interval_str, str_power{ ii }{1} ); % Write the number for time
    Power_interval_str = strcat( Power_interval_str, ', ' ); % Write the delimeter ', '
    
end

Power_interval_str = strcat( Power_interval_str, str_power{ ii + 1 }{1} ); % Write the last time
Power_interval_str = strcat( Power_interval_str , '],[' ); % Write the middle delimeter between time and power

for ii = 1 : ( size( str_power, 1 ) - 1 )
    
    Power_interval_str = strcat( Power_interval_str, str_power{ ii }{2} ); % Write the number for power
    Power_interval_str = strcat( Power_interval_str, ', ' ); % Write the delimeter ', '
    
end

Power_interval_str = strcat( Power_interval_str, str_power{ ii + 1 }{2} ); % Write the last power
Power_interval_str = strcat( Power_interval_str, ']]' ); % Write the ending delimeter
Power_interval_str = strrep( Power_interval_str, ',', ', ');
Power_interval_str = strrep( Power_interval_str, '], [' , '],[');

end