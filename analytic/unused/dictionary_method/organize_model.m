% 1st, it matches the power.
% 2nd, it rotates.
% 3rd, it matches the meshgrid.

function [model_crop] = organize_model ( inputdatavars, all_opt_fig, no_pwr_fig, sim_dim, choice );
% Record the working directory
setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );
n_length = size(all_opt_fig,3);

pwr_hist = str2num(inputdatavars.powerhistory);
for ii = 1:(length(pwr_hist) - 1)  % Write the for loop to iterate through all but the last index of 'pwr_hsitory'
    diff = pwr_hist(ii) - pwr_hist( ii + 1); % Record the difference between neighboring array elements
    if diff >= 0 % If 'diff' is non-negative, it means the times have stopped and the powers are starting.
        break;  % Stop the for-loop.      
    end
end
% Based on immediately previous for-loop, parse the times from powers.
power_log = pwr_hist ( (ii+1):end);
power_log = max(power_log); % Find the maximum power value

clear diff

% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.

% Make the VOI; Note the 'inputdatavars.voi' is from ParaView.
inputdatavars.voi = double( inputdatavars.voi );
VOI.x = inputdatavars.voi(3:4); % The weird index assignment is coz it's from ParaView.
VOI.y = inputdatavars.voi(1:2);
VOI.z = inputdatavars.voi(5:6);

if sum( inputdatavars.UID == '0496' ) ==4
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0402' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 3;
elseif sum( inputdatavars.UID == '0389' ) ==4
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 2;
elseif sum( inputdatavars.UID == '0385' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y - 1;
    inputdatavars.maxheatid = 109;
elseif sum( inputdatavars.UID == '0476' ) ==4
    VOI.x = VOI.x + 5;
    VOI.y = VOI.y - 5;
elseif sum( inputdatavars.UID == '0477' ) ==4
    VOI.x = VOI.x + 3;
    VOI.y = VOI.y - 2;
elseif sum( inputdatavars.UID == '0438' ) ==4
    VOI.x = VOI.x - 1;
    VOI.y = VOI.y + 0;
elseif sum( inputdatavars.UID == '0435' ) ==4
    VOI.x = VOI.x - 3;
    VOI.y = VOI.y + 3;
elseif sum( inputdatavars.UID == '0436' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 6;
    power_log = 10.05;
    inputdatavars.maxheatid = 39;
elseif sum( inputdatavars.UID == '0466' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
    power_log = 12;
elseif sum( inputdatavars.UID == '0468' ) ==4
    VOI.x = VOI.x - 3;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0471' ) ==4
    VOI.x = VOI.x + 3;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0447' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0457' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 0;
    power_log = 10.2;
elseif sum( inputdatavars.UID == '0453' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0451' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0418' ) ==4
    VOI.x = VOI.x + 0;
    VOI.y = VOI.y + 0;
    %inputdatavars.maxheatid = 85;
elseif sum( inputdatavars.UID == '0409' ) ==4
    VOI.x = VOI.x - 2;
    VOI.y = VOI.y + 11;
elseif sum( inputdatavars.UID == '0414' ) ==4
    VOI.x = VOI.x + 2;
    VOI.y = VOI.y + 1;
elseif sum( inputdatavars.UID == '0415' ) ==4
    VOI.x = VOI.x + 1;
    VOI.y = VOI.y + 1;
end

% Make the correct power
if choice == 1

    no_pwr = repmat(no_pwr_fig,[1 1 n_length]);
    model_temp = all_opt_fig .* power_log - (power_log-1).*no_pwr;
    
elseif choice ==4
    
    w_Num = size( no_pwr_fig,3);
    mu_Num = size( all_opt_fig,3) / w_Num;
    
    Temp_sz = size(all_opt_fig);
    no_pwr = zeros( Temp_sz(1), Temp_sz(2), Temp_sz(3) );
    
    for ii = 1:w_Num
        
        ix = ii : w_Num : Temp_sz(3);
        no_pwr(:,:,ix) = repmat( no_pwr_fig(:,:,ii), [1 1 mu_Num]);
    end
    clear ii
    model_temp = all_opt_fig .* power_log - (power_log - 1).*no_pwr;

else
    model_temp = all_opt_fig .* power_log - (power_log-1).*no_pwr_fig;
end

% Rotate and add body temp
model_temp = imrotate(model_temp, str2num( inputdatavars.cv.z_rotate),'crop' );
model_temp = model_temp + 37; % Get to body temp 

% Set up MRTI meshgrid
FOV_l(1) = inputdatavars.spacing(1) .* (inputdatavars.voi(2)-inputdatavars.voi(1)+1);
FOV_l(2) = inputdatavars.spacing(2) .* (inputdatavars.voi(4)-inputdatavars.voi(3)+1);
%[MRXq, MRYq] = meshgrid ( inputdatavars.spacing(1): inputdatavars.spacing(1) : FOV_l(1) , inputdatavars.spacing(2): inputdatavars.spacing(2) : FOV_l(2) );
[MRXq, MRYq] = meshgrid ( -(FOV_l(1)-inputdatavars.spacing(1))/2: inputdatavars.spacing(1) : (FOV_l(1)-inputdatavars.spacing(1))/2 , -(FOV_l(2)-inputdatavars.spacing(2))/2: inputdatavars.spacing(2) : (FOV_l(2)-inputdatavars.spacing(2))/2 );

% Set up model meshgrid
modFOV_l(1) = sim_dim.spacing.x .*sim_dim.mod_point.x .*sim_dim.mod_point.z_subslice;
modFOV_l(2) = sim_dim.spacing.y .*sim_dim.mod_point.y .*sim_dim.mod_point.z_subslice;
[modX, modY] = meshgrid ( -(modFOV_l(1)-sim_dim.spacing.x)/2: sim_dim.spacing.x : (modFOV_l(1)-sim_dim.spacing.x)/2 , -(modFOV_l(2)-sim_dim.spacing.y)/2: sim_dim.spacing.y : (modFOV_l(2)-sim_dim.spacing.y)/2 );

% Do interpolation
n_model = size( model_temp,3);
model_crop = zeros( (inputdatavars.voi(4)-inputdatavars.voi(3)+1), (inputdatavars.voi(2)-inputdatavars.voi(1)+1), n_model);
for ii = 1:n_model
    model_crop(:,:,ii) = qinterp2( modX, modY, model_temp(:,:,ii), MRXq, MRYq);
end
clear ii
% model_isotherm = repmat( model_crop, [1 1 1 15]);
% 
% for ii = 1:15
%     model_isotherm( model_isotherm(:,:,:,ii)) = model_isotherm(:,:,:,ii) >= 50 +ii;
%     
% end


end