% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [total, dice, hd, mutual_threshold, false_pix] = metric_calculator ( MRTI_crop, MRTI_isotherm, MRTI_list,n_MRTI,model_crop,inputdatavars,summary );


% isotherms = 51:65
n_length = size(model_crop,3);
dice = zeros(n_length,15);
hd = dice;
mutual_threshold = dice;
false_pix = zeros(n_length,15,3);
total = zeros(n_length,7);
if summary.choice == 1            % mu
    total(:,1) = summary.mu;
    
elseif summary.choice ==2           % w_perf
    total(:,1) = summary.w_perf;
    
elseif summary.choice ==3           % k_cond
    total(:,1) = summary.k_cond;
    
elseif summary.choice ==4           % w_perf + mu
    
    
end

MRTI_crop = MRTI_crop';
MRTI_isotherm = permute(MRTI_isotherm, [2 1 3]);

base_level= size(MRTI_crop,1).*size(MRTI_crop,2).*37;
for ii = 1:n_length
    
    % Dice
    for kk=1:15
        model_deg_threshold = model_crop(:,:,ii) >= ( 50 + kk);
        [mod_row, mod_column] = find( model_deg_threshold ==1);
        mod_list = [(mod_row .* inputdatavars.spacing(1)) (mod_column .* inputdatavars.spacing(2))];

        n_model = sum( sum( model_deg_threshold ));
        intersection = model_deg_threshold + MRTI_isotherm(:,:,kk);
        intersection = intersection > 1;
        n_intersection = sum( sum( intersection ));
        dice(ii,kk) = 2 * n_intersection / ( n_model + n_MRTI(kk) );  % DSC
        mutual_threshold(ii,kk) = mi( model_deg_threshold, MRTI_isotherm(:,:,kk)); % MI for label map
        false_pix (ii,kk,1) = n_MRTI(kk) - n_intersection;  % False negative
        false_pix (ii,kk,2) = n_model - n_intersection; % False positive
        false_pix(ii,kk,3) = false_pix(ii,kk,1)+false_pix(ii,kk,2); % total false pixels
        if isempty(mod_list)==1
            hd(ii,kk) = 1;
        else
            hd(ii,kk) = HausdorffDist( mod_list, MRTI_list{kk});
        end
    end
    
    % L2
    temperature_diff = model_crop(:,:,ii) - MRTI_crop;
    total(ii,2) = norm ( temperature_diff , 1 ) ; % L1 norm
    
    total(ii,3) = norm ( temperature_diff , 2) ; % L2 norm
    
    total(ii,4) = norm ( temperature_diff , inf); %L_inf norm
    
    total(ii,5) = mi(model_crop(:,:,ii), MRTI_crop); % MI for temperature
    
    % model-base
    total(ii,6) = sum(sum( model_crop(:,:,ii) )) - base_level;
    
    % max temp
    total(ii,7) = max(max( model_crop(:,:,ii) ));
    
    
end

end