% This is the updated Bioheat_script that should be used with DF's DAKOTA
% run. The metric is based on temperature (not dose and isotherms).

function [total, dice, hd, mutual_threshold, false_pix] = metric_calculator_rand ( MRTI_crop, MRTI_isotherm, MRTI_list,n_MRTI,model_crop,inputdatavars,summary );


% isotherms = 51:65
n_length = size(model_crop,3);
dice = zeros(n_length,15);
hd = dice;
mutual_threshold = dice;
false_pix = zeros(n_length,15,3);
total = zeros(n_length,8);
total(:,1) = summary.mu;
total(:,2) = summary.w_perf;

MRTI_crop = MRTI_crop';
MRTI_isotherm = permute(MRTI_isotherm, [2 1 3]);

base_level= size(MRTI_crop,1).*size(MRTI_crop,2).*37;

% n_model= zeros( n_length,1);
% n_intersection = n_model;
% 
% for kk = 1:15
%     model_deg_threshold = model_crop;
%     thresh = 50 +kk;
%     model_deg_threshold (model_deg_threshold < thresh) = 0;
%     model_deg_threshold (model_deg_threshold > thresh) = 1;
% 
%     MRTI_isotherm_rep = repmat(MRTI_isotherm(:,:,kk) ,[1 1 n_length]);
% 
%     intersection = model_deg_threshold + MRTI_isotherm_rep;
%     intersection( intersection > 1)=1;
% 
%     for ii = 1:n_length
%         [mod_row, mod_column] = find( model_deg_threshold(:,:,ii) ==1);
%         mod_list = [(mod_row .* inputdatavars.spacing(1)) (mod_column .* inputdatavars.spacing(2))];
% 
%         n_model(ii) = sum( sum( model_deg_threshold(:,:,ii)));
%         n_intersection(ii) = sum( sum( intersection(:,:,ii)));
% 
%         mutual_threshold(ii,kk) = mi( model_deg_threshold(:,:,ii), MRTI_isotherm_rep(:,:,ii)); % MI for label map
% 
% 
%     end
%     dice(:,kk) = 2 .* n_intersection ./ (n_model+ n_MRTI(kk));
% 
%     false_pix (:,kk,1) = n_MRTI(kk) - n_intersection;  % False negative
%     false_pix (:,kk,2) = n_model - n_intersection; % False positive
%     false_pix(:,kk,3) = false_pix(:,kk,1)+false_pix(:,kk,2); % total false pixels
% 
%     for ii = 1:n_length
%         temperature_diff = model_crop(:,:,ii) - MRTI_crop;
%         total(ii,3) = norm ( temperature_diff , 1 ) ; % L1 norm
% 
%         total(ii,4) = norm ( temperature_diff , 2) ; % L2 norm
% 
%         total(ii,5) = norm ( temperature_diff , inf); %L_inf norm
% 
%         total(ii,6) = mi(model_crop(:,:,ii), MRTI_crop); % MI for temperature
% 
%         % model-base
%         total(ii,7) = sum(sum( model_crop(:,:,ii) )) - base_level;
% 
%         % max temp
%         total(ii,8) = max(max( model_crop(:,:,ii) ));
%     end
% 
% end



for ii = 1:n_length
    
    % Dice
    for kk=1:15
        model_deg_threshold = model_crop(:,:,ii) >= ( 50 + kk);
        [mod_row, mod_column] = find( model_deg_threshold ==1);
        mod_list = [(mod_row .* inputdatavars.spacing(1)) (mod_column .* inputdatavars.spacing(2))];

        n_model = sum( sum( model_deg_threshold));
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
    total(ii,3) = norm ( temperature_diff , 1 ) ; % L1 norm
    
    total(ii,4) = norm ( temperature_diff , 2) ; % L2 norm
    
    total(ii,5) = norm ( temperature_diff , inf); %L_inf norm
    
    total(ii,6) = mi(model_crop(:,:,ii), MRTI_crop); % MI for temperature
    
    % model-base
    total(ii,7) = sum(sum( model_crop(:,:,ii) )) - base_level;
    
    % max temp
    total(ii,8) = max(max( model_crop(:,:,ii) ));
    
    
end
toc
end