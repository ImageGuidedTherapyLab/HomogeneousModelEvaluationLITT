% 
function [] = survival_plot (dice_values, dice_opt);



sizes.mu = size(dice_values,1);
sizes.dice = size(dice_values,2);
is_empty_matrix = zeros(sizes.mu, sizes.dice);

for ii = 1:sizes.mu
    for jj = 1:sizes.dice
        is_empty_matrix(ii,jj) = isempty(dice_values{ii,jj});
    end
end
clear ii jj

thresholds = linspace ( 0.0, 1, 10001);
passes = cell(sizes.mu,sizes.dice);

for ii = 1:sizes.mu
    
    for jj = 1:sizes.dice
        passes{ii,jj} = zeros(10001,1);
        
        if is_empty_matrix(ii,jj) == 0
            
            for kk = 1:10001
                
                passes{ii,jj}(kk) = sum( dice_values{ii,jj} > thresholds(kk));
                
            end
        end
    end
end
clear ii jj kk

passes_opt = zeros(10001,1);
for kk = 1:10001
    
    passes_opt (kk) = sum( dice_opt > thresholds(kk) );
    
end
clear kk

figure(1); plot (thresholds, passes_opt, 'LineWidth',5);
kk = 2;
% ll = 1;
for ii = 1:sizes.mu
    
    if sum(is_empty_matrix(ii,:)) < sizes.dice
        figure(kk);
        
        for jj = 1:sizes.dice
            
            if is_empty_matrix(ii,jj) ==0
                
                plot(thresholds, passes{ii,jj},'LineWidth',5);
                hold all
                
            end
        end
        
        hold off
        kk=kk+1;

    end
%     ll=ll+1;
end

end
% figure(1); plot (thresholds, passes_LOOCV8,'LineWidth',5);
% figure; plot (thresholds, passes_LOOCV7);
% figure; plot (thresholds, passes_LOOCV25,'LineWidth',5);
% figure; plot (thresholds, passes_LOOCV3,'LineWidth',5);