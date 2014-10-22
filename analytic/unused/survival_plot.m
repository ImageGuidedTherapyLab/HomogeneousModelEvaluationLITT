% 
function [] = survival_plot22 (dice_values, dice_opt, fig_labels,run1,run2);

close all

sizes.mu = size(dice_values,1);
sizes.dice = size(dice_values,2);

thresholds = linspace ( 0.0, 1, 10001);
passes = cell(sizes.mu,sizes.dice);

for ii = 1:sizes.mu
    
    for jj = 1:sizes.dice
        passes{ii,jj} = zeros(10001,1);
        
        if run1(ii,jj) >= 2
            
            for kk = 1:10001
                
                passes{ii,jj}(kk) = sum( dice_values{ii,jj} > thresholds(kk));
                
            end
        end
    end
end
clear ii jj kk

passes_opt = zeros(10001,1);
for kk = 1:10001
    
    passes_opt (kk) = sum( dice_opt.values > thresholds(kk) );
    
end
clear kk

figure(1); plot (thresholds, passes_opt, 'LineWidth',5);
title('DSC performance for optimization');
kk = 2;

for ii = 1:sizes.mu
    
    %if sum( run2(ii,:)) < sizes.dice
        figure(kk);
        
        h_title = title( strcat( fig_labels.mu_groups(ii), [' and DSC thresholds of ']));
        hold on
        for jj = 1:sizes.dice
            
            if run1(ii,jj) >=2
                hold all
                origtitle = get(h_title,'String');
                set(h_title, 'String', strcat(origtitle, [', DSC >'], num2str(fig_labels.DSC(jj))));
                
                if run1(ii,jj) == 2
                    plot(thresholds, passes{ii,jj},'LineWidth',5);
                end

                
                if jj == 1
                    legend_string = strcat( 'DSC >',num2str(fig_labels.DSC(jj)));
                    h_legend = legend( legend_string, 'Location','southwest');

                else
                    origlegend = get(h_legend,'String');
                    %set(h_legend, 'String', strcat(origlegend, [' ,  
                end
                
            end
        end
        
        hold off
        kk=kk+1;

    %end

end

end
