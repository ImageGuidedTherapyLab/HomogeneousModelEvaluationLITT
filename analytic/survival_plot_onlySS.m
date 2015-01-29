%
function [] = survival_plot_onlySS (dice_values, naive_values, dice_opt, fig_labels,run1,run2,naive_tag);

%close all

sizes.mu = size(dice_values,1);
sizes.dice = size(dice_values,2);

thresholds = linspace ( 0.0, 1, 10001);
passes = cell(sizes.mu,sizes.dice);
naive_pass = passes;
for ii = 1:sizes.mu
    
    for jj = 1:sizes.dice
        passes{ii,jj} = zeros(10001,1);
        
        if run1(ii,jj) >= 2
            
            for kk = 1:10001
                
                passes{ii,jj}(kk) = sum( dice_values{ii,jj} > thresholds(kk));
                naive_pass{ii,jj}(kk) = sum( naive_values{ii,jj} > thresholds(kk));
                
            end
        end
    end
end
clear ii jj kk

passes_opt = zeros(10001,1);
pass_naive_opt = passes_opt;
for kk = 1:10001
    
    passes_opt (kk) = sum( dice_opt.values > thresholds(kk) );
    pass_naive_opt (kk) = sum ( dice_opt.naive.val > thresholds(kk) );
end
clear kk

figure(1); plot (thresholds, passes_opt, 'LineWidth',5);
title('DSC performance for optimization compared to naive guess');
legend_string = strcat( ['Optimization']);
legend( legend_string, 'Location','southwest');
hold all
legend('-DynamicLegend', 'Location','southwest');
legend_string = strcat(['Naive '], num2str(naive_tag(1)), ' m^-1');
plot (thresholds, pass_naive_opt, 'DisplayName',legend_string, 'LineWidth',5);
hold off

kk = 2;

for ii = 1:sizes.mu
    

    if sum(run2(ii,:)) >=1
        
        figure(kk);
        
        h_title = title( strcat( fig_labels.mu_groups(ii), [' and DSC thresholds of ']));
        hold on
        for jj = 1:sizes.dice
            
            if run1(ii,jj) >=2
                hold all
                origtitle = get(h_title,'String');
                set(h_title, 'String', strcat(origtitle, [', DSC >'], num2str(fig_labels.DSC(jj))));
                
                if run1(ii,jj) == 2
                    
                    if jj == 1
                        plot(thresholds, passes{ii,jj},'LineWidth',5);
                        legend_string = strcat( 'DSC >',num2str(fig_labels.DSC(jj)));
                        legend( legend_string, 'Location','southwest');
                        hold all
                        legend('-DynamicLegend', 'Location','southwest');
                        
                        if naive_tag(2) ==1
                            legend_string = strcat(['Naive '], num2str(naive_tag(1)), ' m^-1');
                            plot (thresholds, naive_pass{ii,jj}, 'DisplayName',legend_string, 'LineWidth',5);
                        end
                        
                    else
                        legend_string = strcat( 'DSC >',num2str(fig_labels.DSC(jj)));
                        plot(thresholds, passes{ii,jj},'DisplayName',legend_string,'LineWidth',5);

                        %legend('-DynamicLegend', 'Location','southwest');
                        %legend(legend_string, 'Location','southwest');
                        %legend(legend_string,'-DynamicLegend', 'Location','southwest');
                        %legend( legend_string, 'Location','southwest');
                        
                        %set(h_legend, 'String', strcat(origlegend, [' ,
                    end
                end
                
                
            end
        end
        
        hold off
        kk=kk+1;
        
    end
    
end

end
