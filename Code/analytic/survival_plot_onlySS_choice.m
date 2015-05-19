%
function [] = survival_plot_onlySS_choice (dice_values, naive_values, run_var_eff_val, dice_opt, fig_labels,run1,run2,naive_tag,choice);

%close all

sizes.var = size(dice_values,1);
sizes.dice = size(dice_values,2);

thresholds = linspace ( 0.0, 1, 10001);
passes = cell(sizes.var,sizes.dice);
naive_pass = passes;
L2_mean_pass =passes;
L2_median_pass = passes;
dice_mean_pass = passes;
dice_median_pass = passes;
for ii = 1:sizes.var
    
    for jj = 1:sizes.dice
        passes{ii,jj} = zeros(10001,1);
        naive_pass{ii,jj} = passes{ii,jj};
%         L2_mean_pass{ii,jj} = passes{ii,jj};
%         L2_median_pass{ii,jj} = passes{ii,jj};
%         dice_mean_pass{ii,jj} = passes{ii,jj};
%         dice_median_pass{ii,jj} = passes{ii,jj};
        
        if run1(ii,jj) >= 2
            
            for kk = 1:10001
                
                passes{ii,jj}(kk) = sum( dice_values{ii,jj} > thresholds(kk));
                naive_pass{ii,jj}(kk) = sum( naive_values{ii,jj} > thresholds(kk));
%                 L2_mean_pass{ii,jj}(kk) = sum( best.group.L2_mean.val{ii,jj} > thresholds(kk));
%                 L2_median_pass{ii,jj}(kk) = sum( best.group.L2_median.val{ii,jj} > thresholds(kk));
%                 dice_mean_pass{ii,jj}(kk) = sum( best.group.dice_mean.val{ii,jj} > thresholds(kk));
%                 dice_median_pass{ii,jj}(kk) = sum( best.group.dice_median.val{ii,jj} > thresholds(kk));
                
            end
        end
    end
end
clear ii jj kk

passes_opt = zeros(10001,1);
pass_naive_opt = passes_opt;
% all_L2_mean_pass = passes_opt;
% all_L2_median_pass = passes_opt;
% all_dice_mean_pass = passes_opt;
% all_dice_median_pass = passes_opt;
for kk = 1:10001
    
    passes_opt (kk) = sum( dice_opt.values > thresholds(kk) );
    pass_naive_opt (kk) = sum ( naive_values{1} > thresholds(kk) );
%     all_L2_mean_pass (kk) = sum( best.all.L2_mean.val > thresholds(kk) );
%     all_L2_median_pass (kk) = sum( best.all.L2_median.val > thresholds(kk) );
%     all_dice_mean_pass (kk) = sum( best.all.dice_mean.val > thresholds(kk) );
%     all_dice_median_pass (kk) = sum( best.all.dice_median.val > thresholds(kk) );
end
clear kk

% figure; h_title = title( 'DSC performance for optimization versus literature and single best guesses'); hold all;
% [h1] = plot (thresholds, [passes_opt pass_naive_opt all_L2_mean_pass all_L2_median_pass all_dice_mean_pass all_dice_median_pass]);
% legend( h1, 'Optimization', strcat(['Literature '], num2str(naive_tag(1))), 'L_2 mean','L_2 median','dice mean','dice median') ;
% legend('-DynamicLegend', 'Location','southwest');hold off;

figure; h_title = title( 'DSC performance for optimization versus literature'); hold all;
[h1] = plot (thresholds, [passes_opt pass_naive_opt],'LineWidth',5);
%legend( h1, 'Optimization', ['Literature ', num2str(naive_tag(1)),' m^{-1}' ]) ;
%legend('-DynamicLegend', 'Location','southwest');hold off;

kk = 2;

for ii = 1:sizes.var
    

    if sum(run2(ii,:)) >=1
        
        figure;
        
        h_title = title( strcat( fig_labels.var_groups(ii), [' and DSC thresholds of ']));
        hold on
        for jj = 1:sizes.dice
            
            if run1(ii,jj) >=2
                hold all
                origtitle = get(h_title,'String');
                set(h_title, 'String', strcat(origtitle, [', DSC >'], num2str(fig_labels.DSC(jj))));
                set(findobj('type','axes'),'fontsize',13);
                
                if run1(ii,jj) == 2
                    
                    if jj == 1
                        plot(thresholds, passes{ii,jj},'LineWidth',5);
                        legend_string = ['DSC > ', num2str(fig_labels.DSC(jj))];
                        %legend( legend_string, 'Location','southwest');
                        hold all
                        %legend('-DynamicLegend', 'Location','southwest');
                        
                        if naive_tag(2) ==1
                            %legend_string = ['Naive ', num2str(naive_tag(1)), ' m^{-1}'];
                            plot (thresholds, naive_pass{ii,jj}, 'DisplayName',legend_string, 'LineWidth',5);
                        end
                        
                    else
                        %legend_string = strcat( 'DSC >',num2str(fig_labels.DSC(jj)));
                        set(findobj('type','axes'),'fontsize',13);
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


for ii = 1:sizes.var
    

    if sum(run2(ii,:)) >=1
        
        figure;
        
        h_title = title( strcat( fig_labels.var_groups(ii), [' and DSC thresholds of ']));
        hold on
        for jj = 1:sizes.dice
            
            if run1(ii,jj) >=2
                hold all
                origtitle = get(h_title,'String');
                set(h_title, 'String', strcat(origtitle, [', DSC >'], num2str(fig_labels.DSC(jj)),', histogram of run var values'));
                
                if run1(ii,jj) == 2
                    
                    if jj == 1
                        hist( run_var_eff_val{ii,jj} );

                    end
                end
                
                
            end
        end
        
        hold off
        kk=kk+1;
        
    end
    
end


end
