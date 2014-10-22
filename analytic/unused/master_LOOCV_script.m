close all
clear

opttype = 'bestfit50';          % Name the opttype
datafilename = 'datasummaryL2_10sourceNewton50.txt';  % Name the datasummary file
DSC_thresholds = [0.3 0.5 0.7];
mu_thresholds = [ 150 200];  % Intervals: 0 to 100, 100 to 150, 150 to 200, 200
% and above
% DSC_thresholds = [0.3 0.7];   % At least one DSC threshold is required
% mu_thresholds = [175];   % Zero or more mu_threholds is possible. 
%Matlab_flag = 1;         % Matlab_flag must be 0 or 1. 0 runs FEM, 1 runs MATLAB Greens Function kernel

fin  = fopen(strcat('./global.ini'));

while ~feof(fin)
    f_line = fgetl(fin);
    f_line_length = length(f_line);
    if length(f_line) >= 14
        
        TF = strcmp( 'MatlabDriver =', f_line(1:14));
        
        if TF == 1
            match_or_not(1) = sum(f_line(end-3:end)=='True');
            match_or_not(2) = sum(f_line(end-4:end)=='False');
            if match_or_not(1)==4
                Matlab_flag = 1;
            elseif match_or_not(2)==5
                Matlab_flag = 0;
            else
                disp('Invalid Matlab_flag. Open ./global.ini and set "MatlabDriver" to "True" or "False"')
            end
            
            break
        end
    end
    
end



% if Matlab_flag ~=1 && Matlab_flag ~=0
%     disp('Invalid Matlab_flag. Only 0 or 1 allowed')
%     break
% end

[opt, LOOCV, fig_labels] = master_LOOCV ( opttype, datafilename, DSC_thresholds, mu_thresholds, Matlab_flag);

survival_plot (LOOCV.dice.values, opt.dice.all, fig_labels);