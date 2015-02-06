clear
close all
clc

setenv ( 'PATH22' , pwd);
path22 = getenv ( 'PATH22' );

cd ../../../MATLAB/Tests/in_silico/
load ('in_silico_test.mat');
cd (path22);

mu_length = size(total,1);
n_length = size(total{1,3},1);
xx = log10 (total{1,3}(2:end,1));
sum_false = total{1,4}(2:end,7) + total{1,5}(2:end,7);

figure;
h_title = title( 'Convergence of maximum point-wise difference with number of point sources');
hold all;
xlabel('Log_{10} of sources (Unity)');
ylabel('Maximum of point-wise temperature difference (^{o}C)');
[h1] = plot( xx, total{1,3}(2:end,3), '--*');
hold off;

figure;
h_title = title( 'Convergence of DSC with number of point sources');
hold all;
xlabel('Log_{10} of sources (Unity)');
ylabel('DSC (Unity)');
[h1] = plot( xx, total{1,6}(2:end,7), '--*');
hold off;

figure;
h_title = title( 'Reduction of total false pixels with number of point sources');
hold all;
xlabel('Log_{10} of sources (Unity)');
ylabel('FP + FN (Unity)');
[h1] = plot( xx, sum_false, '--*');
hold off;

%ylabel (ax(1), ....)
%ylabel (ax(2), ....)

% legend( h1, 'Optimization', ['Literature ', num2str(naive_tag(1)),' m^{-1}' ]) ;
% legend('-DynamicLegend', 'Location','southwest');hold off;
% 
% 
% plot( xx, total{1,3}(2:end,3));

% 
% for ii = 1:mu_length
%     figure;
%     for jj = 1:n_length
%         
%         plot( xx, 