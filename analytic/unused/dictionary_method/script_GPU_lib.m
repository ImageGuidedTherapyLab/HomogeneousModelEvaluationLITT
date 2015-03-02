choice = 1;

GPU_generate_lib(choice);

metric_generation(choice);

%MC_2dimension_GPU (choice);

% for ii = 2:3
%     GPU_generate_lib(ii);
%     
%     metric_generation(ii);
% end