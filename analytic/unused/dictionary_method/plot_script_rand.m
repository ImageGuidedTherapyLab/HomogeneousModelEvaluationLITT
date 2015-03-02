w_lim = [3.25 16];
w_pix = 100;
w_array = linspace( w_lim(1),w_lim(2),w_pix);
mu_lim = [55 295];
mu_pix = 40;
mu_array = linspace( mu_lim(1),mu_lim(2),mu_pix);
%[paraXq, paraYq] = meshgrid (w_lim(1): 0.25 : w_lim(2), mu_lim(1): 50: mu_lim(2) );
[paraXq, paraYq] = meshgrid ( w_array, mu_array);

%[scatterX, scatterY] = meshgrid ( total{1,2}(:,1),total{1,2}(:,2));

[x y z]=griddata(total{1,2}(:,1),total{1,2}(:,2),total{1,2}(:,4),paraYq,paraXq);
figure; contourf(x,y,z);colorbar;

[x y z]=griddata(total{1,2}(:,1),total{1,2}(:,2),total{1,3}(:,7),paraYq,paraXq);
figure; contourf(x,y,z);colorbar;


para_grid = griddata( total{1,2}(:,1), total{1,2}(:,2), total{1,3}(:,7), paraXq, paraYq);
para_pic = interp2 ( scatterX, scatterY, total{1,3}(:,7), paraXq, paraYq);

figure; imagesc(para_pic);