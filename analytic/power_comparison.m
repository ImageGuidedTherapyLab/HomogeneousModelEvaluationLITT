cd /mnt/FUS4/data2/sjfahrenholtz/MATLAB/Tests/direct_search/troubleshoot/Power

pwr0 = load('pwr0.mat');
pwr1 = load('pwr1.mat');
pwr5 = load('pwr5.mat');
pwr10= load('pwr10.mat');
pwr15= load('pwr15.mat');

neg3 = load('neg3.mat');

hi_pwr = pwr15.Temp(:,:,1)./15 + (14/15).*pwr0.Temp(:,:,1);
diff = pwr1.Temp(:,:,1)-hi_pwr;

figure; imagesc( pwr1.Temp(:,:,1));
figure; imagesc( hi_pwr );
figure; imagesc( diff);
figure; imagesc( pwr15.Temp(:,:,1));

hi_pwr = pwr10.Temp(:,:,1)./10 + (9/10).*pwr0.Temp(:,:,1);
diff = pwr1.Temp(:,:,1)-hi_pwr;

figure; imagesc( pwr1.Temp(:,:,1));
figure; imagesc( hi_pwr );
figure; imagesc( diff);
figure; imagesc( pwr10.Temp(:,:,1));

hi_pwr = pwr10.Temp(:,:,3)./10 + (9/10).*pwr0.Temp(:,:,1);
diff = pwr1.Temp(:,:,3)-hi_pwr;

figure; imagesc( pwr1.Temp(:,:,3));
figure; imagesc( hi_pwr );
figure; imagesc( diff);
figure; imagesc( pwr10.Temp(:,:,3));

hi_pwr = pwr1.Temp(:,:,1).*15 - (15-1).*pwr0.Temp(:,:,1);
diff = pwr15.Temp(:,:,1)-hi_pwr;

figure; imagesc( pwr1.Temp(:,:,1));
figure; imagesc( hi_pwr );
figure; imagesc( diff);
figure; imagesc( pwr15.Temp(:,:,1));

neg_pwr = pwr1.Temp(:,:,1).*(-3) - (-3 - 1).*pwr0.Temp(:,:,1);
diff = neg3.Temp(:,:,1)-neg_pwr;

figure; imagesc( pwr1.Temp(:,:,1));
figure; imagesc( neg_pwr );
figure; imagesc( diff);
figure; imagesc( neg3.Temp(:,:,1));