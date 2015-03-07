%aa=model_temp(143,140,:);
%aa=model_temp(143,135,:);
aa=model_temp(147,135,:);
cc=zeros(150,1);
dd=100;
for ii=1:150
    ix = ii * 100;
    bb=aa(1:ix);
    cc(ii)=var(bb);
    dd(end+1) = 100 * ii;
    if ii == 1
        dd = 100;
    end
   
end

figure; plot( dd,cc);