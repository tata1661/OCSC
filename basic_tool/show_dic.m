function [] = show_dic(d,PARA,save_flag,sort_flag)
s1 = size(d,1);
if (isempty(who('sort_flag')))
    sort_flag = 1;
end

if sort_flag==1
    flat_d = reshape(d,[],size(d,3));
    vars = var(flat_d,0,1);
    [V,indices] = sort(vars,'descend');
end

pd = 1;
sqr_k = ceil(sqrt(size(d,3)));
d_disp = zeros( sqr_k * [s1 + pd, s1 + pd] + [pd, pd]);
for j = 0:size(d,3) - 1
    ind = j+1;
    if sort_flag==1
        ind =indices(j+1);
    end
    d_disp( floor(j/sqr_k) * (s1 + pd) + pd + ...
        (1:s1) , mod(j,sqr_k) * (s1 + pd) + pd + (1:s1) ) = d(:,:,ind); 
end

max_sig = max(reshape(d_disp,[],1));
min_sig = min(reshape(d_disp,[],1));
%Transform and save
new_d_disp = (d_disp - min_sig)/(max_sig - min_sig);

rowN = sqr_k;
clN = size(d,3)/(floor(size(d,3)/rowN));
new_d=new_d_disp;
for j = 1:clN+1
    new_d(:,(s1+pd)*(j-1)+pd) = 1;
end

for j = 1:ceil(size(d,3)/rowN)+1
    new_d((s1+pd)*(j-1)+pd,:) = 1;
end

%figure();
imshow(new_d)
%imshow(new_d_disp)
if save_flag==1
        imwrite(new_d , sprintf('%s/filter_%s.png',PARA.path,PARA.save_name),'bitdepth', 16);%%
end
end



