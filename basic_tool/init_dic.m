function [d_small] = init_dic(p)
s = RandStream.create('mt19937ar','seed',0);
RandStream.setGlobalStream(s);

if (p.precS ==1)
    d_small = single(randn(p.kernel_size)); %***
else
    d_small = randn(p.kernel_size); %***
end
% Initialize the first one to a uniform value
if size(p.psf_s,2) == 2
    d_small(:,:,1) = 1/(p.psf_s(1)*p.psf_s(2))*ones(p.psf_s(1),p.psf_s(2));
else
    d_small(:,:,1) = 1/(p.psf_s^2)*ones(p.psf_s,p.psf_s);  
end
%  Normalize the filter bank
d_small = prox_d_small(d_small);
% for i_filter = 1:p.K
%     d_small(:,:,i_filter) = d_small(:,:,i_filter)/(norm(d_small(:,:,i_filter)));
% end
end