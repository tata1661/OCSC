%Prox for kernel constraints
function [d_out,d_hat_out] = prox_ind( d, par)
%Get support
d_proj = d2dsmall(d,par);

%Normalize
d_proj = prox_d_small(d_proj);
%Now shift back and pad again
d_proj = dsmall2d(d_proj,par);

d_out = d_proj;
d_hat_out = fft2(d_proj);
end