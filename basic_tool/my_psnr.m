function [p,rmse] = my_psnr(b,Dz)
it = size(Dz,1);
tmp =  norm(b(:) - Dz(:));
p = 20*log10(it /tmp);

rmse = sqrt( 1/length(b(:)) * (tmp^2));
end
