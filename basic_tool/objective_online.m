function [f_val,f_z] = objective_online(z_hat,d_hat,x_hat, par)
Dz_hat = sum(d_hat.* z_hat,3);
z = real(ifft2(z_hat));
f_z = par.lambda(1) * 1/2 /(par.size_z(1)*par.size_z(2))* norm(Dz_hat(:) - x_hat(:), 2)^2;

% % equivalent to 
% x = real(ifft2(x_hat));
% Dz = real(ifft2( Dz_hat));
% f_z2 = par.lambda(1) * 1/2 * norm( reshape( Dz - x, [], 1) , 2 )^2;

g_z = par.lambda(2) * sum( abs( z(:) ), 1 );
f_val = f_z + g_z;
end
