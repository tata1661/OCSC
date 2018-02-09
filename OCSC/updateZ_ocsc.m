function [z,z_hat]=updateZ_ocsc(x_hat,para,d_hat,var)
% admm: z, t, u (unscaled)
objective = @(z_hat) objective_online(z_hat,d_hat, x_hat,para );
if para.gpu==1
    if (para.precS ==1)
        z = randn(para.size_z,'single','gpuArray');
        t = zeros(para.size_z,'single','gpuArray');
    else
        z = randn(para.size_z,'gpuArray');
        t = zeros(para.size_z,'gpuArray'); 
    end
else
       z = randn(para.size_z);
       t = zeros(para.size_z);
    if (para.precS ==1)
       t=single(t);
       z = single(z);
    end
end
u = t;
z_hat = fft2(z);
t_hat =t;
u_hat = fft2(u);
rho_z = para.rho_Z;
%pre-stat for solve_conv_term_Z
xhat_flat = reshape(x_hat,1,[]);
PRE = [];
PRE.b = var.dhatT_flat.*xhat_flat;
PRE.sc = 1 ./(para.rho_Z + var.dhatTdhat_flat');
z_length = para.size_z(1)*para.size_z(2)*para.K;
if strcmp(para.verbose, 'inner') || strcmp( para.verbose, 'all')
    h.optval(1) = objective(z_hat);
    fprintf('start update Z \n--> Obj %3.3g \n',h.optval(1))
end
for i_z = 1:para.max_it_z
    z_hat = solve_conv_term_Z(var.dhatT_flat,t_hat, u_hat, para,rho_z,PRE);%%%%%
    z = real(ifft2( z_hat));
    told = t;
    t = prox_l1( z+u/rho_z, para.lambda(2)/rho_z );
    t_hat = fft2(t);
    u = u+rho_z*(z-t);
    u_hat = fft2(u);
    h.optval(i_z+1) = objective(z_hat);
    if strcmp(para.verbose, 'inner') || strcmp( para.verbose, 'all')        
        fprintf('-->inner iter_Z %d, Obj %3.3g \n', i_z, h.optval(i_z+1) )
    end
    % stopping criteria 
    ABSTOL = 1e-3;
    RELTOL = 1e-3;
    h.r_norm(i_z) = norm(z(:)-t(:));
    h.s_norm(i_z) = norm(-rho_z*(t(:)-told(:))); 
          
    h.eps_pri(i_z)=sqrt(z_length)*ABSTOL+RELTOL*max(norm(z(:)),norm(t(:)));
    h.eps_dual(i_z)=sqrt(z_length)*ABSTOL+RELTOL*norm(u(:)); 
    
    r1 = h.r_norm(i_z) < h.eps_pri(i_z);
    s1 = h.s_norm(i_z) < h.eps_dual(i_z);
        
    if r1 && s1
        break;
    end    
end  
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zz_hat = solve_conv_term_Z( dhatT, tt_hat, uu_hat, par,rho_Z,pre)
sy = par.size_z(1); sx = par.size_z(2);
b = pre.b+permute(reshape(rho_Z*tt_hat - uu_hat,sy*sx,[]),[2,1]);
clear tt_hat uu_hat
sc = pre.sc.* sum(conj(dhatT).*b, 1);
sc = dhatT.*sc;
zz_hat = reshape(permute(1/rho_Z *(b-sc), [2,1]), par.size_z);
end