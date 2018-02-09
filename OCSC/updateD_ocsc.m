function [d,d_hat,s,y]=updateD_ocsc(para,Ahis,Bhis,s,y,d_hat)
% admm: d,s, dual y (unscaled)
if isempty(s)||isempty(y)
    if para.gpu==1
        if (para.precS ==1)
            s = zeros(para.size_k_full(1),para.size_k_full(2),para.K,'single','gpuArray');
        else
            s = zeros(para.size_k_full(1),para.size_k_full(2),para.K,'gpuArray'); 
        end
    else
        s = zeros(para.size_k_full(1),para.size_k_full(2),para.K);
        if (para.precS ==1)
           s=single(s);
        end
    end
    y=s;
end
s_hat = fft2(s);
y_hat = fft2(y);
%%
d_length = para.size_k_full(1)*para.size_k_full(2)*para.K;
rho = para.rho_D;
for i_d = 1:para.max_it_d
    d_hat = solve_conv_term_D(s_hat,y_hat,Ahis,Bhis, para,rho);
    d = real(ifft2( d_hat ));
    sold =s;
    s = prox_ind( d + y/rho,para);
    s_hat = fft2(s);
    y = y + rho* (d - s);
    y_hat = fft2(y);
    % stopping criteria
    ABSTOL = 1e-3;
    RELTOL = 1e-3; 
    h.r_norm(i_d) = norm(d(:)-s(:));
    h.s_norm(i_d) = norm(-rho*(s(:)-sold(:)));      
    h.eps_pri(i_d)=sqrt(d_length)*ABSTOL+RELTOL*max(norm(d(:)),norm(s(:)));
    h.eps_dual(i_d)=sqrt(d_length)*ABSTOL+RELTOL*norm(y(:));
        
    r1 = h.r_norm(i_d) < h.eps_pri(i_d);
    s1 = h.s_norm(i_d) < h.eps_dual(i_d);
    
	if  r1 && s1
        break;
    end   
end
end
%%
function dd_hat = solve_conv_term_D(ss_hat, yy_hat, Ahi,Bhi, par,rho_D)
    sy = par.size_z(1); sx = par.size_z(2); k = par.size_z(3);
    left = Bhi+permute( reshape(rho_D * ss_hat-yy_hat, sx * sy, 1, k), [3,2,1] );
    clear Bhi ss_hat yy_hat
    if par.gpu==1
        x2 = pagefun(@mtimes,Ahi,left);
    else
%         for i = 1:size(Ahi,3)
%             x2(:,:,i) = Ahi(:,:,i)*left(:,:,i);
%         end
        x2 = mtimesx(Ahi,left);
    end
    clear Ahi
    dd_hat = reshape(permute(x2,[3,1,2]),sy,sx,[]);
end