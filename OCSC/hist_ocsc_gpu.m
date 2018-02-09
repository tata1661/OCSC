function [hisA_mat,hisB_mat] = hist_ocsc_gpu(bhat,zhat, par,hisA_mat,hisB_mat,b_no)
p = par.size_x(1) * par.size_x(2);
zhat_flat = reshape( zhat, p, [] );
zhatT_flat = conj(zhat_flat.'); % exactly zhat_flat'
bhat_flat = reshape(bhat,1,[]);
new_zb = reshape(zhatT_flat.*bhat_flat,par.K,1,[]);
%%
clear bhat_flat bhat
if (b_no==1) 
   %hisA_mat = zeros(par.K, par.K, p); 
    if (par.precS ==1)
        hisA_mat = zeros(par.K, par.K, par.size_x(1)*par.size_x(2),'single','gpuArray');  %***
    else
        hisA_mat = zeros(par.K, par.K, par.size_x(1)*par.size_x(2),'gpuArray');  %***
    end
    hisB_mat= new_zb;
    clear new_zb
    zhatTzhat_flat = sum(conj(zhat_flat).*zhat_flat,2);
    sc1 = 1 ./(par.rho_D*(par.rho_D + zhatTzhat_flat'));
    clear zhatTzhat_flat clear zhat zhat_flat

    z_inv1 = zhatT_flat.*sc1;
    z_inv2 = conj(zhatT_flat);
    z_inv1= reshape(z_inv1,par.K,1,[]);
    z_inv2 = reshape(z_inv2,1,par.K,[]); 
    clear sc1 zhatT_flat zhat_flat

    if (par.precS ==1)
        it = eye(par.K,'single','gpuArray')/par.rho_D;
    else
        it = eye(par.K,'gpuArray')/par.rho_D;
    end
    
 
    for i = 1:size(z_inv1,3)
        z_inv_sm = z_inv1(:,:,i)*z_inv2(:,:,i);
        z_inv_sm = it-z_inv_sm;
        hisA_mat(:,:,i)= z_inv_sm;
    end
%     z_inv = pagefun(@mtimes,z_inv1,z_inv2); %%%%%
%     hisA_mat =it-z_inv;    
else
    hisB_mat=hisB_mat+ new_zb;
    clear new_zb
    clear zhat zhat_flat
    zhat_f = reshape(zhatT_flat,par.K,1,[]);
    zhatT_f =reshape(conj(zhatT_flat),1,par.K,[]); 
    down = pagefun(@mtimes,zhatT_f,hisA_mat);
    down = pagefun(@mtimes,down,zhat_f);
    down = down+1;
    up1 = pagefun(@mtimes,hisA_mat,zhat_f);
    up2 = pagefun(@mtimes,zhatT_f,hisA_mat);
    clear zhatT_flat zhatT_f zhat_f 
    updown =  up2./down;
    clear up2 down
    hisA_mat= hisA_mat-pagefun(@mtimes,up1, updown);
end
end
