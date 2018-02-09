function [hisA_mat,hisB_mat] = hist_ocsc_cpu(bhat,zhat, par,hisA_mat,hisB_mat,b_no)
p = par.size_x(1) * par.size_x(2);
zhat_flat = reshape( zhat, p, [] );
zhatT_flat = conj(zhat_flat.'); % exactly zhat_flat'
bhat_flat = reshape(bhat,1,[]);
new_zb = reshape(zhatT_flat.*bhat_flat,par.K,1,[]);
zhat_f = reshape(zhatT_flat,par.K,1,[]);
zhatT_f =reshape(conj(zhatT_flat),1,par.K,[]);
%%
if (b_no==1) 
    hisB_mat = new_zb;
    zhatTzhat_flat = sum(conj(zhat_flat).*zhat_flat,2);
    sc1 = 1 ./(par.rho_D*(par.rho_D + zhatTzhat_flat.'));
    clear zhatTzhat_flat clear zhat zhat_flat
    z_inv1= reshape(bsxfun(@times,zhatT_flat,sc1),par.K,1,[]);
    clear sc1 zhatT_flat zhat_flat
%     for i = 1:size(z_inv1,3)
%         z_inv_sm = z_inv1(:,:,i)*zhatT_f(:,:,i);
%         z_inv_sm = eye(par.K)/par.rho_D-z_inv_sm;
%         hisA_mat(:,:,i)=z_inv_sm;
%     end
    z_inv = mtimesx(z_inv1,zhatT_f);    
    hisA_mat =eye(par.K)/par.rho_D-z_inv;
else
hisB_mat=hisB_mat+ new_zb;
clear zhat zhat_flat
down = mtimesx(zhatT_f,hisA_mat);
down = mtimesx(down,zhat_f);
down = down+1;
up1 = mtimesx(hisA_mat,zhat_f);
up2 = mtimesx(zhatT_f,hisA_mat);
clear zhatT_flat zhatT_f zhat_f 
updown  = up2./ down;
clear up2 down
hisA_mat= hisA_mat-mtimesx(up1, updown);% can change into for to save mem
%%**** to save mem
% for i = 1:size(up1,3)
%     temp_sm = up1(:,:,i)*updown(:,:,i);
%     hisA_mat(:,:,i)=hisA_mat(:,:,i)-temp_sm;
% end
end
end
