function [stat_Z] = precompute_H_hat_Z(dhat, par )
stat_Z = [];
dhat_flat = reshape( dhat, par.size_x(1) * par.size_x(2), [] );
stat_Z.dhatTdhat_flat = sum(conj(dhat_flat).*dhat_flat,2);
stat_Z.dhatT_flat = conj(dhat_flat.'); 
end