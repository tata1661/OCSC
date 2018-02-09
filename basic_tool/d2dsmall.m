function d_small = d2dsmall(d,para)
if size(para.psf_s,2) == 2
    d_small = circshift( d, [para.psf_radius,0] ); 
    %d_small = d_small(1:para.psf_radius*2+1,:);
      
    if mod(para.psf_s(1),2)==1
        left_idx = 1:para.psf_radius(1)*2+1;
    else
        left_idx = 1:para.psf_radius(1)*2+1;
    end
        
    if mod(para.psf_s(2),2)==1
        right_idx = 1:para.psf_radius(2)*2+1;
    else
        right_idx = 1:para.psf_radius(2)*2;
    end
    
    d_small = d_small(left_idx,right_idx,:);

else 
    d_small = circshift( d, [para.psf_radius, para.psf_radius,0] ); 
    d_small = d_small(1:para.psf_radius*2+1, 1:para.psf_radius*2+1,:);    
end
end