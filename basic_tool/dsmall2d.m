function d = dsmall2d(d_small,para)
    d = padarray( d_small, [para.size_x(1) - para.kernel_size(1), para.size_x(2) - para.kernel_size(2), 0], 0, 'post');
    
    if size(para.psf_radius,2) == 2
        d = circshift(d, -[para.psf_radius(1), para.psf_radius(2), 0] );      
    else
        d = circshift(d, -[para.psf_radius, para.psf_radius, 0] );
    end
end