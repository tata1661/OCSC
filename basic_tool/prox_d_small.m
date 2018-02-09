function [d] = prox_d_small(d)
bound_norm = 1;
d_norm = repmat( sum(sum(d.^2, 1),2), [size(d,1), size(d,2), 1] );
d( d_norm >= bound_norm ) = d( d_norm >= bound_norm ) ./ sqrt(d_norm( d_norm >= bound_norm ));
end