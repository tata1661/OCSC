function x = prox_l1(v, lambda)
s = sign(v);
v = abs(v);
x = max(0, v - lambda) .* s;
end