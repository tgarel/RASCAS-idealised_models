function r_function,u
 xin =       10.0000
 a =     0.0066690859
 xout=       50.000000
xtab = [xin,xout]
y = 1. / !pi^(3./2.) * exp(-u^2) * (atan((min(xtab)+u)/a)-atan(((max(xtab)-u)/a)))
 return,y
end
