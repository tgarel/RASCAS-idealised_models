function r_function_dipolar,u
xin =       9.00000
a =    0.00047000000
xout=       50.000000
xtab = [xin,xout]
xmin = min(xtab)
xmax = max(xtab)
y = 3. / !pi^(3./2.) / 8. * a * exp(-u^2) * ((((3.*a^4 + 3.*u^4 - u^2*xout^2 - xin^2*(u^2 - 3.*xout^2) - a^2*(3.*xin^2 - 2.*u^2 + 12.*xin*xout + 3.*xout^2))* atan((xmin+u)/a) + a*(((xmin+u) - xout)*(-3.*a^2 + 3.*xin^2 - 2.*u^2 - 3.*xin*(xmin+u) + (xmin+u)^2 + 9.*xin*xout - 2*(xmin+u)*xout + xout^2) + (xin + xout)*(3.*a^2 + u^2 - 3.*xin*xout)*alog(a^2 + (xmin+u)^2)))/ (a*u^4)) - (((3.*a^4 + 3.*u^4 - u^2*xout^2 - xin^2*(u^2 - 3.*xout^2) - a^2*(3.*xin^2 - 2.*u^2 + 12.*xin*xout + 3.*xout^2))* atan((xmax-u)/a) + a*(((xmax-u) - xout)*(-3.*a^2 + 3.*xin^2 - 2.*u^2 - 3.*xin*(xmax-u) + (xmax-u)^2 + 9.*xin*xout - 2*(xmax-u)*xout + xout^2) + (xin + xout)*(3.*a^2 + u^2 - 3.*xin*xout)*alog(a^2 + (xmax-u)^2)))/ (a*u^4)))
 return,y
end
