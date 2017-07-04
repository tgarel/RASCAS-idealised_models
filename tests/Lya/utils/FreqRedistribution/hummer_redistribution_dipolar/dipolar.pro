pro dipolar,xin,xout,a,fu

  ;; x=b=xin
  ;; t=x
  ;; xout=y

  ;; q1 and q2 are the primitives of
  ;; "(3-((b-x)/u)^2-((y-x)/u)^2+3*((b-x)/u)^2*((y-x)/u)^2)/(x^2+a^2)"
  ;; = Eq. 3.12.2 of Hummer62
  ;; done with wolfram integrator and evaluated at the bounds:
  ;; (xmin+u) and (xmax-u)
  
  q1 = ((3.*a^4 + 3.*u^4 - u^2*xout^2 - xin^2*(u^2 - 3.*xout^2) - a^2*(3.*xin^2 - 2.*u^2 + 12.*xin*xout + 3.*xout^2))* atan((xmin+u)/a) + a*(((xmin+u) - xout)*(-3.*a^2 + 3.*xin^2 - 2.*u^2 - 3.*xin*(xmin+u) + (xmin+u)^2 + 9.*xin*xout - 2*(xmin+u)*xout + xout^2) + (xin + xout)*(3.*a^2 + u^2 - 3.*xin*xout)*alog10(a^2 + (xmin+u)^2)))/ (a*u^4)

  
  q2 = ((3.*a^4 + 3.*u^4 - u^2*xout^2 - xin^2*(u^2 - 3.*xout^2) - a^2*(3.*xin^2 - 2.*u^2 + 12.*xin*xout + 3.*xout^2))* atan((xmax-u)/a) + a*(((xmax-u) - xout)*(-3.*a^2 + 3.*xin^2 - 2.*u^2 - 3.*xin*(xmax-u) + (xmax-u)^2 + 9.*xin*xout - 2*(xmax-u)*xout + xout^2) + (xin + xout)*(3.*a^2 + u^2 - 3.*xin*xout)*alog10(a^2 + (xmax-u)^2)))/ (a*u^4)

  fu = q1-q2
  
  ;; ((3*a^4 + 3*u^4 - u^2*y^2 - b^2*(u^2 - 3*y^2) - a^2*(3*b^2 - 2*u^2 + 12*b*y + 3*y^2))* ArcTan[x/a] + a*((x - y)*(-3*a^2 + 3*b^2 - 2*u^2 - 3*b*x + x^2 + 9*b*y - 2*x*y + y^2) + (b + y)*(3*a^2 + u^2 - 3*b*y)*Log[a^2 + x^2]))/ (a*u^4)

  ;;   ((3.*a^4 + 3.*u^4 - u^2*y^2 - b^2*(u^2 - 3.*y^2) - a^2*(3.*b^2 - 2.*u^2 + 12.*b*y + 3.*y^2))* atan(x/a) + a*((x - y)*(-3.*a^2 + 3.*b^2 - 2.*u^2 - 3.*b*x + x^2 + 9.*b*y - 2*x*y + y^2) + (b + y)*(3.*a^2 + u^2 - 3.*b*y)*alog10(a^2 + x^2)))/ (a*u^4)

end
