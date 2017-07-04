pro hummer_integ_dipolar,xin

  ;; Num. integration of Hummer62 redist. function RII-B (Eq. 3.12.2)
  ;; Dipolar scatt., radiation damping, coherence in atom's frame (no recoil)
  
  ;; a = Gamma / 4 / pi / Delta_nuD = 4.7d-4 (T/1.d4)^(-1/2)
   a = 4.7d-4                    ; => T=10^4 K
  ;; a = 4.7d-3                    ; => T=10^2 K

  ;; Assume xout from -50 to 50 for num. convergence issues in the
  ;; (normalisation of the) integral
  ;xout_arr = dindgen(10000)/9999.d0*100.d0-50.d0

    xout_arr = dindgen(10000)/9999.d0*100.d0-50.d0
    if xin gt 5. then begin
       ;xout_arr = dindgen(1500)/1499.d0*100.d0-50.d0
       xout_arr = dindgen(600)/599.d0*100.d0-50.d0
    endif else begin
       xout_arr = dindgen(1000)/999.d0*100.d0-50.d0
    endelse
    
  for i=0L,n_elements(xout_arr)-1L do begin
     xout = xout_arr[i]
     openw,11,'r_function_dipolar.pro'
     printf,11,'function r_function_dipolar,u'
     printf,11,'xin = ',xin
     printf,11,'a = ',a
     printf,11,'xout=',xout
     printf,11,'xtab = [xin,xout]'
     printf,11,'xmin = min(xtab)'
     printf,11,'xmax = max(xtab)'  
     printf,11,'y = 3. / !pi^(3./2.) / 8. * a * exp(-u^2) * ((((3.*a^4 + 3.*u^4 - u^2*xout^2 - xin^2*(u^2 - 3.*xout^2) - a^2*(3.*xin^2 - 2.*u^2 + 12.*xin*xout + 3.*xout^2))* atan((xmin+u)/a) + a*(((xmin+u) - xout)*(-3.*a^2 + 3.*xin^2 - 2.*u^2 - 3.*xin*(xmin+u) + (xmin+u)^2 + 9.*xin*xout - 2*(xmin+u)*xout + xout^2) + (xin + xout)*(3.*a^2 + u^2 - 3.*xin*xout)*alog(a^2 + (xmin+u)^2)))/ (a*u^4)) - (((3.*a^4 + 3.*u^4 - u^2*xout^2 - xin^2*(u^2 - 3.*xout^2) - a^2*(3.*xin^2 - 2.*u^2 + 12.*xin*xout + 3.*xout^2))* atan((xmax-u)/a) + a*(((xmax-u) - xout)*(-3.*a^2 + 3.*xin^2 - 2.*u^2 - 3.*xin*(xmax-u) + (xmax-u)^2 + 9.*xin*xout - 2*(xmax-u)*xout + xout^2) + (xin + xout)*(3.*a^2 + u^2 - 3.*xin*xout)*alog(a^2 + (xmax-u)^2)))/ (a*u^4)))'
     printf,11,' return,y'
     printf,11,'end'
     close,11
     RESOLVE_ROUTINE, 'r_function_dipolar', /IS_FUNCTION
     xtab = [xin,xout]
     x1 = abs(max(xtab)-min(xtab))/2.
     if xin gt 5. then begin
        vv = qromo('r_function_dipolar',x1,/midexp,eps=1.d-7)
     endif else begin
        vv = qromo('r_function_dipolar',x1,/midexp,eps=1.d-7)
     endelse
     
     ;; Sum of R(xout,xin) over xout = Voigt(a,xin) (Furlanetto+06,
     ;; Eq. 39)  
     phi = voigt(a,xin) / sqrt(!pi)
     vv = vv / phi
     
     if i eq 0 then begin
        p_xout = vv
     endif else begin
        p_xout = [p_xout,vv]
     endelse
    
     
  endfor
  xout = xout_arr

  if xin lt 10 then begin
     save,filename='hummer_DIPOLAR_xin'+strtrim(string(xin,format='(i1)'),2)+'_temp10000K_eps7.sav',xout,p_xout
  endif else begin
     save,filename='hummer_DIPOLAR_xin'+strtrim(string(xin,format='(i2)'),2)+'_temp10000K_eps7.sav',xout,p_xout
  endelse

  plot,xout,p_xout,xr=[-5.,13.],yr=[0.,0.95]
  
end
