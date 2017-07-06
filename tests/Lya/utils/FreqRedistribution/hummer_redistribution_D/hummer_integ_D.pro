pro hummer_integ_D,xin

  ;; Num. integration of Hummer62 redist. function RII-A (Eq. 3.12.1)
  ;; Isotropic scatt., radiation damping, coherence in atom's frame (no recoil)
  
  ;; a = Gamma / 4 / pi / Delta_nuD = 4.7d-4 (T/1.d4)^(-1/2)
  ;; a = 4.7d-4 ; => T=10^4 K
  ;; a = 4.7d-3  ; => T=10^2 K

  T = 1.d2
  kb                     = 1.38064852d-16
  mp                     = 1.672621898d-24         ; [g] proton mass
  sqrt_H2Deut_mass_ratio = 0.7071067811865d0
  gamma_over_fourpi      = 1215.34d0 
  vth                    = sqrt(2.d0*kb*T/mp) ; cm/s
  vth                    = vth * sqrt_H2Deut_mass_ratio ; cm/s

  lambda_0               = 1215.34d0  
  lambda_0_cm            = lambda_0 / 1.d8
  delta_nu_D             = vth / lambda_0_cm
  a                      = 6.265d8 / 4.d0 / !pi / delta_nu_D
  
  ;; Assume xout from -50 to 50 for num. convergence issues in the
  ;; (normalisation of the) integral
  ;xout_arr = dindgen(10000)/9999.d0*100.d0-50.d0

  xout_arr = dindgen(10000)/9999.d0*100.d0-50.d0
  ;xout_arr = dindgen(10000)/9999.d0*300.d0-150.d0

  for i=0L,n_elements(xout_arr)-1L do begin
     xout = xout_arr[i]
     openw,11,'r_function.pro'
     printf,11,'function r_function,u'
     printf,11,' xin = ',xin
     printf,11,' a = ',a
     printf,11,' xout=',xout
     printf,11,'xtab = [xin,xout]'
     printf,11,'y = 1. / !pi^(3./2.) * exp(-u^2) * (atan((min(xtab)+u)/a)-atan(((max(xtab)-u)/a)))'
     printf,11,' return,y'
     printf,11,'end'
     close,11
     RESOLVE_ROUTINE, 'r_function', /IS_FUNCTION
     xtab = [xin,xout]
     x1 = abs(max(xtab)-min(xtab))/2.
     vv = qromo('r_function',x1,/midexp,eps=1.d-10)

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
     save,filename='hummer_xin'+strtrim(string(xin,format='(i1)'),2)+'_temp100K_D.sav',xout,p_xout
  endif else begin
     save,filename='hummer_xin'+strtrim(string(xin,format='(i2)'),2)+'_temp100K_D.sav',xout,p_xout
  endelse
  
 ; save,filename='hummer_xin'+strtrim(string(xin,format='(i2)'),2)+'.sav',xout,p_xout
 ; save,filename='hummer_xin'+strtrim(string(xin,format='(i3)'),2)+'_max150.sav',xout,p_xout

end
