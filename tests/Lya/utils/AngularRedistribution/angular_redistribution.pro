pro angular_redistribution,iso=iso,dip=dip,ray=ray

  if keyword_set(iso) then choice = 'isotropic'
  if keyword_set(dip) then choice = 'anisotropic'
  if keyword_set(ray) then choice = 'Rayleigh'

  temp = '1.E+04'
  startpdf,'plots/tests_hi_test_angular_redist_T'+temp+'K_'+choice+'_bis'
  loadct,39

  ;; mycol = [cgcolor('cg3'),cgcolor('chartreuse'),cgcolor('ygb3'),210,cgcolor('red4'),cgcolor('violet'),cgcolor('gray')]
  ;; xin_arr = [0,1,2,3,4,5];,8]

    mycol = [cgcolor('cg3'),cgcolor('chartreuse'),cgcolor('ygb3'),210,cgcolor('red4'),cgcolor('violet')];,cgcolor('gray')]
    xin_arr = [0,1,2,3,4,5];,8]
    
  !Y.Margin = [5, 5]
  !X.Margin = [7, 7]

  !p.multi=[0,4,3]
  !p.thick=5
  !x.thick=4
  !y.thick=4
  
  for j=0,n_elements(xin_arr)-1 do begin
     x_in = xin_arr(j)
     
     nphotons = 0L
     xx  = -1.
     kk1 = -1.
     kk2 = -1.
     kk3 = -1.
        
     if keyword_set(iso) then begin
        outfile = '../../LyaTests_output_data/SingleScatt_H/tests_hi_xin'+strtrim(string(x_in,format='(i11)'),2)+'_T'+strtrim(temp,2)+'_Isotropic_NoRecoil.dat'
     endif
     if keyword_set(dip) then begin
        outfile = '../../LyaTests_output_data/SingleScatt_H/tests_hi_xin'+strtrim(string(x_in,format='(i1)'),2)+'_T'+strtrim(temp,2)+'_Dipolar_NoRecoil.dat'
     endif
     if keyword_set(ray) then begin
        outfile = '../../LyaTests_output_data/SingleScatt_H/tests_hi_xin'+strtrim(string(x_in,format='(i1)'),2)+'_T'+strtrim(temp,2)+'_Rayleigh_NoRecoil.dat'
     endif
     
     openr,11,outfile
     readf,11, nphotons, x_in, T, nu_D
     print,'Nphotons = ',nphotons
     xout   = dblarr(nphotons)
     kout1  = xout
     kout2  = xout
     kout3  = xout
     for i=0L,nphotons-1L do begin
        readf,11,xx,kk1,kk2,kk3
        xout[i]  = xx
        kout1[i] = kk1
        kout2[i] = kk2
        kout3[i] = kk3
     endfor
     close,11
     
     phi       = atan(kout2/kout1)
     cos_phi   = cos(phi)
     cos_theta = kout3
     theta     = acos(cos_theta)
     

     histo,phi,-2,2.,130.,h
     plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[-2.,2.],/nodata,xtitle=greek('Phi'),ytitle='P('+greek('Phi')+')',title='x!din!n = '+string(x_in,format='(f4.1)')+' - '+choice+' case',charsize=0.8,yr=[0.,max(h.dn/total(h.dn)/h.dx)*1.2],thick=4 ;,yr=[0.0001,1.1],/ylog
     oplot,h.x,h.dn/total(h.dn)/h.dx,psym=10,color=mycol[j],thick=4
     ;; ,ytitle='dN/N!dtot!n/bin'
     
     histo,cos_theta,-1.5,1.5,120.,h
     plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[-1.5,1.5],/nodata,xtitle='cos '+greek('theta'),ytitle='P(cos '+greek('theta')+')',title='x!din!n = '+string(x_in,format='(f4.1)')+' - '+choice+' case',charsize=0.8 ,yr=[0.,max(h.dn/total(h.dn)/h.dx)*1.5],thick=4 ;,yr=[0.0001,1.1],/ylog
     oplot,h.x,h.dn/total(h.dn)/h.dx,psym=10,color=mycol[j],thick=4
     ;; if j eq 0 then begin
     ;;    legendold,['Rayleigh scattering','Core scattering'],thick=4 ,charthick=6,charsize=0.4,box=0,linestyle=[2,1],linsize=1.5,/top,/left
     ;; endif
     
    ; print,total(h.dn/total(h.dn))

     if keyword_set(ray) then begin
        ;; Normalization of phase functions:
        ;; 0.5 * INT_-1^1 Phi(mu)dmu = 1 where mu = cos(theta)
        ;; Rayleigh
        cost = dindgen(1000) / 999. * 2. - 1.
        pp = 3./4. * (1.+cost*cost) / 2. ; / !pi
        oplot,cost,pp,linestyle=2,thick=6
       ; legendold,['Function phase for ','Rayleigh scattering'],linsize=0.6,colors=[0,0],textcolors=[0,0],thick=8 ,charthick=7,charsize=0.7,box=0,/top,/center,linestyle=[2,-8]
     endif

     if keyword_set(dip) then begin
        ;; core
        cost = dindgen(1000) / 999. * 2. - 1.
        pp = (11./12. + 3./12. * cost*cost) / 2. ; / !pi
        oplot,cost,pp,linestyle=1,thick=10
        ;legendold,['Function phase for ','Rayleigh scattering','core scattering'],linsize=0.6,colors=[0,0,0],textcolors=[0,0,0],thick=8 ,charthick=7,charsize=1.4,box=0,linestyle=[-8,2,1],position=[-0.7,0.15]
        ;; Also show Rayleigh for anisotr. case
         ;; 0.5 * INT_-1^1 Phi(mu)dmu = 1 where mu = cos(theta)
        ;; Rayleigh
        cost = dindgen(1000) / 999. * 2. - 1.
        pp = 3./4. * (1.+cost*cost) / 2. ; / !pi
       ; oplot,cost,pp,linestyle=2,thick=6
     endif
     
  endfor
  
  
  !p.multi=0
  
  
  stoppdf

end 
