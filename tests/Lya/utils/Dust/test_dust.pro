pro test_dust

  ;; In case you dont have some of the functions called in this routine....:
  ;; IDL_PATH=$IDL_PATH:+":+/gridsoft/ipnls/idl/idl/idl81/lib/:+/home/cral/garel/idl_routines/cosmo/:+/home/cral/garel/idl_routines/Esrg/:+/home/cral/garel/idl_routines/ajs/:+/home/cral/garel/idl_routines/astron/:+/home/cral/garel/idl_routines/coyote/:+/home/cral/garel/idl_routines/cosmo/:+/home/cral/garel/idl_routines/:+/cral2/garel/lyon/Lya/sam_lya_model/:"

  temp = 1.d4
  if temp eq 1.d4 then temp_st = '_temp1d4'
  if temp eq 1.d2 then temp_st = '_temp1d2'

  startpdf,'plots/test_dust'+temp_st

  albedo         = 0.46
  g_dust         = 0.73

  !Y.Margin = [6, 3]
  !X.Margin = [7, 7]  
  !p.multi=[0,3,3]

  !p.thick=5
  !x.thick=4
  !y.thick=4
  
  loadct,39

 ; mycol = [cgcolor('cg3'),cgcolor('chartreuse'),cgcolor('ygb3'),210,cgcolor('red4'),cgcolor('violet'),cgcolor('gray')]
  mycol = [cgcolor('cg3'),210,cgcolor('violet')]
    
  xin_arr = [0,3,5] ;,1,2,3,4,5]
  
  for j=0,n_elements(xin_arr)-1 do begin
     x_in = xin_arr(j)
     
     nphotons = 0L
     xx  = -99.
     kk1 = -99.
     kk2 = -99.
     kk3 = -99.
     iabs = -99
     outfile = '../../LyaTests_output_data/SingleScatt_dust/tests_dust_xin'+strtrim(string(x_in,format='(i1)'),2)+'_HG41.dat'

     openr,11,outfile
     readf,11, nphotons, x_in, T, nu_D
     print,'Nphotons = ',nphotons
     xout   = dblarr(nphotons)
     kout1  = xout
     kout2  = xout
     kout3  = xout
     ilost  = intarr(nphotons)
     
     for i=0L,nphotons-1L do begin
        readf,11,xx,kk1,kk2,kk3,iabs
        xout[i]  = xx
        kout1[i] = kk1
        kout2[i] = kk2
        kout3[i] = kk3
        ilost[i] = iabs
     endfor
     close,11
     
     oo = where(ilost eq 0,no)
     print,'Fesc = ',double(no)/double(nphotons),' Albedo = ',albedo

     mm = where(ilost gt 1 and ilost lt 0,nm)
     if nm ne 0 then stop
     
     histo,xout[oo],-10.,10.,2000,h
    
     plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[x_in*(-1.)-1.,x_in+1.],/nodata,xtitle='x!dout!n',ytitle='R!ddust!n(x!din!n,x!dout!n)',charsize=0.9,yr=[0.,max(h.dn/total(h.dn)/h.dx)*1.6],title='x!din!n = '+string(x_in,format='(f4.1)') ; ,yr=[0.0001,1.1],/ylog;,yr=[0.,0.99]
     oplot,h.x,h.dn/total(h.dn)/h.dx,psym=10,color=mycol[j],thick=5
     legendold,['Escaped photons'],thick=6 ,charthick=6,charsize=0.6,box=0,linestyle=[0],linsize=0.7,/top,/left,colors=[mycol[j]],textcolors=[mycol[j]]
     legendold,['g!ddust!n = '+string(g_dust,format='(f5.2)'),'albedo = '+string(albedo,format='(f5.2)'),'F!desc!n = '+string(double(no)/double(nphotons),format='(f7.3)')],thick=6 ,charthick=6,charsize=0.6,box=0,/center,/left,spacing=1.3
     
     phi       = atan(kout2/kout1)
     cos_phi   = cos(phi)
     cos_theta = kout3
     theta     = acos(cos_theta)
     
     histo,phi[oo],-2,2.,130.,h
     plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[-2.,2.],/nodata,xtitle=greek('Phi'),ytitle='P('+greek('Phi')+')',title='x!din!n = '+string(x_in,format='(f4.1)'),charsize=0.9,yr=[0.,max(h.dn/total(h.dn)/h.dx)*1.3] ;,yr=[0.0001,1.1],/ylog
     oplot,h.x,h.dn/total(h.dn)/h.dx,psym=10,color=mycol[j],thick=6
     legendold,['Escaped photons'],thick=6 ,charthick=6,charsize=0.6,box=0,linestyle=[0],linsize=0.7,/top,/left,colors=[mycol[j]],textcolors=[mycol[j]]

     
     histo,cos_theta[oo],-1.5,1.5,300.,h
     plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[-1.5,1.5],/nodata,xtitle='cos '+greek('theta'),ytitle='P(cos '+greek('theta')+')',title='x!din!n = '+string(x_in,format='(f4.1)'),charsize=0.9 ,yr=[0.01,max(h.dn/total(h.dn)/h.dx)*5.5],/ylog
     oplot,h.x,h.dn/total(h.dn)/h.dx,psym=10,color=mycol[j],thick=6   
     legendold,['Escaped photons','HG41 phase function'],thick=6 ,charthick=6,charsize=0.6,box=0,linestyle=[0,0],linsize=0.5,/top,/left,colors=[mycol[j],0],textcolors=[mycol[j],0],spacing=1.1

     ;; HG phase function: P(cost) = (1.-g_dust^2.) / (1.+g_dust^2. - 2.*g_dust*cost)^1.5 / 4. / !pi
     ;; Int_(-1)^1 pp dcos(theta) = 1/(2pi)
     cost = dindgen(10000) / 9999. * 2. - 1.
    ; cost = dindgen(120) / 119. * 3. - 1.5
     pp = (1.-g_dust^2.) / (1.+g_dust^2. - 2.*g_dust*cost)^1.5 / 4. / !pi * 2. * !pi 
     oplot,cost,pp,linestyle=2,thick=4

     print,total(h.dn)
     
     ;print,total(h.dn/total(h.dn)/h.dx*h.dx),total(pp*(cost(1)-cost(0)))
    ; print,h.dx

     ;print,min(cos_theta),max(cos_theta)
  endfor

  !p.multi=0

  stoppdf
 ; print,1./total(h.dn)/h.dx

end 
