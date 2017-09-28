pro plot_sphere_spectra_fixedtau

  ;; tau0_arr_st = ['4','5','6','7'];,'8']
  ;; tau0_arr    = [1.d4,1.d5,1.d6,1.d7];,1.d8]

  tau0_arr_st = ['6','6','6','6','6','6','6','6','6']      ;,'8']
  tau0_arr    = [1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6]  ;,1.d8]
   
  vth         = 128500.0d0       ; cm/s

  path = '/Users/tgarel/Rascas_tests/output/sphere_dom/'
  
  filearr = ['sphere_np1e6_T1e2_ndust0.0_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','sphere_np1e6_T1e2_taudust01_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/','sphere_np1e6_T1e2_taudust1_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/','sphere_np1e6_T1e2_taudust5_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/']

  ;xin_arr = [ 0,1,2,3,4,5,10,15,20]
  taud_arr = [0,0.1,1,5]

  PS_Start, File='plots/sphere_spectraDIJK_temp100K_tau6_vary_taudust.ps',nomatch=1,font=0

  device, helvetica=1,/bold
  device, isolatin1=1,/bold

  ct = 39

  loadct,ct 

  !x.style=1
  !y.style=1
  !x.thick=7
  !y.thick=7
  !x.charsize=1.3
  !y.charsize=1.4
  !p.thick=7
  !p.charthick=20

  for i=0,n_elements(filearr)-1 do begin

     myfile = path+filearr[i]+'photons_out.dat'
     
     print,myfile
     tau0 = tau0_arr[i]
 

     temp        = (vth / 12.85d5)^2. * 1.d4
     a           = 4.71d-4 * (temp/1.d4)^(-0.5) ;(12.85 / b)
     b           = vth / 1.d5                   ; km/s
     print,temp,vth
     
     read_photons,myfile,vth,x_out,tau0,status,xlast,nb_abs,time,ids
     nphotons = n_elements(x_out)
     ;; Doppler units - x
     
     xx = dindgen(4000) / 3999. * 400. - 200.
     ii = where(status eq 1,ni)  ;; escape
     qq = where(status eq 2,nq)  ;; killed by dust

     print,ni,nq,double(ni)/double(nphotons)
     histo,x_out[ii],-250.,250.,1000,h

     ;if i eq 0 then begin
        plot,h.x,h.dn/h.dx/total(h.dn),xtitle='x',ytitle='J(x,'+greek('tau')+'!d0!n)',charthick=5,xr=[-40.,40.],/xs,thick=8,yr=[0.0001,max(h.dn/h.dx/total(h.dn))*1.3],/nodata ;,/ylog
        ;legendold,['Iso, no Recoil, D/H=0','Iso, Recoil, D/H=0','Aniso, no Recoil, D/H=0','Iso, no Recoil, D/H=3e-5'],textcolors=indgen(4)*80,box=0,/left,charsize=0.9,spacing=1.2
     ;endif

        legendold,[greek('tau')+'!ddust!n = '+string(taud_arr[i],format='(f5.3)'),'f!desc!n = '+string(double(ni)/double(nphotons),format='(f5.3)')],textcolors=[i*80,i*80],box=0,/left,charsize=1.3,spacing=1.9
                
       ; legendold,['x!din!n = '+string(xin_arr[i])],textcolors=i*30,box=0,/left,charsize=1.3,spacing=1.2
        ;; if i eq 0 then legendold,['Iso, no Recoil, D/H=0'],textcolors=i*80,box=0,/left,charsize=1.3,spacing=1.2
        ;; if i eq 1 then legendold,['Iso, Recoil, D/H=0'],textcolors=i*80,box=0,/left,charsize=1.3,spacing=1.2
        ;; if i eq 2 then legendold,['Aniso, no Recoil, D/H=0'],textcolors=i*80,box=0,/left,charsize=1.3,spacing=1.2
        ;; if i eq 3 then legendold,['Iso, no Recoil, D/H=3e-5'],textcolors=i*80,box=0,/left,charsize=1.3,spacing=1.2

     oplot,h.x,h.dn/h.dx/total(h.dn),thick=5,color=i*80,psym=10 ;,linestyle=i

     ;; Exact Eq. 9 from Dijkstra06 - his tau0 is same as us
     j_x = sqrt(!pi) / sqrt(24.) / a / tau0 * xx * xx * (1.d0 / (1. + cosh(sqrt(2.*!pi^3./27.) * abs(xx*xx*xx)/ a / tau0)))
     x_peak = 0.92 * (a * tau0)^(1./3.)
     oplot,xx,j_x*2.*!pi,color=0,linestyle=2,thick=3 ;,psym=10
     ;; plots,[-1.*x_peak,-1.*x_peak],[0.,1.*max(j_x*2.*!pi)],linestyle=2.,color=254
     ;; plots,[1.*x_peak,1.*x_peak],[0.,1.*max(j_x*2.*!pi)],linestyle=2,color=254


     print,total(j_x*2.*!pi*(xx[1]-xx[0]))
     
     
     ;; Velocity
  
     ;; v_out = b * x_out * (-1.)
     ;; v     = dindgen(100000) / 99999. * 10000. - 5000. ; km/s
     ;; dv    = v(2)-v(1)
     
     ;; ;window,1
     ;; histo,v_out,-2000.,2000.,25000,h
     ;; plot,h.x,h.dn/h.dx/total(h.dn),psym=10,xr=[-50.,50.],/xs
     ;; j_v    = sqrt(!pi) / sqrt(24.) * v * v / a / b / b / tau0 * 1.d0 / (1.0 + cosh(sqrt(2.*!pi^3./27.) * abs(v*v*v)/a / b^3. / tau0))
     ;; v_peak = b * 0.92 * (a * tau0)^(1./3.)
     
     
     ;; oplot,v*(-1.),j_v*2.*!pi / b,color=254,linestyle=2,psym=10
     ;; plots,[-1.*v_peak,-1.*v_peak],[0.,1.*max(j_v*2.*!pi / b)],linestyle=2.,color=254
     ;; plots,[1.*v_peak,1.*v_peak],[0.,1.*max(j_v*2.*!pi / b)],linestyle=2,color=254
     
     
  endfor

  legendold,[greek('tau')+'!d0!n = 10!u6!n','T=100 K'],box=0,/right,charsize=1.3,spacing=1.8

 ; legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u5!n',greek('tau')+'!d0!n = 10!u6!n',greek('tau')+'!d0!n = 10!u7!n'],textcolors=(1+indgen(4))*50,box=0,/right,charsize=1.4

  
  ;legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u5!n',greek('tau')+'!d0!n = 10!u6!n',greek('tau')+'!d0!n = 10!u7!n',greek('tau')+'!d0!n = 10!u8!n'],textcolors=(1+indgen(5))*50,box=0,/right,charsize=1.4
  
  PS_End, /PDF, /Delete_PS
  
  
end