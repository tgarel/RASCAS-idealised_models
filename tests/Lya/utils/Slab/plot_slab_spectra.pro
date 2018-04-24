pro plot_slab_spectra

  tau0_arr_st = ['4','5','6','7']
  tau0_arr    = [1.d4,1.d5,1.d6,1.d7]

  ;; tau0_arr_st = ['4','7']      ;,'8']
  ;; tau0_arr    = [1.d4,1.d7]  ;,1.d8]
   
 ;; vth         = 128500.0d0       ; cm/s
  vth = 40635.269d0             ; cm/s

  PS_Start, File='plots/tauHall_neufeld_temp1d1_xin0_oneplot.ps',nomatch=1,font=0
 ; PS_Start, File='plots/test_thickness.ps',nomatch=1,font=0
  
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

  col = [70,254,cgcolor('lime green'),210];,30]
    
  xxi = 0. ; input frequency
  
  for i=0,n_elements(tau0_arr)-1 do begin

     ;if i eq 0 then myfile = '/Users/tgarel/Rascas/output/Lya_tests/slab/tauH'+tau0_arr_st[i]+'_temp1d2_xin0/photons_out.dat'
     myfile = '/Users/tgarel/Rascas/output/Lya_tests/slab/tauH'+tau0_arr_st[i]+'_temp1d1_xin0_hotfix_width01/photons_out.dat'
     
     print,myfile
     tau0 = tau0_arr[i]
 
     temp        = (vth / 12.85d5)^2. * 1.d4
     a           = 4.71d-4 * (temp/1.d4)^(-0.5) ;(12.85 / b)
     b           = vth / 1.d5                   ; km/s
     print,temp,vth
     
     read_photons,myfile,vth,x_out,tau0,status,xlast,nb_abs,time,ids
     
     ;; Doppler units - x
     
     xx = dindgen(4000) / 3999. * 400. - 200.
     histo,x_out,-250.,250.,480,h

     if i eq 0 then begin
        plot,h.x,h.dn/h.dx/total(h.dn),xtitle='x',ytitle='J(x,'+greek('tau')+'!d0!n)',charthick=5,xr=[-115.,115.],/xs,thick=8,yr=[0.0001,max(h.dn/h.dx/total(h.dn))*1.1],/nodata ;,/ylog
        ;;legendold,['R!dslab!n within 1 cell'],box=0,/left,charsize=1.3,spacing=2.1
     endif

     oplot,h.x,h.dn/h.dx/total(h.dn),thick=6,color=col[i],psym=10; ,linestyle=i

     ;; Exact Eq. 39 from Smith15 - his tau0 is same as us
     j_x = 1. / 4. / sqrt(6.*!pi) / a / tau0 * xx * xx * (1.d0 / (cosh(sqrt(!pi^3./54.) * abs(xx*xx*xx - xxi*xxi*xxi)/ a / tau0)))

     x_peak = 1.066 * (a * tau0)^(1./3.)
     oplot,xx,j_x*4.*!pi,color=0,linestyle=2,thick=3 ;,psym=10
     ;; plots,[-1.*x_peak,-1.*x_peak],[0.,1.*max(j_x*2.*!pi)],linestyle=2.,color=254
     ;; plots,[1.*x_peak,1.*x_peak],[0.,1.*max(j_x*2.*!pi)],linestyle=2,color=254

     print,total(j_x*4.*!pi*(xx[1]-xx[0]))
     
    ; legendold,[greek('tau')+'!d0!n = 10!u'+tau0_arr_st[i]+'!n'],textcolors=col[i],box=0,/right,charsize=1.6
   ;  legendold,['Hydrogen','T=10!u2!n K'],/center,/right,charthick=7,charsize=1.3,box=0,spacing=1.9
  endfor


    legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u5!n',greek('tau')+'!d0!n = 10!u6!n',greek('tau')+'!d0!n = 10!u7!n'],textcolors=col,box=0,/right,charsize=1.7

    
 ;; legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u5!n',greek('tau')+'!d0!n = 10!u6!n',greek('tau')+'!d0!n = 10!u7!n',greek('tau')+'!d0!n = 10!u8!n'],textcolors=col,box=0,/right,charsize=1.4
  
  PS_End, /PDF, /Delete_PS
  
  
end
