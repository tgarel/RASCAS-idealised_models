pro mean_free_path

  ;; To test MFP, 'kill' photons after their 1st HI scattering
  ;; l. 411 in module_photon.f90:
  ;; call gas_scatter(scatter_flag, cell_gas, nu_cell, k, nu_ext, iran) 
  ;; !! TIBO: test mfp
  ;;  scatter_flag=4
  ;; !! OBIT   
  
  tau0_arr_st = ['4','5','6','7','8']
  tau0_arr    = [1.d4,1.d5,1.d6,1.d7,1.d8]
  nhi_arr     = [1.,10.,100.,1000.,10000.]
  nhi_arr_st  = ['1e0','1e1','1e2','1e3','1e4']
  
  ;; tau0_arr_st = ['4']      ;,'8']
  ;; tau0_arr    = [1.d4]  ;,1.d8]
   
  vth         = 128500.0d0     ; cm/s

  PS_Start, File='plots/sphere_np1e6_T1e2_ndust0.0_DH0.0_tauH1e4-8_vexp0_NoRecoil_iso_MFP.ps',nomatch=1,font=0

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


  mean_arr = dblarr(5)
  dev_arr  = dblarr(5)
  mfp_expected_Arr = dblarr(5)

  
  for i=0,n_elements(tau0_arr)-1 do begin

      myfile = '/Users/tgarel/Rascas_tests/output/sphere_dom/sphere_np1e6_T1e2_ndust0.0_DH0.0_tauH1e'+tau0_arr_st[i]+'_vexp0_NoRecoil_iso_MFP/photons_out.dat'

     print,myfile
     tau0 = tau0_arr[i]

     temp        = (vth / 12.85d5)^2. * 1.d4
     a           = 4.71d-4 * (temp/1.d4)^(-0.5) ;(12.85 / b)
     b           = vth / 1.d5                   ; km/s
     print,temp,vth
    
     read_photons,myfile,vth,x_out,tau0,status,xlast,nb_abs,time,ids

     ;; MFP = 1 / (sigma_0 * nhi)
     clight  =  2.9979250d10 ; cm/s
     nhi     = nhi_arr[i]                ; cm^-3
     sigma_0 = 5.898d-14 * (temp/1.d4)^(-0.5) ; cm^2
     mfp_expected = 1. / (sigma_0 * nhi) ; cm

     distance = time * clight  ; cm

     print,mean(distance),mfp_expected


     mean_arr[i] = mean(distance)
     dev_arr[i]  = STDDEV(distance)
     mfp_expected_Arr[i] = mfp_expected
     
     ii = where(status ne 2,ni)
     print,ni

     ;; plot,ids,alog10(distance),psym=3,xtitle='photon ID',ytitle='log  [D/cm]'
     ;; toto = dindgen(n_elements(time))
     ;; toto[*] = alog10(mfp_expected)
     ;; oplot,ids,toto,linestyle=2,color=254
     ;; legendold,['n!dHI!n = '+nhi_arr_st[i]+' cm!u-3!n','Expected '+greek('lambda')+'!dMFP!n'],box=0,/right,charsize=1.4,textcolors=[0,254],/bottom,spacing=2.2
     
     histo,alog10(distance),1.,20.,2000,h
     plot_io,h.x,h.dn/total(h.dn),psym=10,/ys,yr=[1.d-10,0.1],xr=[1.,14.],/xs,xtitle='log [D/cm]',ytitle='N/N!dtot!n',thick=3
     oplot,[alog10(mfp_expected),alog10(mfp_expected)],[1.d-10,10.],linestyle=2,color=254
     legendold,['n!dHI!n = '+nhi_arr_st[i]+' cm!u-3!n','Expected '+greek('lambda')+'!dMFP!n'],box=0,/left,charsize=1.4,textcolors=[0,254],spacing=2.2
     
     oo = where(nb_abs ne 1,no)
     print,no
     
  endfor
  
  ;legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u5!n',greek('tau')+'!d0!n = 10!u6!n',greek('tau')+'!d0!n = 10!u7!n',greek('tau')+'!d0!n = 10!u8!n'],textcolors=(1+indgen(5))*50,box=0,/right,charsize=1.4

 
  plot,findgen(5),alog10(mean_arr),psym=sym(1),xtitle='log n!dHI!n [cm!u-3!n]',ytitle='log [<D>/cm]',xr=[-0.5,4.5],/xs,yr=[7.,13.]
  errplot,findgen(5),alog10(mean_arr-dev_arr),alog10(mean_arr+dev_arr),thick=3
  legendold,['RASCAS','Expected '+greek('lambda')+'!dMFP!n'],box=0,/right,charsize=1.4,textcolors=[0,254],spacing=2.2

  oplot,findgen(5),alog10(mfp_expected_Arr),color=254,thick=3
  
  PS_End, /PDF, /Delete_PS
  
  
end
