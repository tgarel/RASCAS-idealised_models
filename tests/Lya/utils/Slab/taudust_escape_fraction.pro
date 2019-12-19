pro taudust_escape_fraction

  tau0_arr_st = ['6','6','6','6','6','6','6','6','6','6']           ;,'8']
  tau0_arr    = [1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6]  ;,1.d8]
   
  vth         = 128500.0d0       ; cm/s

  path = '/Users/tgarel/Rascas/output/Lya_tests/slab/'
  
  filearr = ['slab_np1e6_T1e2_ndust0.0_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust00001_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust0001_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust001_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust01_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust05_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust1_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust2_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust5_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust10_DH0.0_tauH1e6_vexp0_NoRecoil_iso/']

  taud_arr = [0,0.0001,0.001,0.01,0.1,0.5,1,2,5,10]
  fesc = dblarr(n_elements(filearr))
  xxx  = dblarr(n_elements(filearr))
  
  PS_Start, File='plots/slab_fesc_temp100K_tau6_taudust_paperVersion.ps',nomatch=1,font=0

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

  for i=1,n_elements(filearr)-1 do begin

     myfile = path+filearr[i]+'photons_out.dat'
     
     print,myfile
     tau0 = tau0_arr[i]
 

     temp        = (vth / 12.85d5)^2. * 1.d4
     a           = 4.71d-4 * (temp/1.d4)^(-0.5) ;(12.85 / b)
     b           = vth / 1.d5                   ; km/s
     print,temp,vth
     
     read_photons_old,myfile,vth,x_out,tau0,status,xlast,nb_abs,time,ids

     nphotons = n_elements(x_out)
     
     ii = where(status eq 1,ni)  ;; escape
     qq = where(status eq 2,nq)  ;; killed by dust

    ; print,ni,nq,double(ni)/double(nphotons)
    
     fesc[i] = double(ni)/double(nphotons)

     tau0fix = 1.d6
     xxx[i] = (a*tau0fix)^(1./3.)*taud_arr[i]*0.54
     print,'a = ',alog10(xxx[i]),alog10(fesc[i]),taud_arr[i],double(ni)
     
  endfor

  plot,alog10(xxx),alog10(fesc),xtitle='log[(a'+greek('tau')+'!dHI!n)!u1/3!n '+greek('tau')+'!ddust!n]',ytitle='log f!desc!n',charthick=5,xr=[-4.,2.],/xs,thick=8,yr=[-4.5,1.5],psym=sym(1),symsize=1.4 ;,/ylog

  tau0fix = 1.d6
  ;x_neuf = alog10((a*tau0fix)^(1./3.)*taud_arr*0.54)
  x_neuf = dindgen(1000)/999. * 8.-5. ; =  log10[(a*tau0)^(1./3.)*taud_arr[i]]
  y_neuf2 = 1.d0 / cosh(sqrt(3.)/0.5/!pi^(5./12.)*sqrt(10.^(x_neuf)))
  ;y_neuf2 = 1.d0 / cosh(sqrt(3.)/0.525/!pi^(5./12.)*sqrt(10.^(x_neuf)))
  oplot,x_neuf,alog10(y_neuf2),linestyle=2,thick=6,color=254
   
  legendold,[greek('tau')+'!dHI!n = 10!u6!n','T=100 K'],box=0,/left,/center,charsize=1.3,spacing=1.8

  legendold,['RASCAS','Neufeld90'],box=0,/left,/bottom,charsize=1.5,spacing=2.1,colors=[0,254],textcolors=[0,254]

 
  PS_End, /PDF, /Delete_PS
  


  
end
