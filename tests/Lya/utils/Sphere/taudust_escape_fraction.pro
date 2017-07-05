pro taudust_escape_fraction

  ;; tau0_arr_st = ['4','5','6','7'];,'8']
  ;; tau0_arr    = [1.d4,1.d5,1.d6,1.d7];,1.d8]

  tau0_arr_st = ['6','6','6','6','6','6','6','6','6']      ;,'8']
  tau0_arr    = [1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6]  ;,1.d8]
   
  vth         = 128500.0d0       ; cm/s

  path = '/Users/tgarel/Rascas_tests/output/sphere_dom/'
  
  filearr = ['sphere_np1e6_T1e2_ndust0.0_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','sphere_np1e6_T1e2_taudust001_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/','sphere_np1e6_T1e2_taudust01_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/','sphere_np1e6_T1e2_taudust1_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/','sphere_np1e6_T1e2_taudust2_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/','sphere_np1e6_T1e2_taudust5_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin0/']

  ;xin_arr = [ 0,1,2,3,4,5,10,15,20]
  taud_arr = [0,0.01,0.1,1,2,5]
  fesc = dblarr(n_elements(filearr))
  xxx  = dblarr(n_elements(filearr))
  
  PS_Start, File='plots/sphere_fesc_temp100K_tau6_taudustBIS.ps',nomatch=1,font=0

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
     
     ii = where(status eq 1,ni)  ;; escape
     qq = where(status eq 2,nq)  ;; killed by dust

     print,ni,nq,double(ni)/double(nphotons)
    
     fesc[i] = double(ni)/double(nphotons)

     xxx[i] = (a*tau0)^(1./3.)*taud_arr[i]*0.54
     
  endfor

  plot,alog10(xxx),alog10(fesc),xtitle='log[(a'+greek('tau')+'!d0!n)!u1/3!n '+greek('tau')+'!ddust!n]',ytitle='log f!desc!n',charthick=5,xr=[-2.,2.],/xs,thick=8,yr=[-3.,0.5],psym=sym(1) ;,/ylog

  x_neuf = dindgen(1000)/999. * 4.-2. ; =  log10[(a*tau0)^(1./3.)*taud_arr[i]]
  tau0fix = 1.d6
  y_neuf = 1.d0 / cosh(0.565*sqrt(3.)/!pi^(5./12.)*sqrt(10.^(x_neuf)))
  oplot,x_neuf,alog10(y_neuf),linestyle=2,thick=6,color=254
  
  legendold,[greek('tau')+'!d0!n = 10!u6!n','T=100 K'],box=0,/left,/center,charsize=1.3,spacing=1.8

  legendold,['RASCAS','Neufeld90'],box=0,/left,/bottom,charsize=1.5,spacing=2.1,colors=[0,254],textcolors=[0,254]

  PS_End, /PDF, /Delete_PS
  
  print,xxx
  
end
