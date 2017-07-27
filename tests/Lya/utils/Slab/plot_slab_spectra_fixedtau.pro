pro plot_slab_spectra_fixedtau

  ;; tau0_arr_st = ['4','5','6','7'];,'8']
  ;; tau0_arr    = [1.d4,1.d5,1.d6,1.d7];,1.d8]

  tau0_arr_st = ['6','6','6','6','6','6','6','6','6','6']      ;,'8']
  tau0_arr    = [1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6,1.d6]  ;,1.d8]
   
  vth         = 128500.0d0       ; cm/s

  path = '/Users/tgarel/Rascas_tests/output/slab_dom/'
  
  filearr = ['slab_np1e6_T1e2_ndust0.0_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust00001_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust0001_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust001_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust01_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust05_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust1_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust2_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust5_DH0.0_tauH1e6_vexp0_NoRecoil_iso/','slab_np1e6_T1e2_taudust10_DH0.0_tauH1e6_vexp0_NoRecoil_iso/']

  ;xin_arr = [ 0,1,2,3,4,5,10,15,20]
  taud_arr = [0,0.0001,0.001,0.01,0.1,0.5,1,2,5,10]

  PS_Start, File='plots/slab_spectraDIJK_temp100K_tau6_vary_taudust.ps',nomatch=1,font=0

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
     xxi = 0.d0

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
     plot,h.x,h.dn/h.dx/double(nphotons),xtitle='x',ytitle='J(x,'+greek('tau')+'!d0!n)',charthick=5,xr=[-40.,40.],/xs,thick=8,yr=[0.0001,max(h.dn/h.dx/double(nphotons))*1.3],/nodata ;,/ylog
      
     ;endif

     legendold,[greek('tau')+'!ddust!n = '+string(taud_arr[i],format='(f8.5)'),'f!desc!n = '+string(double(ni)/double(nphotons),format='(f8.5)')],textcolors=[i*25,i*25],box=0,/left,charsize=1.3,spacing=1.9
     legendold,[greek('tau')+'!d0!n = 10!u6!n','T=100 K'],box=0,/right,charsize=1.3,spacing=1.8

     if ni ne 0L then oplot,h.x,h.dn/h.dx/double(nphotons),thick=5,color=i*25,psym=10 ;,linestyle=i

    ;; Exact Eq. 39 from Smith15 - his tau0 is same as us
     j_x = 1. / 4. / sqrt(6.*!pi) / a / tau0 * xx * xx * (1.d0 / (cosh(sqrt(!pi^3./54.) * abs(xx*xx*xx - xxi*xxi*xxi)/ a / tau0)))
     x_peak = 1.066 * (a * tau0)^(1./3.)
     oplot,xx,j_x*4.*!pi,color=0,linestyle=2,thick=3 ;,psym=10
     print,total(j_x*4.*!pi*(xx[1]-xx[0]))
     
  endfor

;  legendold,[greek('tau')+'!d0!n = 10!u6!n','T=100 K'],box=0,/right,charsize=1.3,spacing=1.8
  
  PS_End, /PDF, /Delete_PS
  
  
end
