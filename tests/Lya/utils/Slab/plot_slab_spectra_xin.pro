pro plot_slab_spectra_xin

  ;; tau0_arr_st = ['4','5','6','7'];,'8']
  ;; tau0_arr    = [1.d4,1.d5,1.d6,1.d7];,1.d8]

 
   
  vth         = 128500.0d0       ; cm/s

  PS_Start, File='plots/slab_spectraNEUF_temp100K_tau6_vary_xin.ps',nomatch=1,font=0

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

  xin    = [0., 1., 2., 3., 4., 5., 10., 20.] ; input frequency
  xin_st = ['0', '1', '2', '3', '4', '5', '10', '20'] ; input frequency

  tau0_arr_st    = strarr(n_elements(xin))
  tau0_arr_st[*] = '6'   
  tau0_arr       = fltarr(n_elements(xin))
  tau0_arr[*]    = 1.d6
  
  for i=0,n_elements(xin)-1 do begin

     xxi = xin[i]

     if i eq 0 then begin
        myfile = '/Users/tgarel/Rascas_tests/output/slab_dom/slab_np1e6_T1e2_ndust0.0_DH0.0_tauH1e6_vexp0_NoRecoil_iso/photons_out.dat'
     endif else begin
        myfile = '/Users/tgarel/Rascas_tests/output/slab_dom/slab_np1e6_T1e2_ndust0.0_DH0.0_tauH1e6_vexp0_NoRecoil_iso_xin'+xin_st[i]+'/photons_out.dat'        
     endelse
     

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
     histo,x_out,-250.,250.,500,h

    ; if i eq 0 then begin
        plot,h.x,h.dn/h.dx/double(nphotons),xtitle='x',ytitle='J(x,'+greek('tau')+'!d0!n)',charthick=5,xr=[-60,60.],/xs,thick=8,yr=[0.0001,max(h.dn/h.dx/double(nphotons))*1.45],/nodata ;,/ylog
       ; legendold,['R!dsphere!n = 0.3 L!dbox!n','R!dsphere!n within 1 cell'],textcolors=col,box=0,/left,charsize=1.3,spacing=2.1
   ;  endif

     oplot,h.x,h.dn/h.dx/double(nphotons),thick=5,color=i*30,psym=10 ;,linestyle=i

     ;; Exact Eq. 39 from Smith15 - his tau0 is same as us
     j_x = 1. / 4. / sqrt(6.*!pi) / a / tau0 * xx * xx * (1.d0 / (cosh(sqrt(!pi^3./54.) * abs(xx*xx*xx - xxi*xxi*xxi)/ a / tau0)))

     x_peak = 1.066 * (a * tau0)^(1./3.)
     oplot,xx,j_x*4.*!pi,color=0,linestyle=2,thick=3 ;,psym=10
     ;; plots,[-1.*x_peak,-1.*x_peak],[0.,1.*max(j_x*2.*!pi)],linestyle=2.,color=254
     ;; plots,[1.*x_peak,1.*x_peak],[0.,1.*max(j_x*2.*!pi)],linestyle=2,color=254

     print,total(j_x*4.*!pi*(xx[1]-xx[0]))
     
     legendold,['x!din!n = '+xin_st[i]],textcolors=30*i,box=0,/left,charsize=1.6
     legendold,['Hydrogen','T=10!u2!n K'],/left,/center,charthick=7,charsize=1.2,box=0,spacing=1.9
     legendold,['RASCAS','Neufeld90'],box=0,/right,charsize=1.2,spacing=2.1,colors=[30*i,0],textcolors=[30*i,0],linestyle=[0,2],linsize=0.3
  endfor

  ;legendold,['Nxbin = 2000','Nxbin=500'],textcolors=col,box=0,/right,charsize=1.4

 ; legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u7!n'],textcolors=col,box=0,/right,charsize=1.4

  
 ; legendold,[greek('tau')+'!d0!n = 10!u4!n',greek('tau')+'!d0!n = 10!u5!n',greek('tau')+'!d0!n = 10!u6!n',greek('tau')+'!d0!n = 10!u7!n'],textcolors=col,box=0,/right,charsize=1.4
  
  PS_End, /PDF, /Delete_PS
  
  
end
