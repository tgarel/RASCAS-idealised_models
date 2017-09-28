pro nbscatt
  
  tau0_arr_st = ['4','5','6','7','8']
  tau0_arr    = [1.d4,1.d5,1.d6,1.d7,1.d8]

  tau0_arr_st = ['6']      ;,'8']
  tau0_arr    = [1.d6]  ;,1.d8]
   
  vth         = 128500.0d0     ; cm/s

  PS_Start, File='plots/slab_np1e6_T1e2_ndust0.0_DH0.0_tauH1e4-8_vexp0_NoRecoil_iso_NBscatt.ps',nomatch=1,font=0

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

     ; myfile = '/Users/tgarel/Rascas/output/Lya_tests/slab/tauH'+tau0_arr_st[i]+'_temp1d2_xin0/photons_out.dat'
      myfile = '/Users/tgarel/Rascas/output/Lya_tests/slab/tauH'+tau0_arr_st[i]+'_temp1d2_xin0/photons_out.dat'

     print,myfile
     tau0 = tau0_arr[i]

     temp        = (vth / 12.85d5)^2. * 1.d4
     a           = 4.71d-4 * (temp/1.d4)^(-0.5) ;(12.85 / b)
     b           = vth / 1.d5                   ; km/s
     print,temp,vth
    
     read_photons,myfile,vth,x_out,tau0,status,xlast,nb_abs,time,ids,k_ext

     ;; MFP = 1 / (sigma_0 * nhi)
     clight  =  2.9979250d10 ; cm/s  
     distance = time * clight  ; cm

     mean_arr[i] = mean(nb_abs)
     dev_arr[i] = stddev(nb_abs)

     histo,alog10(nb_abs),-1.,9.,500,h
     plot_io,h.x,h.dn/total(h.dn),psym=10,/ys,yr=[1.d-7,0.1],xr=[0.,9.],/xs,xtitle='log N!dscatt!n',ytitle='N/N!dtot!n',thick=3

     oo = where(time lt 0.,no)
     time2 = -1. * time[oo]
     print,no
     if no gt 1 then begin
        histo,alog10(time2),min(alog10(time2)),max(alog10(time2)),1000,h
        plot_io,h.x,h.dn/total(h.dn),psym=10,xtitle='-log time',ytitle='N/N!dtot!n',thick=3
     endif
     legendold,[greek('tau')+'!d0!n = 10!u'+strtrim(string(i+4,format='(i1)'),2)+'!n'],textcolors=[0],box=0,/left,charsize=1.7     

     print,'min/max time = ',min(time),max(time)
     print,no,n_elements(x_out)
     
  endfor


  plot,4+indgen(5),alog10(mean_arr),psym=sym(1),/xs,xr=[3.5,8.5],yr=[3.,10.],/ys,ytitle='log <N!dscatt!n>',xtitle='log '+ greek('tau')+'!d0!n',title='Sphere',symsize=1.5
  errplot,4+indgen(5),alog10(mean_arr-dev_arr),alog10(mean_arr+dev_arr)
  ;oplot,4+indgen(4),alog10(1.*tau0_arr),linestyle=2,color=254
 ; oplot,3+indgen(6),alog10(1.612*10.d0^(3+findgen(6))),linestyle=2,color=254

 ; legendold,['RASCAS','<N!dscatt!n> = 1.612'+ greek('tau')+'!d0!n'],textcolors=[0,254],box=0,/left,charsize=1.4,spacing=2.4

 ; legendold,['RASCAS','Expected '+greek('lambda')+'!dMFP!n'],box=0,/right,charsize=1.4,textcolors=[0,254],spacing=2.2

  ;; plot,10.^(4+indgen(4)),(mean_arr),psym=sym(1),/xs,xr=[1.d3,1.d8],yr=[1.d3,1.d8],/ys,/ylog,/xlog
  ;; errplot,10.^(4+indgen(4)),(mean_arr-dev_arr),(mean_arr+dev_arr)
  ;; oplot,10.^(4+indgen(4)),(1.*tau0_arr),linestyle=2,color=254

  ;legendold,['RASCAS','N!dscatt!n = '+ greek('tau')+'!d0!n'],textcolors=[0,254],box=0,/left,charsize=1.4,spacing=2.4
  
  PS_End, /PDF, /Delete_PS
  
  
end
