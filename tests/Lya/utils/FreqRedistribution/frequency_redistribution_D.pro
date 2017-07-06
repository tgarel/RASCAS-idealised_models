pro frequency_redistribution_D,log_plots=log_plots

  xin_arr   = [0,1,2,3,4,5,8]
  temp_arr  = ['1.E+02','1.E+04']
  ;temp_arr  = ['1.E+04']
  temp_arr_bis = temp_arr
  for k=0,n_elements(temp_arr_bis)-1 do begin
     aaa = temp_arr_bis[k]
     temp_arr_bis[k] = strsplit(aaa,/extract,'1.E+0')
  endfor
  
     
  isotropic = 'true'
  recoil    = 'false'

  if not keyword_set(log_plots) then begin
     scale = 'linear'
  endif else begin
     scale = 'log'
  endelse
  
  if (isotropic eq 'true') then begin
     if (recoil eq 'true') then begin
        myplot_title = 'plots/tests_D_FreqRedist_isotropic_Recoil_'+scale
     endif else begin
        myplot_title = 'plots/tests_D_FreqRedist_isotropic_NoRecoil_'+scale
     endelse
  endif else begin
     if (recoil eq 'true') then begin
        myplot_title = 'plots/tests_D_FreqRedist_dipolar_Recoil_'+scale
     endif else begin
       ;; myplot_title = 'plots/tests_D_FreqRedist_dipolar_NoRecoil_'+scale
        myplot_title = 'plots/tests_D_FreqRedist_Rayleigh_NoRecoil_'+scale
     endelse
  endelse

  startpdf,myplot_title
  loadct,39

  mycol = [cgcolor('cg3'),cgcolor('chartreuse'),cgcolor('ygb3'),210,cgcolor('red4'),cgcolor('violet'),cgcolor('gray')]
  xin_arr = [0,1,2,3,4,5,8]

  for p=0,n_elements(temp_arr)-1 do begin

     temp = temp_arr[p]
     
     for j=0,n_elements(xin_arr)-1 do begin

        x_in = xin_arr(j)
        
        nphotons = 0L
        xx  = -1.
        kk1 = -1.
        kk2 = -1.
        kk3 = -1.
        
        if (isotropic eq 'true') then begin
           if (recoil eq 'true') then begin
              outfile = '../../LyaTests_output_data/SingleScatt_D/tests_D_xin'+strtrim(string(x_in,format='(i1)'),2)+'_T'+strtrim(temp,2)+'_Isotropic_Recoil.dat'
           endif else begin
              outfile = '../../LyaTests_output_data/SingleScatt_D/tests_D_xin'+strtrim(string(x_in,format='(i11)'),2)+'_T'+strtrim(temp,2)+'_Isotropic_NoRecoil.dat'
           endelse
        endif else begin
           if (recoil eq 'true') then begin
              outfile = '../../LyaTests_output_data/SingleScatt_D/tests_D_xin'+strtrim(string(x_in,format='(i1)'),2)+'_T'+strtrim(temp,2)+'_Dipolar_Recoil.dat'
           endif else begin
             ; outfile = '../../LyaTests_output_data/SingleScatt_D/tests_D_xin'+strtrim(string(x_in,format='(i1)'),2)+'_T'+strtrim(temp,2)+'_Dipolar_NoRecoil.dat'
              outfile = '../../LyaTests_output_data/SingleScatt_D/tests_D_xin'+strtrim(string(x_in,format='(i1)'),2)+'_T'+strtrim(temp,2)+'_Rayleigh_NoRecoil.dat'
           endelse
        endelse

        ;; print,outfile
        
        openr,11,outfile
        readf,11, nphotons, x_in, T, nu_D
        ;; print,'Nphotons = ',nphotons
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

        mytitle = 'Deuterium - T=10!u'+strtrim(temp_arr_bis[p],2)+'!n K'
        histo,xout,-20.,20.,1600,h
        if not keyword_set(log_plots) then begin
           if j eq 0 then plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[-5.,13.],/nodata,xtitle='x!dout!n',ytitle='R(x!din!n,x!dout!n)',title=mytitle,charsize=1.3,yr=[0.,0.95]
        ;; ytitle='R!dII-A!n(x!din!n,x!dout!n)'
        endif else begin
           if j eq 0 then plot,h.x,h.dn/total(h.dn)/h.dx,psym=10,xr=[-5.,13.],/nodata,xtitle='x!dout!n',ytitle='R(x!din!n,x!dout!n)',title=mytitle,charsize=1.3,yr=[1.d-4,9.9],/ylog
        endelse
        
        oplot,h.x,h.dn/total(h.dn)/h.dx,psym=10,color=mycol[j],thick=8

     endfor


     ;; so far, I have only computed Hummer62 R functions for T=10^4 K
     if temp_arr[p] eq '1.E+04' then begin
        for j=0,n_elements(xin_arr)-1 do begin
           x_in = xin_arr(j)
           if (isotropic eq 'true') then begin
              restore,'hummer_redistribution_D/hummer_xin'+strtrim(string(x_in,format='(i1)'),2)+'_temp10000K_D.sav'
           endif else begin
              restore,'hummer_redistribution_dipolar_D/hummer_DIPOLAR_xin'+strtrim(string(x_in,format='(i1)'),2)+'_temp10000K_D.sav'
           endelse
           oplot,xout,p_xout,linestyle=2,color=0,thick=6  
        endfor
     endif

     if temp_arr[p] eq '1.E+02' then begin
        for j=0,n_elements(xin_arr)-1 do begin
           x_in = xin_arr(j)
           if (isotropic eq 'true') then begin
              restore,'hummer_redistribution_D/hummer_xin'+strtrim(string(x_in,format='(i1)'),2)+'_temp100K_D.sav'
           endif else begin
              restore,'hummer_redistribution_dipolar_D/hummer_DIPOLAR_xin'+strtrim(string(x_in,format='(i1)'),2)+'_temp100K.sav'
           endelse
           oplot,xout,p_xout,linestyle=2,color=0,thick=6  
        endfor
     endif
      
     
     my_legend = strarr(n_elements(xin_arr))
     for i=0,n_elements(xin_arr)-1 do begin
        my_legend[i] = 'x!din!n = '+strtrim(string(xin_arr[i],format='(i1)'),2)
     endfor
     legendold,[my_legend],linsize=0.45,colors=[mycol],textcolors=[mycol],thick=8 ,charthick=7,charsize=1.2,box=0,/top,/right,linestyle=[0,0,0,0,0,0,0]
     
     
     
     if (isotropic eq 'true') then begin
        if (recoil eq 'true') then begin
           legendold,['isotropic case','with recoil'],colors=[0,0],thick=8 ,charthick=7,charsize=1.2,box=0,/top,/left
        endif else begin
           legendold,['isotropic case','no recoil'],colors=[0,0],thick=8 ,charthick=7,charsize=1.2,box=0,/top,/left        
        endelse
     endif else begin
        if (recoil eq 'true') then begin
           legendold,['anisotropic case','with recoil'],colors=[0,0],thick=8 ,charthick=7,charsize=1.2,box=0,/top,/left
        endif else begin
           ;legendold,['anisotropic case','no recoil'],colors=[0,0],thick=8 ,charthick=7,charsize=1.2,box=0,/top,/left
           legendold,['Rayleigh','no recoil'],colors=[0,0],thick=8 ,charthick=7,charsize=1.2,box=0,/top,/left                 
        endelse
     endelse
     
  endfor
  
  stoppdf
  
end 
