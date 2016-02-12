plot,r,r,tit='!6initialize font'
loadct,5
col=-1
col2=230
ntimes=10
ntfrac=0.1*indgen(ntimes)
colfrac=0.1*(1+indgen(ntimes))
print,colfrac
;
nt=150
print,'nt=',nt
t1_Gyr=0.8
t2_Gyr=1.2
t1=t1_Gyr/t0_Gyr
t2=t2_Gyr/t0_Gyr
;
plotdirec=strjoin([s1,'plots/'  ,s0])
;
r_kpc=r*r_disk_kpc
ts_rmax_kpc=ts_rmax*r_disk_kpc
ts_delta_r_kpc=ts_delta_r*r_disk_kpc
;
ts_Lam=        ts_h/ts_l
ts_Phi=        ts_om*ts_tau
ts_VVV=        ts_Uz/ts_v
ts_q=         -ts_G/ts_om
ts_D=          9*ts_G*ts_om*ts_h^2/ts_v^2
ts_Dc=        -!pi/32*(!pi^2 +ts_Uz*ts_h/ts_etat)^2
ts_DoverDc=    ts_D/ts_Dc
ts_pkin_rad=  -sqrt(2./!pi/ts_q)/ts_Lam    
ts_pkin_deg=   180./!pi*ts_pkin_rad
ts_Etot=       total(ts_Br^2+ts_Bp^2+ts_Bzmod^2,2)  ;total magnetic energy
ts_Btot=       sqrt(ts_Etot)
;
ts_Br_mkG=     ts_Br*B0_mkG
ts_Bp_mkG=     ts_Bp*B0_mkG
ts_Bzmod_mkG=  ts_Bzmod*B0_mkG
ts_Btot_mkG=   ts_Btot*B0_mkG
ts_alp_k_kms=  ts_alp_k*h0_kpc/t0_kpcskm
ts_alp_m_kms=  ts_alp_m*h0_kpc/t0_kpcskm
ts_alp_kms=    ts_alp  *h0_kpc/t0_kpcskm
ts_p_rad=atan(ts_Br/ts_Bp)
ts_p_deg=180./!pi*ts_p_rad
ts_h_kpc=   ts_h*h0_kpc
ts_l_kpc=   ts_l*h0_kpc
ts_v_kms=   ts_v*h0_kpc/t0_kpcskm
ts_n_cm3=   ts_n*n0_cm3
ts_Uz_kms=  ts_Uz*h0_kpc/t0_kpcskm
ts_om_Gyr=  ts_om/t0_Gyr
ts_G_Gyr=   ts_G/t0_Gyr
ts_t_Gyr=   ts_t*t0_Gyr
ts_Uphi_kms=dblarr(n1+1,nx)
for i=0,nx-1 do begin
  ts_Uphi_kms[*,i]=ts_om[*,i]*r[i]*r_disk_kpc/t0_kpcskm
endfor
ts_etat_kmskpc= 1.d0/3*ts_l_kpc*ts_v_kms
ts_etat_cm2s=   1.d0/3*ts_l*ts_v*etat0_cm2s/1.d26
ts_tau_Gyr=     ts_l/ts_v*t0_Gyr
ts_tau_Myr=     ts_tau_Gyr*1000
ts_Beq_mkG=     sqrt(4*!pi*ts_n)*ts_v*B0_mkG
;
ts_gamma= deriv(ts_t,alog(ts_Btot))
range_gamma=[where(ts_t ge t1 and ts_t le t2)]
gamma_mean=mean(ts_gamma[range_gamma])
ts_gammatau=gamma_mean*ts_tau
;
ts_gamma_Gyr=  ts_gamma/t0_Gyr
gamma_mean_Gyr=gamma_mean/t0_Gyr
;
!p.multi=0
;
;
;1) Bphi, Br vs r
for iplot=0,1 do begin
  if (iplot eq 0) then begin
    set_plot,'ps'
    cd,plotdirec
    xsize=10. & ysize=6.
    ntstr=strjoin(['nt',strtrim(string(nt),2)])
    device,/encapsul,filename=strjoin([s0,s2,'_B_',ntstr,'.eps']),xsize=xsize,ysize=ysize,bits=8,/color
    char=0.5
    th=1.
    !p.charsize=char
  endif else begin
    device,/close
    set_plot,'x'
    !p.charsize=1.5
    xsize=420 & ysize=300 & xpos=0 & ypos=800
    window,0,xsize=3*xsize,ysize=3*ysize,xpos=xpos,ypos=ypos
  endelse
  ;
  ysp=0.5
  xr=[0,r_disk_kpc]
  yr=[(min(ts_Bp_mkG[nt,nxghost:nx-1-nxghost])<min(ts_Br_mkG[nt,nxghost:nx-1-nxghost]))-ysp*(max(ts_Bp_mkG[nt,nxghost:nx-1-nxghost])>max(ts_Br_mkG[nt,nxghost:nx-1-nxghost])),(1.+ysp)*(max(ts_Bp_mkG[nt,nxghost:nx-1-nxghost])>max(ts_Br_mkG[nt,nxghost:nx-1-nxghost]))]
  xticks=3 & xminor=5
  yticks=4 & yminor=5
  ;
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Bp_mkG   [       nt,nxghost:nx-1-nxghost],li=0,th=th,xr=xr,yr=yr,xtit='!8r !6[kpc]',tit='!8B!D!7u!N!6 (__), !8B!Dr!N!6 (_ _), !3|!8B!Dz!N!3|!6 (_._) [mkG]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Bp_mkG   [       nt,nxghost:nx-1-nxghost],li=0,th=th,col=col
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Br_mkG   [       nt,nxghost:nx-1-nxghost],li=2,th=th,col=col
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Bzmod_mkG[       nt,nxghost:nx-1-nxghost],li=3,th=th,col=col
  ;oplot,[ts_rmax_kpc[nt],ts_rmax_kpc[nt]],[-1000,1000],li=0,col=col
  oplot,[ts_rmax_kpc[nt]-2*ts_delta_r_kpc[nt]/2,ts_rmax_kpc[nt]],[-1000,1000],li=0,col=col
  oplot,[ts_rmax_kpc[nt]+2*ts_delta_r_kpc[nt]/2,ts_rmax_kpc[nt]],[-1000,1000],li=0,col=col
  loadct,38
  for i=0,ntimes-1 do begin
    oplot,r_kpc[nxghost:nx-1-nxghost],ts_Bp_mkG   [ntfrac[i]*nt,nxghost:nx-1-nxghost],li=0,th=th,col=colfrac[i]*col2
    oplot,r_kpc[nxghost:nx-1-nxghost],ts_Br_mkG   [ntfrac[i]*nt,nxghost:nx-1-nxghost],li=2,th=th,col=colfrac[i]*col2
    oplot,r_kpc[nxghost:nx-1-nxghost],ts_Bzmod_mkG[ntfrac[i]*nt,nxghost:nx-1-nxghost],li=3,th=th,col=colfrac[i]*col2
    ;oplot,[ts_rmax_kpc[ntfrac[i]*nt],ts_rmax_kpc[ntfrac[i]*nt]],[-10000,10000],li=0,col=ntfrac[i]*col2-1
    oplot,[ts_rmax_kpc[ntfrac[i]*nt]-2*ts_delta_r_kpc[ntfrac[i]*nt]/2,ts_rmax_kpc[ntfrac[i]*nt]],[-1000,1000],li=0,col=colfrac[i]*col2
    oplot,[ts_rmax_kpc[ntfrac[i]*nt]+2*ts_delta_r_kpc[ntfrac[i]*nt]/2,ts_rmax_kpc[ntfrac[i]*nt]],[-1000,1000],li=0,col=colfrac[i]*col2
  endfor
  loadct,5
  oplot,xr,[0,0],li=1,th=th
  ;
  axis,xaxis=0,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,xaxis=1,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,yaxis=0,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
  axis,yaxis=1,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
  ;
  xyouts,0.5*xr[1],0.9*yr[1],'t(Gyr)='          ,col=col
  xyouts,0.6*xr[1],0.9*yr[1],ts_t_Gyr[       nt],col=col
  loadct,38
  for i=0,ntimes-1 do begin
    xyouts,0.5*xr[1],ntfrac[i]*0.9*yr[1],'t(Gyr)='             ,col=colfrac[i]*col2-1
    xyouts,0.6*xr[1],ntfrac[i]*0.9*yr[1],ts_t_Gyr[ntfrac[i]*nt],col=colfrac[i]*col2-1
  endfor
endfor
;
;cd,programdirec
;stop
;
loadct,5
;2) alp_m, alp_k vs r
for iplot=0,1 do begin
  if (iplot eq 0) then begin
    set_plot,'ps'
    cd,plotdirec
    xsize=10. & ysize=6.
    ntstr=strjoin(['nt',strtrim(string(nt),2)])
    device,/encapsul,filename=strjoin([s0,s2,'_alp_',ntstr,'.eps']),xsize=xsize,ysize=ysize;,bits=8;,/color
    char=0.5
    th=1.
    !p.charsize=char
  endif else begin
    device,/close
    set_plot,'x'
    !p.charsize=1.5
    xsize=420 & ysize=300 & xpos=420 & ypos=800
    window,1,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos
  endelse
  ;
  xr=[0,r_disk_kpc]
  yr=[(min(ts_alp_k_kms[nt,nxghost:nx-1-nxghost])<min(ts_alp_m_kms[nt,nxghost:nx-1-nxghost]))-ysp*(max(ts_alp_k_kms[nt,nxghost:nx-1-nxghost])>max(ts_alp_m_kms[nt,nxghost:nx-1-nxghost])),(1.+ysp)*(max(ts_alp_k_kms[nt,nxghost:nx-1-nxghost])>max(ts_alp_m_kms[nt,nxghost:nx-1-nxghost]))]
  ;yr=[(min(ts_alp_kms[nt,nxghost:nx-1-nxghost]))-0.05*(max(ts_alp_kms[nt,nxghost:nx-1-nxghost])),1.05*(max(ts_alp_kms[nt,nxghost:nx-1-nxghost]))]
  ;yr=[0.,10]
  xticks=3 & xminor=5
  yticks=4 & yminor=5
  ;
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_alp_m_kms[nt,nxghost:nx-1-nxghost],th=th,xr=xr,yr=yr,xtit='!8r !6[kpc]',tit='!7a!D!6m!N (__), !7a!D!6k!N (_ _), !7a!6 (_._) [km s!E-1!N]'
  ;plot ,r_kpc[nxghost:nx-1-nxghost],ts_alp_kms[nt,nxghost:nx-1-nxghost]/ts_alp_k_kms[nt,nxghost:nx-1-nxghost],th=th,xr=xr,yr=yr,xtit='!8r !6[kpc]',tit='!7a!D!6m!N (__), !7a!D!6k!N (_ _), !7a!6 (_._) [km s!E-1!N]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_alp_m_kms[nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_alp_k_kms[nt,nxghost:nx-1-nxghost],th=th,col=col,li=2
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_alp_kms  [nt,nxghost:nx-1-nxghost],th=th,col=col,li=3
  oplot,xr,[0,0],li=1,th=th
  ;
  axis,xaxis=0,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,xaxis=1,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,yaxis=0,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
  axis,yaxis=1,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
endfor
;
;3) -p vs r
for iplot=0,1 do begin
  if (iplot eq 0) then begin
    set_plot,'ps'
    cd,plotdirec
    xsize=10. & ysize=6.
    ntstr=strjoin(['nt',strtrim(string(nt),2)])
    device,/encapsul,filename=strjoin([s0,s2,'_p_',ntstr,'.eps']),xsize=xsize,ysize=ysize
    char=0.5
    th=1.
    !p.charsize=char
  endif else begin
    device,/close
    set_plot,'x'
    !p.charsize=1.5
    xsize=420 & ysize=300 & xpos=840 & ypos=800
    window,2,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos
  endelse
  ;
  xr=[0,r_disk_kpc]
  yr=[min(abs(ts_p_deg[nt,nxghost+1:nx-2-nxghost]))-0.03*max(abs(ts_p_deg[nt,nxghost+1:nx-2-nxghost])),1.03*max(abs(ts_p_deg[nt,nxghost+1:nx-2-nxghost]))]
  plot ,r_kpc[nxghost+1:nx-2-nxghost],-ts_p_deg   [nt,nxghost+1:nx-2-nxghost],th=th,xr=xr,yr=yr,xtit='!8r !6[kpc]',tit='!3|!8p!3| !6(__), !3|!8p!D!6kin!N!3| (_ _) !6[!9%!6]',charth=th 
  oplot,r_kpc[nxghost+1:nx-2-nxghost],-ts_p_deg   [nt,nxghost+1:nx-2-nxghost],th=th,col=col
  oplot,r_kpc[nxghost+1:nx-2-nxghost],-ts_pkin_deg[nt,nxghost+1:nx-2-nxghost],th=th,col=col,li=2
  ;
  axis,xaxis=0,xth=th,xticklen=0.03,xtickname=replicate(' ',30)
  axis,xaxis=1,xth=th,xticklen=0.03,xtickname=replicate(' ',30)
  axis,yaxis=0,yth=th              ,ytickname=replicate(' ',30)
  axis,yaxis=1,yth=th              ,ytickname=replicate(' ',30)
endfor
;
;4) B vs t
!p.multi=[0,1,2]
for iplot=0,1 do begin
  if (iplot eq 0) then begin
    set_plot,'ps'
    cd,plotdirec
    xsize=10. & ysize=6.
    ntstr=strjoin(['nt',strtrim(string(nt),2)])
    device,/encapsul,filename=strjoin([s0,s2,'_t_',ntstr,'.eps']),xsize=xsize,ysize=ysize;,bits=8;,/color
    char=0.5
    th=1.
    !p.charsize=char
  endif else begin
    device,/close
    set_plot,'x'
    !p.charsize=1.5
    xsize=420 & ysize=400 & xpos=1260 & ypos=700
    window,3,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos
  endelse
  ;
  xr=[0,max(ts_t_Gyr)]
  plot ,ts_t_Gyr[0:nt],ts_Btot_mkG[0:nt],th=th,xtit='!8t !6[Gyr]',tit='Log!D10!N(!8E!D!6tot!N!E1/2!N) !6[!7l!6G]',charth=th,ylog=1
  oplot,ts_t_Gyr[0:nt],ts_Btot_mkG[0:nt],th=th,col=col
  oplot,[t1_Gyr,t1_Gyr],[1.d-5,1.d5],th=th,li=1
  oplot,[t2_Gyr,t2_Gyr],[1.d-5,1.d5],th=th,li=1
  ;
  axis,xaxis=0,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,xaxis=1,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,yaxis=0,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
  axis,yaxis=1,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
  ;
  plot ,ts_t_Gyr[0:nt],ts_gamma_Gyr[0:nt],th=th,xtit='!8t !6[Gyr]',tit='!7C !6[Gyr!E-1!N]!6',charth=th
  oplot,[0,100000],[0,0],th=th/2,li=2,col=100
  oplot,ts_t_Gyr[0:nt],ts_gamma_Gyr[0:nt],th=th,col=col
  oplot,[t1_Gyr,t1_Gyr],[1.d-5,1.d5],th=th/2,li=1
  oplot,[t2_Gyr,t2_Gyr],[1.d-5,1.d5],th=th/2,li=1
  oplot,[0,100000],[gamma_mean_Gyr,gamma_mean_Gyr],th=th/2,li=1
  ;
  axis,xaxis=0,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,xaxis=1,xth=th,xticklen=0.03,xtickname=replicate(' ',30);,xticks=xticks,xminor=xminor
  axis,yaxis=0,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
  axis,yaxis=1,yth=th              ,ytickname=replicate(' ',30);,yticks=yticks,yminor=yminor
endfor
!p.multi=0
;
;5) h,l,v,n,Uphi,Uz,etat,tau,Beq
for iplot=0,1 do begin
  if (iplot eq 0) then begin
    set_plot,'ps'
    cd,plotdirec
    xsize=10. & ysize=6.
    ntstr=strjoin(['nt',strtrim(string(nt),2)])
    device,/encapsul,filename=strjoin([s0,s2,'_profiles_',ntstr,'.eps']),xsize=xsize,ysize=ysize;,bits=8;,/color
    char=0.5
    th=1.
    !p.charsize=char
  endif else begin
    device,/close
    set_plot,'x'
    !p.charsize=3.
    xsize=1680 & ysize=500 & xpos=0 & ypos=100
    window,4,xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos
  endelse
  ;
  !p.multi=[0,6,3]
  ;a) h
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_h_kpc      [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8h !6[kpc]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_h_kpc      [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;b) l
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_l_kpc      [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8l !6[kpc]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_l_kpc      [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;c) v
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_v_kms      [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8v !6[km s!E-1!N]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_v_kms      [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;d) n
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_n_cm3      [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8n !6[cm!E-3!N]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_n_cm3      [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;e) Uphi
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Uphi_kms   [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8U!D!7u!N !6[km s!E-1!N]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Uphi_kms   [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;f) Uz
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Uz_kms     [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8U!Dz!N !6[km s!E-1!N]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Uz_kms     [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;g) etat
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_etat_cm2s  [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!7g!D!6t!N [10!E26!Ncm!E2!N s]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_etat_cm2s  [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;h) tau
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_tau_Myr    [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!7s !6[Myr]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_tau_Myr    [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;i) Beq
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Beq_mkG    [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8B!D!6eq!N [!7l!6G]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Beq_mkG    [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;j) G
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_G_Gyr      [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8S !6[Gyr!E-1!N]'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_G_Gyr      [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;k) D
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_D          [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8D!6'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_D          [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;l) Dcrit
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Dc         [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8D!D!6c!N'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Dc         [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;m) Lambda
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Lam        [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!7K!6'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Lam        [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;n) Phi
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_Phi        [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!7U!6'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_Phi        [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;o) VVV
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_VVV        [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8V!6'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_VVV        [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  ;
  ;p) D/Dcrit
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_DoverDc    [nt,nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!8D/D!D!6c!N!6'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_DoverDc    [nt,nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  oplot,[0,1000],[1,1],li=1
  ;
  ;q) gamma*tau
  plot ,r_kpc[nxghost:nx-1-nxghost],ts_gammatau   [nxghost:nx-1-nxghost],th=th,xtit='!8r !6[kpc]',tit='!7C!D!6mean!N!7s!6'
  oplot,r_kpc[nxghost:nx-1-nxghost],ts_gammatau   [nxghost:nx-1-nxghost],th=th,col=col
  oplot,[r_sol_kpc,r_sol_kpc],[-1000,1000],li=1
  oplot,[0,1000],[1,1],li=1
  ;
endfor
cd,programdirec
;
end
