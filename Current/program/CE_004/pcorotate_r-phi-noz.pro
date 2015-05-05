@colorbar
set_plot,'x'
xsizeps=5.
xsize_screen=220.
ysize_screen=xsize_screen
xpos_screen=1200.
ypos_screen=500.
charfac=1.5*xsizeps/12
xchar=1.5*xsizeps/12
cth=1.5
Corotate=0
Windup=1
Plot_part=0
diamond=0	;whether to plot diamonds at local maxima of delta
Arrow=0
ndig=6	;number of characters in time in Gyr (including decimal point)
window,xsize=xsize_screen,ysize=ysize_screen,xpos=xpos_screen,ypos=ypos_screen
plot,r(*,0),r(*,0),tit='!6Initialize font!6'
col0=-1
col1=255
col2=95
col3=250
col4=20
col5=80
!p.multi=[0,5,4]
plotdirec=strjoin([s1,'/fortran_pde/2D/telegraph/r_phi_noz/plots/',s0])
cd,plotdirec
xr=[0,0.5]*r_disk_kpc
yr=[0,2*!pi]*180./!pi
;rc=0.4
;nrc=rc*nx
nrsamp=0.4*nx
;nphisamp=where(alp_k[nrc,*] eq max(alp_k[nrc,*]))
nphisamp=0
;
B=sqrt(ts_Br^2+ts_Bp^2)
for j=0,ny-1 do begin
for i=0,nx-1 do begin
  delta[*,i,j]=(ts_Bp[*,i,j]-ts_Bzero[*,1,i])/ts_Bzero[*,1,i]
endfor
endfor
delta[*,nxghost   ,*]=0.
delta[*,nx-1-nxghost,*]=0.
;
pB=-180./!pi*atan(ts_Br/ts_Bp)
pB[*,nxghost   ,*]=0.
pB[*,nx-1-nxghost,*]=0.
Bcol=5;16;11
if (nvar eq 2 or nvar eq 4) then begin
  Bmax=max(B/Beq0)
endif else begin
  Bmax=0.85;0.5
  Bmaxlog=0.20;0.14
endelse
deltamax=0.6;0.6
alpmax=0.9*max(alp_k)
alpnormmax=0.9;0.60;0.75
Uzmax=0.96
alppartnormmax=0.95
Uzpartnormmax=0.95
pBmax= 25.
pBmin= 0.
;divcol=6
ncol=24;42	;should be divisible by divcol if keyword bottom used in colorbar
VBlog=Bmaxlog*findgen(ncol+1)/ncol
VB=Bmax*findgen(ncol+1)/ncol
;VB=Bmax*(findgen(ncol+ncol/divcol+1)-float(ncol/divcol))/ncol		;start at color_id /= 0
;bottom=fix(float(ncol/divcol)/(ncol+ncol/divcol)*256)
VpB=2*pBmax*(findgen(ncol+1)-ncol/2)/ncol
VpB=(pBmax -pBmin)*findgen(ncol+1)/ncol +pBmin
Vdelta=deltamax*(findgen(ncol+1)-ncol/2)*2./ncol
Vdeltalog=deltamax*(findgen(ncol+1)-ncol/2)*2./ncol
;
rnew=  dblarr(nx,ny+1)
phinew=dblarr(nx,ny+1)
rnew  [*,0:ny-1]= r
phinew[*,0:ny-1]= phi
rnew  [*,ny]=   rnew[*,0]	;to avoid leaving out a wedge
phinew[*,ny]= 2*!pi		;to avoid leaving out a wedge
;
Brnew=   dblarr(nx,ny+1)
Bpnew=   dblarr(nx,ny+1)
Bznew=   dblarr(nx,ny+1)
Bnew=    dblarr(nx,ny+1)
deltanew=dblarr(nx,ny+1)
pBnew=   dblarr(nx,ny+1)
r2=reform(r[*,nphisamp],nx)
; 
;xlogr=[0,1]
xlogr=[0.60,0.95]
box=9./r_disk_kpc*[-1,1]
;
r_disk=1.
rclong_kpc =dblarr(nspiral,1000001)
xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
;
r2=reform(  r[   *,nphisamp],nx)
;
ts_meanB=fltarr(n1+1)
for i=0,n1 do begin
  ts_meanB[i]=total(r[0:nx/2,*,*]*B[i,0:nx/2,*,*])/total(r[0:nx/2,*,*])	;multiply by r since grid is non-uniform and cell-size is proportional to r
endfor
;
;
!p.multi=0
Gam=    deriv(ts_t,alog(+ts_meanB))	;m=0 tau=tau1	growth rate of B_phi(m=0) averaged over all radii
Gam_dec=deriv(ts_t,alog(-ts_meanB))	;		decay  rate of B_phi(m=0) averaged over all radii (if Gamma<0)
good=where(ts_t gt t1 and ts_t lt t2)
Gamma=    mean(Gam(good))			;		time-averaged growth rate
Gamma_dec=mean(Gam_dec(good))			;		time-averaged  decay rate
;
window,1,xsize=800,ysize=500
plot,ts_t[0:100],ts_meanB[0:100]
;cd,programdirec
;stop
;
;
asp=1.
nvert=1
!p.multi=0
sym1=1
sym2=4
symcol=col0
symth=12.*xsizeps/12;3.5
symsize=8.*xsizeps/12;2
;
ntmin=410
ntmax=410
;
for nt=ntmin,ntmax do begin
  print,'nt=',nt
  ;
  if (Corotate ne 1) then begin
    for isp=0,nspiral-1 do begin
      rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
    endfor
  endif
  xcoord=dblarr(nx,ny+1)
  ycoord=dblarr(nx,ny+1)
  ;
  Bnew=    dblarr(nx,ny+1)
  deltanew=dblarr(nx,ny+1)
  pBnew=   dblarr(nx,ny+1)
  Brnew   [*,0:ny-1]=ts_Br[nt,*,*]
  Bpnew   [*,0:ny-1]=ts_Bp[nt,*,*]
  Bznew   [*,0:ny-1]=ts_Bzmod[nt,*,*]
  Bnew    [*,0:ny-1]=    B[nt,*,*]
  deltanew[*,0:ny-1]=delta[nt,*,*]
  pBnew   [*,0:ny-1]=   pB[nt,*,*]
  Brnew   [*,ny]=    Brnew[*,0]	;TO AVOID LEAVING OUT A WEDGE
  Bpnew   [*,ny]=    Bpnew[*,0]	;TO AVOID LEAVING OUT A WEDGE
  Bznew   [*,ny]=    Bznew[*,0]	;TO AVOID LEAVING OUT A WEDGE
  Bnew    [*,ny]=     Bnew[*,0]	;TO AVOID LEAVING OUT A WEDGE
  deltanew[*,ny]= deltanew[*,0]	;TO AVOID LEAVING OUT A WEDGE
  pBnew   [*,ny]=    pBnew[*,0]	;TO AVOID LEAVING OUT A WEDGE
  phinewcorotate=  dblarr(nx,ny+1)
  Bnewcorotate=    dblarr(nx,ny+1)
  Brnewcorotate=   dblarr(nx,ny+1)
  Bpnewcorotate=   dblarr(nx,ny+1)
  Bznewcorotate=   dblarr(nx,ny+1)
  deltanewcorotate=dblarr(nx,ny+1)
  pBnewcorotate=   dblarr(nx,ny+1)
  B2=reform(B[nt,*,nphisamp],nx)
  pB2=reform(pB[nt,*,nphisamp],nx)
  r2=reform(r[*,nphisamp],nx)
  Bnew    [*,0:ny-1]=    B[nt,*,*]
  deltanew[*,0:ny-1]=delta[nt,*,*]
  pBnew   [*,0:ny-1]=   pB[nt,*,*]
  Bnew    [*,ny]=    B[nt,*,0]
  deltanew[*,ny]=delta[nt,*,0]
  pBnew   [*,ny]=   pB[nt,*,0]
  ;
  if (Corotate eq 1) then begin
    phinewcorotate=phinew mod (2*!pi) +2*!pi ;Override the transformation to the frame which corotates with the spiral pattern (since there is no such frame)
  endif else begin
    phinewcorotate= (phinew -ts_Omega[nt,0]*nt*0.025/td_Gyr) mod (2*!pi) +2*!pi
  endelse
  phinewcorotate[*,ny]= phinewcorotate[0] +2*!pi
  wherephizero= where(phinewcorotate[0,*] eq min(phinewcorotate[0,*]))
  wherephizero= wherephizero[0]
  ;
  for i=0,nx-1 do begin
    phinewcorotate  [i,0:ny-1]= shift(phinewcorotate[i,0:ny-1],-wherephizero)
    Bnewcorotate    [i,0:ny-1]= shift(Bnew          [i,0:ny-1],-wherephizero)
    Brnewcorotate   [i,0:ny-1]= shift(Brnew         [i,0:ny-1],-wherephizero)
    Bpnewcorotate   [i,0:ny-1]= shift(Bpnew         [i,0:ny-1],-wherephizero)
    Bznewcorotate   [i,0:ny-1]= shift(Bznew         [i,0:ny-1],-wherephizero)
    deltanewcorotate[i,0:ny-1]= shift(deltanew      [i,0:ny-1],-wherephizero)
    pBnewcorotate   [i,0:ny-1]= shift(pBnew         [i,0:ny-1],-wherephizero)
  endfor
  phinewcorotate  [*,ny]= phinewcorotate  [*,0] +2*!pi
  Bnewcorotate    [*,ny]= Bnewcorotate    [*,0]
  Brnewcorotate   [*,ny]= Brnewcorotate   [*,0]
  Bpnewcorotate   [*,ny]= Bpnewcorotate   [*,0]
  Bznewcorotate   [*,ny]= Bznewcorotate   [*,0]
  deltanewcorotate[*,ny]= deltanewcorotate[*,0]
  pBnewcorotate   [*,ny]= pBnewcorotate   [*,0]
  ;
  xcoord= rnew*cos(phinewcorotate)
  ycoord= rnew*sin(phinewcorotate)
  ; 
  alp_knew=     dblarr(nx,ny+1)
  alp_knew_max= dblarr(nx     )
  alp_knew_norm=dblarr(nx,ny+1)
  ;
  Uznew=     dblarr(nx,ny+1)
  Uznew_max= dblarr(nx     )
  Uznew_norm=dblarr(nx,ny+1)
  Uznewdelta=dblarr(nx,ny+1)
  ;
  alp_k_partnew=     dblarr(nx,ny+1,nspiral)
  alp_k_partnew_max= dblarr(nx     ,nspiral)
  alp_k_partnew_norm=dblarr(nx,ny+1,nspiral)
  ;
  Uz_partnew=     dblarr(nx,ny+1,nspiral)
  Uz_partnew_max= dblarr(nx     ,nspiral)
  Uz_partnew_norm=dblarr(nx,ny+1,nspiral)
  ;
  alp_knewcorotate=     dblarr(nx,ny+1)
  alp_knew_maxcorotate= dblarr(nx     )
  alp_knew_normcorotate=dblarr(nx,ny+1)
  ;
  Uznewcorotate=     dblarr(nx,ny+1)
  Uznew_maxcorotate= dblarr(nx     )
  Uznew_normcorotate=dblarr(nx,ny+1)
  ;
  alp_k_partnewcorotate=     dblarr(nx,ny+1,nspiral)
  alp_k_partnew_maxcorotate= dblarr(nx     ,nspiral)
  alp_k_partnew_normcorotate=dblarr(nx,ny+1,nspiral)
  ;
  Uz_partnewcorotate=     dblarr(nx,ny+1,nspiral)
  Uz_partnew_maxcorotate= dblarr(nx     ,nspiral)
  Uz_partnew_normcorotate=dblarr(nx,ny+1,nspiral)
  ;
  alp_knew     [*,0:ny-1  ]=  ts_alp_k     [nt,*,*  ]
  alp_knew     [*,ny      ]=  alp_knew     [   *,0  ]
  alp_k_partnew[*,0:ny-1,*]=  ts_alp_k_part[nt,*,*,*]
  alp_k_partnew[*,ny    ,*]=  alp_k_partnew[   *,0,*]
  ;
  Uznew     [*,0:ny-1]=  ts_Uz[nt,*,*]
  Uznew     [*,ny    ]=  Uznew[*,0]
  Uz_partnew[*,0:ny-1,*]=  ts_Uz_part[nt,*,*,*]
  Uz_partnew[*,ny    ,*]=  Uz_partnew[   *,0,*]
  ;
  for i=0,nx-1 do begin
    alp_knew_max [i  ]=  max(alp_knew[i,*])
    alp_knew_norm[i,*]= alp_knew[i,*]/alp_knew_max[i]
    for isp=0,nspiral-1 do begin
      alp_k_partnew_max [i  ,isp]= max(alp_k_partnew[i,*,isp])
      alp_k_partnew_norm[i,*,isp]= alp_k_partnew[i,*,isp]/alp_k_partnew_max[i,isp]
    endfor
  endfor
  ;
  for i=0,nx-1 do begin
    Uznew_max [i  ]=  max(Uznew[i,*])
    Uznew_norm[i,*]= Uznew[i,*]/Uznew_max[i]
    Uznewdelta[i,*]= (Uznew[i,*] -mean(Uznew[i,*]))/mean(Uznew[i,*])
    if (i ge 1 and max(Uznewdelta[i,*]) ge 0.1 and max(Uznewdelta[i-1,*]) lt 0.1) then begin 
      ixU=i
    endif
    for isp=0,nspiral-1 do begin
      Uz_partnew_max [i  ,isp]= max(Uz_partnew[i,*,isp])
      Uz_partnew_norm[i,*,isp]= Uz_partnew[i,*,isp]/Uz_partnew_max[i,isp]
    endfor
  endfor
  ;
  alp_knew_maxcorotate= alp_knew_max
  alp_k_partnew_maxcorotate= alp_k_partnew_max
  for i=0,nx-1 do begin
    alp_knewcorotate          [i,0:ny-1]= shift(alp_knew          [i,0:ny-1],-wherephizero)
    alp_knew_normcorotate     [i,0:ny-1]= shift(alp_knew_norm     [i,0:ny-1],-wherephizero)
    for isp=0,nspiral-1 do begin
      alp_k_partnewcorotate     [i,0:ny-1,isp]= shift(alp_k_partnew     [i,0:ny-1,isp],-wherephizero)  
      alp_k_partnew_normcorotate[i,0:ny-1,isp]= shift(alp_k_partnew_norm[i,0:ny-1,isp],-wherephizero)
    endfor
  endfor
  alp_knewcorotate          [*,ny  ]= alp_knewcorotate          [*,0  ]
  alp_knew_normcorotate     [*,ny  ]= alp_knew_normcorotate     [*,0  ]
  alp_k_partnewcorotate     [*,ny,*]= alp_k_partnewcorotate     [*,0,*]
  alp_k_partnew_normcorotate[*,ny,*]= alp_k_partnew_normcorotate[*,0,*]
  ;
  Uznew_maxcorotate= Uznew_max
  Uz_partnew_maxcorotate= Uz_partnew_max
  for i=0,nx-1 do begin
    Uznewcorotate          [i,0:ny-1]= shift(Uznew          [i,0:ny-1],-wherephizero)
    Uznew_normcorotate     [i,0:ny-1]= shift(Uznew_norm     [i,0:ny-1],-wherephizero)
    for isp=0,nspiral-1 do begin
      Uz_partnewcorotate     [i,0:ny-1,isp]= shift(Uz_partnew     [i,0:ny-1,isp],-wherephizero)  
      Uz_partnew_normcorotate[i,0:ny-1,isp]= shift(Uz_partnew_norm[i,0:ny-1,isp],-wherephizero)
    endfor
  endfor
  Uznewcorotate          [*,ny  ]= Uznewcorotate          [*,0  ]
  Uznew_normcorotate     [*,ny  ]= Uznew_normcorotate     [*,0  ]
  Uz_partnewcorotate     [*,ny,*]= Uz_partnewcorotate     [*,0,*]
  Uz_partnew_normcorotate[*,ny,*]= Uz_partnew_normcorotate[*,0,*]
  ;
  ;
  deltawhere=where(deltanew eq max(deltanew))
  deltawherecorotate=where(deltanewcorotate gt 0.99999*max(deltanewcorotate))
  deltawherecorotate_half1=where(deltanewcorotate eq max(deltanewcorotate[*,0     :ny/2]))
  deltawherecorotate_half2=where(deltanewcorotate eq max(deltanewcorotate[*,ny/2+1:ny  ]))
  for jw=0,ny do begin
    wherex=where(deltanew[*,jw] eq max(deltanew))
    sizewherex=size(wherex)
      if (wherex[0] ne -1) then begin
        deltawherex=wherex
      endif
  endfor
  if (nspiral eq 1 and n_arm[0] eq 3) then begin
    deltawherex=[deltawherex,deltawherex,deltawherex]
  endif else begin
    if (nspiral eq 1 and n_arm[0] eq 2) then begin
      deltawherex=[deltawherex,deltawherex]
    endif
  endelse
  ;    
  for iw=0,nx-1 do begin
    wherey=where(deltanew[iw,*] eq max(deltanew))
    whereycorotate=where(deltanewcorotate[iw,*] gt 0.99999*max(deltanewcorotate))
    whereycorotate_half1=where(deltanewcorotate[iw,*] eq max(deltanewcorotate[*,0     :ny/2]))
    whereycorotate_half2=where(deltanewcorotate[iw,*] eq max(deltanewcorotate[*,ny/2+1:ny  ]))
    if (wherey[0] ne -1) then begin
      deltawherey=wherey
    endif
    if (whereycorotate[0] ne -1) then begin
      deltawhereycorotate=whereycorotate
    endif
    if (whereycorotate_half1[0] ne -1) then begin
      deltawhereycorotate_half1=whereycorotate_half1
    endif
    if (whereycorotate_half2[0] ne -1) then begin
      deltawhereycorotate_half2=whereycorotate_half2
    endif
  endfor
    ;
;  if (nspiral gt 1) then begin
    for jw=0,ny/2 do begin
      wherexcorotate_half1=where(deltanewcorotate[*,jw] eq max(deltanewcorotate[*,0:ny/2   ]))
      if (wherexcorotate_half1[0] ne -1) then begin
        deltawherexcorotate_half1=wherexcorotate_half1
      endif
    endfor
    for jw=ny/2+1,ny do begin
      wherexcorotate_half2=where(deltanewcorotate[*,jw] eq max(deltanewcorotate[*,ny/2+1:ny]))
      if (wherexcorotate_half2[0] ne -1) then begin
        deltawherexcorotate_half2=wherexcorotate_half2
      endif
    endfor
;  endif
  ;
  ;  plot delta in separate file
  ;
  ; O) Contour delta (x-y)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_deltacartcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
  ;
    loadct,5,/silent
    contour,deltanewcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Vdelta,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,/fill,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    ;contour,deltanewcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=0.1,/overplot,col=col1,c_th=4;,/cell_fill	;overplot orig alpha contour (max of matter spiral)
    ;
    if (Plot_part eq 1) then begin
      if (n_arm[0] eq 2) then begin
        loadct,38,/silent
        contour,alp_k_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        loadct,3,/silent
      endif else begin
        contour,Uz_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endelse
      if (nspiral gt 1) then begin
        contour,Uz_partnew_normcorotate[*,*,1],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endif
    endif
    loadct,5,/silent
    contour,alp_knew_normcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alpnormmax*max(alp_knew_norm),/overplot,col=col3,c_li=0,c_th=2*cth,xstyle=1,ystyle=1
    ;contour,Uznew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,Uznew_normcorotate[ixU:nx-1-nxghost,*],xcoord[ixU:nx-1-nxghost,*]*r_disk_kpc,ycoord[ixU:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm[ixU:nx-1-nxghost,*]),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      endif
    endif else begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
      endif
    endelse
    if (Corotate ne 1) then begin
      rclong_kpc =dblarr(nspiral,1000001)
      for isp=0,nspiral-1 do begin
        rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
      endfor
      xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
      for isp=0,nspiral-1 do begin
        oplot,xcoordlong_kpc,+sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
        oplot,xcoordlong_kpc,-sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
      endfor
    endif
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    ;xyouts,0.40,0.85,/normal,'deltamax='                               ,col=col1,charth=cth
    ;xyouts,0.55,0.85,/normal,deltamax                                     ,col=col1,charth=cth
    ;xyouts,0.47,0.78,/normal,s0                                       ,col=col1,charth=cth
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=-1
    ;
  endfor
  ;
  ; I) Contour delta (phi-log10(r))
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_deltacorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    ;
  ;
    loadct,5,/silent
    contour,transpose(deltanewcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=ncol,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',levels=Vdeltalog,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],/cell_fill,charth=cth,col=col0,xstyle=1,ystyle=1
    ;
    if (ts_t[nt] gt ti_spiral[0] and ts_t[nt] lt tf_spiral[0]) then begin
      if (Corotate ne 1 and Windup ne 1) then begin
        for isp=0,nspiral-1 do begin
          oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(r_disk_kpc*ts_rc[nt,isp]),li=1,col=col1,th=4./3*cth
        endfor
      endif
      ;contour,transpose(alp_knew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpnormmax*max(alp_knew_norm),/overplot,c_th=6,col=col0	;overplot orig alpha contour (max of matter spiral)
      contour,transpose(Uznew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=0.50*max(Uznew_norm),/overplot,c_th=2*cth,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      if (Plot_part eq 1) then begin
        if (n_arm[0] eq 2) then begin
          loadct,38,/silent
          contour,transpose(Uz_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(Uznew_norm),/overplot,c_th=cth,c_li=5,col=col2,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        endif else begin
          loadct,3,/silent
          contour,transpose(Uz_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=3,col=col5,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        endelse
        if (nspiral gt 1) then begin
          loadct,3,/silent
          contour,transpose(Uz_partnew_normcorotate[*,*,1]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(Uznew_norm),/overplot,c_th=cth,c_li=3,col=col5,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        endif
      endif
    endif
    !p.charsize=xchar
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    ;if (deltawherex[0] ne -1 and deltawhereycorotate[0] ne -1 and sizewherex[1] lt 4) then begin
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
      endif
    endelse
    ;endif
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
  endfor
  ;
  ; II) Contour B/B0 (phi-log10(r))
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_Bcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
  ;
    loadct,5,/silent
    contour,transpose(Bnewcorotate)/Beq0,transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=100,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',/cell_fill,levels=VBlog,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth
    ;
    if (ts_t[nt] gt ti_spiral[0] and ts_t[nt] lt tf_spiral[0]) then begin
      if (Corotate ne 1 and Windup ne 1) then begin
        for isp=0,nspiral-1 do begin
          oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(r_disk_kpc*ts_rc[nt,isp]),li=1,col=col1,th=4./3*cth
        endfor
      endif
      contour,transpose(Uznew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=Uzmax*max(Uznew_norm),/overplot,c_th=2*cth,col=col0,xstyle=1,ystyle=1,/fill	;overplot orig alpha contour (max of matter spiral)
      ;if (n_arm[0] eq 2) then begin
      ;  loadct,38,/silent
      ;  contour,transpose(alp_k_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=3,c_li=5,col=col2	;overplot orig alpha contour (max of matter spiral)
      ;endif else begin
      ;  loadct,3,/silent
      ;  contour,transpose(alp_k_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=3,c_li=3,col=col5	;overplot orig alpha contour (max of matter spiral)
      ;endelse
      ;if (nspiral gt 1) then begin
      ;  loadct,3,/silent
      ;  contour,transpose(alp_k_partnew_normcorotate[*,*,1]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=3,c_li=3,col=col5	;overplot orig alpha contour (max of matter spiral)
      ;endif
    endif
    !p.charsize=xchar
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    ;if (deltawherex[0] ne -1 and deltawherey[0] ne -1 and sizewherex[1] lt 4) then begin
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
      endif
    endelse
    ;endif
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
  endfor
  ;
  ; III) Fourier mode plot
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_fourier_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    E0=(Bzero[  0,*]^2+Bzero[  1,*]^2)/8/!pi
    E1=(Bmmax[0,0,*]^2+Bmmax[0,1,*]^2)/8/!pi/2	;Note extra factor of 2 in denominator compared to m=0 E0 due to <cos^2(phi)>=1/2
    E2=(Bmmax[1,0,*]^2+Bmmax[1,1,*]^2)/8/!pi/2
    E3=(Bmmax[2,0,*]^2+Bmmax[2,1,*]^2)/8/!pi/2
    E4=(Bmmax[3,0,*]^2+Bmmax[3,1,*]^2)/8/!pi/2
    E5=(Bmmax[4,0,*]^2+Bmmax[4,1,*]^2)/8/!pi/2
    E6=(Bmmax[5,0,*]^2+Bmmax[5,1,*]^2)/8/!pi/2
    ;
    loadct,39,/silent
    xrfourier=[0.20,0.50]*r_disk_kpc
    yrfourier=[0,0.55]
    plot,r[*,0]*r_disk_kpc,sqrt(E1/E0),xr=xrfourier,yr=yrfourier,th=8./3*cth,li=0,xtit='!8r!6 [kpc]',ytit='!8B!E(m)!N/B!E(0)!N!6',charth=cth,xminor=4,yminor=4,pos=[0.20*asp,0.15,0.95*asp,0.90]
    oplot,r[*,0]*r_disk_kpc,sqrt(E2/E0),th=8./3*cth,li=3,col=col4+10
    oplot,r[*,0]*r_disk_kpc,sqrt(E3/E0),th=8./3*cth,li=2,col=col4+210
    oplot,r[*,0]*r_disk_kpc,sqrt(E4/E0),th=8./3*cth,li=5,col=col4+60
    oplot,r[*,0]*r_disk_kpc,sqrt(E5/E0),th=8./3*cth,li=0,col=col4+110
    oplot,r[*,0]*r_disk_kpc,sqrt(E6/E0),th=8./3*cth,li=4,col=col4+190
    ;oplot,r[*,0]*r_disk_kpc,sqrt(BB1r(time_mod,*)^2+BB1p(time_mod,*)^2)/sqrt(BB0r(time_mod,*)^2+BB0p(time_mod,*)^2),th=8,li=0,col=col0
    ;
    axis,xaxis=0,xtickname=replicate(' ',30),xth=5./3*cth,color=-1,xminor=4
    axis,xaxis=1,xtickname=replicate(' ',30),xth=5./3*cth,color=-1,xminor=4
    axis,yaxis=0,ytickname=replicate(' ',30),yth=5./3*cth,color=-1,yminor=4
    axis,yaxis=1,ytickname=replicate(' ',30),yth=5./3*cth,color=-1,yminor=4
    ;
    for isp=0,nspiral-1 do begin
      oplot,r[*,0]*0+r_disk_kpc*ts_rc[nt,isp],r[*,0]*r_disk_kpc,th=4./3*cth,li=1,col=col0
    endfor
    !p.charsize=xchar
    ;
    xyouts,0.48,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    loadct,5,/silent
    ;
  endfor
  ;
  ;  ARROW PROPERTIES 
  n_in=0.20;0.999	;Sets fraction of radius below which to use nsparse_phi_in (need not be an integer)
  n_mid=0.999;0.999	;Sets fraction of radius below which to use n_sparse_phi_mid and above which to use n_sparse_phi_out (need not be an integer)
  nsparse_r=12;4;12	;Should divide nx; The larger this number, the smaller the number of arrows in r: N_r= nx/nsparse_r-1
  nsparse_phi_in=10;2 	;Should divide ny; The larger this number, the smaller the number of arrows in phi: N_phi= ny/Nsparse_phi
  nsparse_phi_mid=10 	;Should divide ny; The larger this number, the smaller the number of arrows in phi: N_phi= ny/Nsparse_phi
  nsparse_phi_out=1	;Should divide ny; The larger this number, the smaller the number of arrows in phi: N_phi= ny/Nsparse_phi
  N_r=nx/nsparse_r-1	;Number of arrows in r
  N_rmax=N_r		;Cut off outer part of disc for plot
  N_phi_in=ny/nsparse_phi_in	;Number of arrows in phi
  N_phi_mid=ny/nsparse_phi_mid	;Number of arrows in phi
  N_phi_out=ny/nsparse_phi_out	;Number of arrows in phi
  nr_in=fix(N_r*n_in)		;Sets radius below which nsparse_phi_in is used
  nr_out=fix(N_r*n_mid)	;Sets radius below which nsparse_phi_mid is used and above which nsparse_phi_out is used
  frac=0.1;0.06	;The larger this number, the longer the arrows
  ;
  r_plot=0.45		;Sets the fraction of the disc radius up to which to plot (in cartesian coordinates plot from -r_plot to r_plot for x, y)
  xcoord_prime=dblarr(nx,ny)
  ycoord_prime=dblarr(nx,ny)
  xcoord_sparse_in=dblarr(nx,ny)
  ycoord_sparse_in=dblarr(nx,ny)
  xcoord_sparse_mid=dblarr(nx,ny)
  ycoord_sparse_mid=dblarr(nx,ny)
  xcoord_sparse_out=dblarr(nx,ny)
  ycoord_sparse_out=dblarr(nx,ny)
  xcoord_prime_sparse_in=dblarr(nx,ny)
  ycoord_prime_sparse_in=dblarr(nx,ny)
  xcoord_prime_sparse_mid=dblarr(nx,ny)
  ycoord_prime_sparse_mid=dblarr(nx,ny)
  xcoord_prime_sparse_out=dblarr(nx,ny)
  ycoord_prime_sparse_out=dblarr(nx,ny)
  Bx=dblarr(nx,ny)
  By=dblarr(nx,ny)
  B_max=max(sqrt(Brnew^2+Bpnew^2))
  Bx= Brnewcorotate*cos(phinewcorotate[*,0:ny-1]) -Bpnewcorotate*sin(phinewcorotate[*,0:ny-1])
  By= Brnewcorotate*sin(phinewcorotate[*,0:ny-1]) +Bpnewcorotate*cos(phinewcorotate[*,0:ny-1])
  xcoord_prime= xcoord +frac*Bx/B_max
  ycoord_prime= ycoord +frac*By/B_max
  for k=0,nx/nsparse_r-2 do begin
  for j=0,ny/nsparse_phi_in-1 do begin
    xcoord_sparse_in(k,j)=xcoord(nsparse_r*(k+1),nsparse_phi_in*j)
    ycoord_sparse_in(k,j)=ycoord(nsparse_r*(k+1),nsparse_phi_in*j)
    if (abs(xcoord_sparse_in(k,j)) lt r_disk*r_plot and abs(ycoord_sparse_in(k,j)) lt r_disk*r_plot) then begin	;to avoid plotting arrows outside of the plotting region
      xcoord_prime_sparse_in(k,j)=xcoord_prime(nsparse_r*(k+1),nsparse_phi_in*j)
      ycoord_prime_sparse_in(k,j)=ycoord_prime(nsparse_r*(k+1),nsparse_phi_in*j)
    endif else begin
      xcoord_sparse_in(k,j)=0.		;to avoid plotting any arrow at all
      ycoord_sparse_in(k,j)=0.		;to avoid plotting any arrow at all
      xcoord_prime_sparse_in(k,j)=xcoord_sparse_in(k,j)		;to avoid plotting any arrow at all
      ycoord_prime_sparse_in(k,j)=xcoord_sparse_in(k,j)		;to avoid plotting any arrow at all
    endelse
  endfor
  for j=0,ny/nsparse_phi_mid-1 do begin
    xcoord_sparse_mid(k,j)=xcoord(nsparse_r*(k+1),nsparse_phi_mid*j)
    ycoord_sparse_mid(k,j)=ycoord(nsparse_r*(k+1),nsparse_phi_mid*j)
    if (abs(xcoord_sparse_mid(k,j)) lt r_disk*r_plot and abs(ycoord_sparse_mid(k,j)) lt r_disk*r_plot) then begin
      xcoord_prime_sparse_mid(k,j)=xcoord_prime(nsparse_r*(k+1),nsparse_phi_mid*j)
      ycoord_prime_sparse_mid(k,j)=ycoord_prime(nsparse_r*(k+1),nsparse_phi_mid*j)
    endif else begin
      xcoord_sparse_mid(k,j)=1.e6	;arbitrary large number (don't set to zero because then would plot an arrow at the origin)
      ycoord_sparse_mid(k,j)=1.e6 ;arbitrary large number (don't set to zero because then would plot an arrow at the origin)
      xcoord_prime_sparse_mid(k,j)=xcoord_sparse_mid(k,j)
      ycoord_prime_sparse_mid(k,j)=ycoord_sparse_mid(k,j)
    endelse
  endfor
  for j=0,ny/nsparse_phi_out-1 do begin
    xcoord_sparse_out(k,j)=xcoord(nsparse_r*(k+1),nsparse_phi_out*j)
    ycoord_sparse_out(k,j)=ycoord(nsparse_r*(k+1),nsparse_phi_out*j)
    if (abs(xcoord_sparse_out(k,j)) lt r_disk*r_plot and abs(ycoord_sparse_out(k,j)) lt r_disk*r_plot) then begin
      xcoord_prime_sparse_out(k,j)=xcoord_prime(nsparse_r*(k+1),nsparse_phi_out*j)
      ycoord_prime_sparse_out(k,j)=ycoord_prime(nsparse_r*(k+1),nsparse_phi_out*j)
    endif else begin
      xcoord_sparse_out(k,j)=1.e6	;arbitrary large number (don't set to zero because then would plot an arrow at the origin)
      ycoord_sparse_out(k,j)=1.e6	;arbitrary large number (don't set to zero because then would plot an arrow at the origin)
      xcoord_prime_sparse_out(k,j)=xcoord_sparse_out(k,j)
      ycoord_prime_sparse_out(k,j)=ycoord_sparse_out(k,j)
    endelse
  endfor
  endfor
  ;
  ; IV) Contour B/B0 (x-y)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_Bcartcorotate_',ntstr,'.eps']),xsize=xsize,ysize=ysize,bits=8,/color
      char=float(nvert)*charfac
      !p.charsize=char
      hsize_default=!d.x_size/64
      hsize=-0.75;0.75*hsize_default;-0.5;-0.25	;negative number makes head size scale with arrow length
      ;hsize=0.25		;negative number makes head size scale with arrow length
      hthick=3.5;2.0;1.0;0.5	;changed July 26, 2012 from 1.0 to 2.0
      thick=-3.5;-2.0;-1.0;2.0;1.0	;changed July 26, 2012 from -1.0 to -2.0
      arrowcol=col0;col1;-1
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
      hsize=2.
      hthick=1.
      thick=0.5
      arrowcol=col0
    endelse
  ;
    print,'max(B)=',max(Bnewcorotate[nxghost:nx-1-nxghost,*]/Beq0)
    print,'max(delta)=',max(deltanewcorotate[nxghost:nx-1-nxghost,*])
    print,'min(delta)=',min(deltanewcorotate[nxghost:nx-1-nxghost,*])
    print,'rmax=',rnew[deltawherexcorotate_half1,0]*r_disk_kpc
    print,'phimax=',180./!pi*phinewcorotate[0,deltawhereycorotate_half1]
    print,''
    nxphys=nx-2*nxghost
    print,'max(Bzmod(r>3kpc))=',max(Bznewcorotate[nxghost+nxphys/5:nx-1-nxghost,*]/Beq0)
    print,'max(Br(r<3kpc))='   ,max(Brnewcorotate[nxghost+nxphys/5:nx-1-nxghost,*]/Beq0)

    loadct,Bcol,/silent
    contour,Bnewcorotate[nxghost:nx-1-nxghost,*]/Beq0,xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=VB,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,/fill,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    if (Plot_part eq 1) then begin
      if (n_arm[0] eq 2) then begin
        contour,alp_k_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
      endif else begin
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endelse
      if (nspiral gt 1) then begin
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,1],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endif
    endif
    loadct,5,/silent
    ;contour,Uznew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,Uznew_normcorotate[ixU:nx-1-nxghost,*],xcoord[ixU:nx-1-nxghost,*]*r_disk_kpc,ycoord[ixU:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm[ixU:nx-1-nxghost,*]),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,alp_knew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(alp_knew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1
    if (Arrow eq 1) then begin
      arrow,xcoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,xcoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,/data,thick=thick,hsize=hsize,hthick=hthick,col=arrowcol
    endif
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      endif
    endif else begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
      endif
    endelse
    if (Corotate ne 1) then begin
    rclong_kpc =dblarr(nspiral,1000001)
      for isp=0,nspiral-1 do begin
        ;rclong_kpc[isp,*]=rc_kpc[isp]
        rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
      endfor
      xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
      for isp=0,nspiral-1 do begin
        oplot,xcoordlong_kpc,+sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
        oplot,xcoordlong_kpc,-sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
      endfor
    endif
    ;
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                    ,col=col0,charth=cth
    ;xyouts,0.40,0.85,/normal,'Bmax/B0='                               ,col=col1,charth=cth
    ;xyouts,0.55,0.85,/normal,Bmax                                     ,col=col1,charth=cth
    ;xyouts,0.47,0.78,/normal,s0                                       ,col=col1,charth=cth
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col1
    ;
  endfor
  ;
  ; IV) Contour pB (normalized, x-y)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_pBcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    ;
    contour,pBnewcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=VpB,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,/fill,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    if (Plot_part eq 1) then begin
      if (n_arm[0] eq 2) then begin
        loadct,38,/silent
        contour,alp_k_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        loadct,3,/silent
      endif else begin
        contour,alp_k_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endelse
      if (nspiral gt 1) then begin
        contour,alp_k_partnew_normcorotate[*,*,1],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[*,*,1],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endif
    endif
    loadct,5,/silent
    if (Corotate ne 1) then begin
      rclong_kpc =dblarr(nspiral,1000001)
      for isp=0,nspiral-1 do begin
        ;rclong_kpc[isp,*]=rc_kpc[isp]
        rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
      endfor
      xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
      for isp=0,nspiral-1 do begin
        oplot,xcoordlong_kpc,+sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
        oplot,xcoordlong_kpc,-sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
      endfor
    endif
    contour,alp_knew_normcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alpnormmax*max(alp_knew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1
    ;contour,Uznew_normcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=Uzmax*max(Uznew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
      endif
    endelse
    ;
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                    ,col=col0,charth=cth
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col0
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col0
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col0
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col0
    ;
  endfor
  ;
  ; V) Contour pB (phi-log10(r))
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_pBlogcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
  ;
    loadct,5,/silent
    contour,transpose(pBnewcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=100,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',/cell_fill,levels=VpB,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    ;
    if (ts_t[nt] gt ti_spiral[0] and ts_t[nt] lt tf_spiral[0]) then begin
      if (Corotate ne 1 and Windup ne 1) then begin
        for isp=0,nspiral-1 do begin
          oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(r_disk_kpc*ts_rc[nt,isp]),li=1,col=col1,th=4./3*cth
        endfor
      endif
      contour,transpose(alp_knew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpnormmax*max(alp_knew_norm),/overplot,c_th=2*cth,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      if (n_arm[0] eq 2) then begin
        loadct,38,/silent
        contour,transpose(alp_k_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=5,col=col2,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      endif else begin
        loadct,3,/silent
        contour,transpose(alp_k_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=3,col=col5,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      endelse
      if (nspiral gt 1) then begin
        loadct,3,/silent
        contour,transpose(alp_k_partnew_normcorotate[*,*,1]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=3,col=col5,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      endif
    endif
    !p.charsize=xchar
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    ;if (deltawherex[0] ne -1 and deltawherey[0] ne -1 and sizewherex[1] lt 4) then begin
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
      endif
    endelse
    ;endif
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
  endfor
  ;
  ; VI) Contour alpha (normalized, x-y)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_alphacorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    ;
    loadct,5,/silent
    ;contour,alp_knew_normcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=VB/Bmax,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,/fill,c_th=6,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=3
    contour,alp_knew_normcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=0.50,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,c_th=2*cth,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    ;contour,alp_knew_normcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alpnormmax*max(alp_knew_norm),/overplot,col=col0,c_li=0,c_th=6
    if (n_arm[0] eq 2) then begin
      loadct,38,/silent
      contour,alp_k_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,col=col2,c_li=5,c_th=cth,xstyle=1,ystyle=1;,/closed;,/cell_fill	;overplot orig alpha contour (max of matter spiral)
      loadct,39,/silent
    endif else begin
      contour,alp_k_partnew_normcorotate[*,*,0],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,col=col3,c_li=3,c_th=cth,xstyle=1,ystyle=1;,/closed;,/cell_fill	;overplot orig alpha contour (max of matter spiral)
    endelse
    if (nspiral gt 1) then begin
      contour,alp_k_partnew_normcorotate[*,*,1],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,col=col3,c_li=3,c_th=cth,xstyle=1,ystyle=1;,/closed;,/cell_fill	;overplot orig alpha contour (max of matter spiral)
      if (nspiral gt 2) then begin
        contour,alp_k_partnew_normcorotate[*,*,2],xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,col=col3,c_li=4,c_th=2./3*cth,xstyle=1,ystyle=1;,/closed;,/cell_fill	;overplot orig alpha contour (max of matter spiral)
      endif
    endif
    if (Corotate ne 1) then begin
      rclong_kpc =dblarr(nspiral,1000001)
      for isp=0,nspiral-1 do begin
        ;rclong_kpc[isp,*]=rc_kpc[isp]
        rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
      endfor
      xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
      for isp=0,nspiral-1 do begin
        oplot,xcoordlong_kpc,+sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col0
        oplot,xcoordlong_kpc,-sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col0
      endfor
    endif
    ;contour,deltanewcorotate,xcoord*r_disk_kpc,ycoord*r_disk_kpc,levels=0.1,/overplot,col=col0,c_th=4;,/cell_fill	;overplot orig alpha contour (max of matter spiral)
    ;arrow,xcoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,xcoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,/data,thick=thick,hsize=hsize,hthick=hthick,col=arrowcol
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
      endif
    endelse
    ;
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=-1
  endfor
    ;
  ; VII) Contour normalized alpha (phi-log10(r))
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_alphanormlogcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
  ;
    loadct,5,/silent
    ;contour,transpose(alp_knew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=100,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',/cell_fill,levels=VBlog/Bmaxlog,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],charth=3
    contour,transpose(alp_knew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=100,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',levels=0.50,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    ;contour,transpose(alp_knew_normcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpnormmax*max(alp_knew_norm),/overplot,c_th=6,c_li=0,col=col0	;overplot orig alpha contour (max of matter spiral)
    ;
    if (ts_t[nt] gt ti_spiral[0] and ts_t[nt] lt tf_spiral[0]) then begin
      if (Corotate ne 1 and Windup ne 1) then begin
        for isp=0,nspiral-1 do begin
          ;oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(rc_kpc[isp]),li=1,col=col0,th=2
          oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(r_disk_kpc*ts_rc[nt,isp]),li=1,col=col0,th=4./3*cth
        endfor
      endif
      ;contour,transpose(deltanewcorotate),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=0.1,/overplot,c_th=4,col=col1	;overplot orig alpha contour (max of matter spiral)
      if (n_arm[0] eq 2) then begin
        loadct,38,/silent
        contour,transpose(alp_k_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=5,col=col2,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        loadct,39,/silent
      endif else begin
        loadct,39,/silent
        contour,transpose(alp_k_partnew_normcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=3,col=col3,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      endelse
      if (nspiral gt 1) then begin
        contour,transpose(alp_k_partnew_normcorotate[*,*,1]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=cth,c_li=3,col=col3,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        if (nspiral gt 2) then begin
          contour,transpose(alp_k_partnew_normcorotate[*,*,2]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alppartnormmax*max(alp_knew_norm),/overplot,c_th=2./3*cth,c_li=4,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        endif
      endif
    endif
    !p.charsize=xchar
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    ;if (deltawherex[0] ne -1 and deltawhereycorotate[0] ne -1 and sizewherex[1] lt 4) then begin
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
      endif
    endelse
    ;endif
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
  endfor
  ;
  ; VIII) Contour alpha (phi-log10(r))
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_alphalogcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    ;
    loadct,5,/silent
    if (ts_t[nt] gt ti_spiral[0] and ts_t[nt] lt tf_spiral[0]) then begin
      if (Corotate ne 1 and Windup ne 1) then begin
        contour,transpose(alp_knew),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=100,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',/cell_fill,levels=alpmax*VBlog/Bmaxlog,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
        for isp=0,nspiral-1 do begin
          ;oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(rc_kpc[isp]),li=1,col=col0,th=2
          oplot,phinewcorotate[0,*]*180./!pi,phinewcorotate[0,*]*0+alog10(r_disk_kpc*ts_rc[nt,isp]),li=1,col=col1,th=4./3*cth
        endfor
        contour,transpose(alp_knew),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),nlevels=100,xr=yr,yr=xlogr,xtit='!7u!6 [degrees]',ytit='!6Log!D10!N!8r!6 [kpc]',/cell_fill,xticks=4,xminor=3,yticks=7,yminor=5,pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
      endif else begin
      endelse
      contour,transpose(Bnewcorotate)/Beq0,transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=0.3*max(Bnew/Beq0),/overplot,c_th=4./3*cth,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      if (n_arm[0] eq 2) then begin
        loadct,38,/silent
        contour,transpose(alp_k_partnewcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpmax*max(alp_knew),/overplot,c_th=cth,c_li=5,col=col2,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        loadct,5,/silent
      endif else begin
        contour,transpose(alp_k_partnewcorotate[*,*,0]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpmax*max(alp_knew),/overplot,c_th=cth,c_li=3,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
      endelse
      if (nspiral gt 1) then begin
        contour,transpose(alp_k_partnewcorotate[*,*,1]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpmax*max(alp_knew),/overplot,c_th=cth,c_li=3,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        if (nspiral gt 2) then begin
          contour,transpose(alp_k_partnewcorotate[*,*,2]),transpose(phinewcorotate*180./!pi),transpose(alog10(rnew*r_disk_kpc)),levels=alpmax*max(alp_knew),/overplot,c_th=2./3*cth,c_li=4,col=col0,xstyle=1,ystyle=1	;overplot orig alpha contour (max of matter spiral)
        endif
      endif
    endif
    !p.charsize=xchar
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                  ,col=col0,charth=cth
    ;if (deltawherex[0] ne -1 and deltawhereycorotate[0] ne -1 and sizewherex[1] lt 4) then begin
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      oplot,180./!pi*phinew[0,deltawhereycorotate],alog10(rnew[deltawherex,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
    endif else begin
      if (diamond eq 1) then begin
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half1],alog10(rnew[deltawherexcorotate_half1,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherexcorotate_half1,deltawhereycorotate_half1]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
        oplot,180./!pi*phinewcorotate[0,deltawhereycorotate_half2],alog10(rnew[deltawherexcorotate_half2,0]*r_disk_kpc),psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherexcorotate_half2,deltawhereycorotate_half2]
      endif
    endelse
    ;endif
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,xth=5./3*cth,color=-1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=4,xminor=3,yticks=7,yminor=5,yth=5./3*cth,color=-1
  endfor
;
  ; IX) Contour |Bz|/B0 (x-y)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_Bzcartcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
      hsize_default=!d.x_size/64
      hsize=-0.75;0.75*hsize_default;-0.5;-0.25	;negative number makes head size scale with arrow length
      ;hsize=0.25		;negative number makes head size scale with arrow length
      hthick=3.5;2.0;1.0;0.5	;changed July 26, 2012 from 1.0 to 2.0
      thick=-3.5;-2.0;-1.0;2.0;1.0	;changed July 26, 2012 from -1.0 to -2.0
      arrowcol=col0;col1;-1
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
      hsize=2.
      hthick=1.
      thick=0.5
      arrowcol=col0
    endelse
    for ix=nxghost,nx-1-nxghost do begin
      ;print,'r, Bz/Br, Bz/B',r[ix,0]*r_disk_kpc,Bznewcorotate[ix,0]/Brnewcorotate[ix,0],Bznewcorotate[ix,0]/Bnewcorotate[ix,0]
    endfor
  ;
    loadct,Bcol,/silent
    contour,Bznewcorotate[nxghost:nx-1-nxghost,*]/Beq0,xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=VB,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,/fill,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    if (Plot_part eq 1) then begin
      if (n_arm[0] eq 2) then begin
        contour,alp_k_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
      endif else begin
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endelse
      if (nspiral gt 1) then begin
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,1],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endif
    endif
    loadct,5,/silent
    ;contour,Uznew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,Uznew_normcorotate[ixU:nx-1-nxghost,*],xcoord[ixU:nx-1-nxghost,*]*r_disk_kpc,ycoord[ixU:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm[ixU:nx-1-nxghost,*]),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,alp_knew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(alp_knew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1
    if (Arrow eq 1) then begin
      arrow,xcoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,xcoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,/data,thick=thick,hsize=hsize,hthick=hthick,col=arrowcol
    endif
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      endif
    endif else begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
      endif
    endelse
    if (Corotate ne 1) then begin
    rclong_kpc =dblarr(nspiral,1000001)
      for isp=0,nspiral-1 do begin
        ;rclong_kpc[isp,*]=rc_kpc[isp]
        rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
      endfor
      xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
      for isp=0,nspiral-1 do begin
        oplot,xcoordlong_kpc,+sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
        oplot,xcoordlong_kpc,-sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
      endfor
    endif
    ;
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                    ,col=col0,charth=cth
    ;xyouts,0.40,0.85,/normal,'Bmax/B0='                               ,col=col1,charth=cth
    ;xyouts,0.55,0.85,/normal,Bmax                                     ,col=col1,charth=cth
    ;xyouts,0.47,0.78,/normal,s0                                       ,col=col1,charth=cth
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col1
    ;
  endfor
  ;
  ; X) Contour |Br|/B0 (x-y)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=xsizeps & ysize=float(nvert)*xsizeps*asp
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_Brcartcorotate_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      char=float(nvert)*charfac
      !p.charsize=char
      hsize_default=!d.x_size/64
      hsize=-0.75;0.75*hsize_default;-0.5;-0.25	;negative number makes head size scale with arrow length
      ;hsize=0.25		;negative number makes head size scale with arrow length
      hthick=3.5;2.0;1.0;0.5	;changed July 26, 2012 from 1.0 to 2.0
      thick=-3.5;-2.0;-1.0;2.0;1.0	;changed July 26, 2012 from -1.0 to -2.0
      arrowcol=col0;col1;-1
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=xsize_screen & ysize=float(nvert)*xsize*asp
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
      hsize=2.
      hthick=1.
      thick=0.5
      arrowcol=col0
    endelse
  ;
    loadct,Bcol,/silent
    contour,abs(Brnewcorotate[nxghost:nx-1-nxghost,*])/Beq0,xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=VB,xr=box*r_disk_kpc,xticks=6,xminor=3,yr=box*r_disk_kpc,yticks=6,yminor=3,/fill,xtit='!8x !6(kpc)',ytit='!8y !6(kpc)',pos=[0.20*asp,0.15,0.95*asp,0.90],charth=cth,xstyle=1,ystyle=1
    if (Plot_part eq 1) then begin
      if (n_arm[0] eq 2) then begin
        contour,alp_k_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=alppartnormmax*max(alp_knew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=5,c_th=cth,col=col2,xstyle=1,ystyle=1
      endif else begin
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,0],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endelse
      if (nspiral gt 1) then begin
        contour,Uz_partnew_normcorotate[nxghost:nx-1-nxghost,*,1],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzpartnormmax*max(Uznew_norm),/overplot,c_li=3,c_th=cth,col=col5,xstyle=1,ystyle=1
      endif
    endif
    loadct,5,/silent
    ;contour,Uznew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,Uznew_normcorotate[ixU:nx-1-nxghost,*],xcoord[ixU:nx-1-nxghost,*]*r_disk_kpc,ycoord[ixU:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(Uznew_norm[ixU:nx-1-nxghost,*]),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1,/fill
    contour,alp_knew_normcorotate[nxghost:nx-1-nxghost,*],xcoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,ycoord[nxghost:nx-1-nxghost,*]*r_disk_kpc,levels=Uzmax*max(alp_knew_norm),/overplot,col=col0,c_li=0,c_th=2*cth,xstyle=1,ystyle=1
    if (Arrow eq 1) then begin
      arrow,xcoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,xcoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,ycoord_prime_sparse_mid[nr_in:nr_out-1,0:N_phi_mid-1]*r_disk_kpc,/data,thick=thick,hsize=hsize,hthick=hthick,col=arrowcol
    endif
    ;
    if (nspiral eq 1 and ts_t[nt] gt ti_spiral[0]) then begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
        oplot,xcoord[deltawherecorotate]*r_disk_kpc,ycoord[deltawherecorotate]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate[0]],symsize=symsize*deltanewcorotate[deltawherecorotate[0]]
      endif
    endif else begin
      if (diamond eq 1) then begin
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half1]*r_disk_kpc,ycoord[deltawherecorotate_half1]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half1],symsize=symsize*deltanewcorotate[deltawherecorotate_half1]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym1,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
        oplot,xcoord[deltawherecorotate_half2]*r_disk_kpc,ycoord[deltawherecorotate_half2]*r_disk_kpc,psym=sym2,col=symcol,th=symth*deltanewcorotate[deltawherecorotate_half2],symsize=symsize*deltanewcorotate[deltawherecorotate_half2]
      endif
    endelse
    if (Corotate ne 1) then begin
    rclong_kpc =dblarr(nspiral,1000001)
      for isp=0,nspiral-1 do begin
        ;rclong_kpc[isp,*]=rc_kpc[isp]
        rclong_kpc[isp,*]=r_disk_kpc*ts_rc[nt,isp]
      endfor
      xcoordlong_kpc=findgen(1000001)*2*r_disk_kpc/1000001-r_disk_kpc
      for isp=0,nspiral-1 do begin
        oplot,xcoordlong_kpc,+sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
        oplot,xcoordlong_kpc,-sqrt(rclong_kpc[isp,*]^2-xcoordlong_kpc^2),th=4./3*cth,li=1,col=col1
      endfor
    endif
    ;
    xyouts,0.44,0.92,/normal,strmid(strtrim(string(nt*0.025),2),0,ndig),col=col0,charth=cth
    xyouts,0.62,0.92,/normal,'Gyr'                                    ,col=col0,charth=cth
    ;xyouts,0.40,0.85,/normal,'Bmax/B0='                               ,col=col1,charth=cth
    ;xyouts,0.55,0.85,/normal,Bmax                                     ,col=col1,charth=cth
    ;xyouts,0.47,0.78,/normal,s0                                       ,col=col1,charth=cth
    !p.charsize=char
    axis,xaxis=0,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col1
    axis,xaxis=1,xtickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,xth=5./3*cth,color=col1
    axis,yaxis=0,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col1
    axis,yaxis=1,ytickname=replicate(' ',30),xticks=6,xminor=3,yticks=6,yminor=3,yth=5./3*cth,color=col1
    ;
  endfor
  ;
  ;
  divs1=5
  divs2=2
  divs3=4
  divspB=5
  colr_ratio=5.5
  colr_charsize=0.6
  ;
  ;  COLORBARI) COLOUR PLOT OF SQUARE ROOT OF MAGNETIC ENERGY (cartesian)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=6. & ysize=xsize/colr_ratio
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_colr_Bcart_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      charfac=1.5
      char=float(nvert)*charfac
      xchar=1.5
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=250. & ysize=xsize/5
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    ;
    loadct,Bcol,/silent
    colorbar, color=0, minrange=0., maxrange=Bmax, pos=[0.20*asp,0.35,0.95*asp,0.50], divisions=divs1,format='(f5.2)',/pscolor,charsize=colr_charsize,ncolors=256,tit='!8B/B!D!60!N',thick=2,charth=1.5
  endfor
  
  ;  COLORBARII) COLOUR PLOT OF SQUARE ROOT OF MAGNETIC ENERGY (log r vs phi plot)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=6. & ysize=xsize/colr_ratio
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_colr_B_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      charfac=1.5
      char=float(nvert)*charfac
      xchar=1.5
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=250. & ysize=xsize/5
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    loadct,5,/silent
    colorbar, color=0, minrange=0., maxrange=Bmaxlog, pos=[0.20*asp,0.60,0.95*asp,0.75], divisions=divs2,format='(f4.2)',/pscolor,charsize=1.5,ncolors=256,tit='!8B/B!D!60!N',thick=2,charth=3
  endfor
  
  ;  COLORBARIII) COLOUR PLOT OF delta (log r vs phi plot)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=6. & ysize=xsize/colr_ratio
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_colr_delta_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      charfac=1.5
      char=float(nvert)*charfac
      xchar=1.5
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=250. & ysize=xsize/5
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    loadct,5,/silent
    colorbar, color=0, minrange=-deltamax, maxrange=deltamax, pos=[0.20*asp,0.35,0.95*asp,0.50], divisions=divs3,format='(f4.1)',/pscolor,charsize=colr_charsize,ncolors=256,tit='!7d!6',thick=2,charth=1.5
  endfor
  ;
  ;  COLORBARIV) COLOUR PLOT OF pB (cartesian)
  for iplot=0,1 do begin
    if (iplot eq 0) then begin
      set_plot,'ps'
      cd,plotdirec
      xsize=12. & ysize=xsize/colr_ratio
      ntstr=strjoin(['nt',strtrim(string(nt),2)])
      device,filename=strjoin([s0,s2,'_colr_pitch_',ntstr,'.eps']),/color,xsize=xsize,ysize=ysize,bits=8
      charfac=1.5
      char=float(nvert)*charfac
      xchar=1.5
      !p.charsize=char
    endif else begin
      device,/close
      set_plot,'x'
      !p.charsize=1.0
      xsize=250. & ysize=xsize/5
      window,8,xsize=xsize,ysize=ysize,xpos=xpos_screen,ypos=ypos_screen
    endelse
    loadct,5,/silent
    colorbar, color=0, minrange=min(VpB), maxrange=max(VpB), pos=[0.20*asp,0.60,0.95*asp,0.75], divisions=divspB,format='(I3)',/pscolor,charsize=1.5,ncolors=256,tit='!8-p!DB!N!6 [degrees]',thick=5,charth=3
  endfor
endfor
cd,programdirec
end
