s0='Uz_155'
s1='/Users/luke' ;'/home/test' ;'/data1/student/cluke'
s2=''
datadirec=strjoin([s1,'/fortran_pde/2D/telegraph/r_phi_noz/data/',s0])
programdirec=strjoin([s1,'/fortran_pde/2D/telegraph/r_phi_noz/program/',s0])
cd,datadirec
openr,10,strjoin(['param',s2,'.out'])
openr,11,strjoin(['init' ,s2,'.out'])
openr,12,strjoin(['run'  ,s2,'.out'])
openr,13,strjoin(['ts'   ,s2,'.out'])
xchar=1.4
t1=8.	;time 1 for calculating growth rate Gamma
t2=9.5	;time 2 for calculating growth rate Gamma
loadct,5
nspiral=1
n_arm=        intarr(nspiral)
rc0_kpc=      fltarr(nspiral)
kspiral=      fltarr(nspiral)
wamp=         fltarr(nspiral)
twon_fac=     fltarr(nspiral)
eps_alp=      fltarr(nspiral)
eps_Uz=       fltarr(nspiral)
eps_Beq=      fltarr(nspiral)
ti_spiral=    fltarr(nspiral)
tf_spiral=    fltarr(nspiral)
Omega_kmskpc= fltarr(nspiral)
Omega=        fltarr(nspiral)
rc_kpc=       fltarr(nspiral)
rc=           fltarr(nspiral)
readf,10,nx, ny, h0_kpc, h0_km, td_s, r_disk_kpc, Beq0, r_quench, tau, dt, tsnap, nvar, n1, n2, U0_kms, l0_kpc, n_arm, rc0_kpc, kspiral, wamp, twon_fac, eps_alp, eps_Uz, eps_Beq, ti_spiral, tf_spiral, mmax, C_alp, td_Gyr, C_etat, v_turb_kms, v_turb, Omega_kmskpc, Omega, r_om_kpc, r_om, om0_kmskpc, om0, rc_kpc, rc, Uz_thresh,nxghost
print,'td_Gyr=',td_Gyr
nx=fix(nx) & ny=fix(ny) & nvar=fix(nvar) & n1=fix(n1) & n2=fix(n2) & n_arm=fix(n_arm) & mmax=fix(mmax)
r=            fltarr(nx,ny)
phi=          fltarr(nx,ny)
h=            fltarr(nx,ny)
om=           fltarr(nx,ny)
Uz=           fltarr(nx,ny)
alp_k=        fltarr(nx,ny)
alp_m=        fltarr(nx,ny)
Beq=          fltarr(nx,ny)
Br_init=      fltarr(nx,ny)
Bp_init=      fltarr(nx,ny)
Br=           fltarr(nx,ny)
Bp=           fltarr(nx,ny)
Bzero=        fltarr(3,nx)
Bmmax=        fltarr(mmax,3,nx)
ts_t=         fltarr(n1+1) 
ts_rc=        fltarr(n1+1,nspiral) 
ts_Omega=     fltarr(n1+1,nspiral) 
ts_Br=        fltarr(n1+1,nx,ny)
ts_Bp=        fltarr(n1+1,nx,ny)
delta=        fltarr(n1+1,nx,ny)
ts_alp_m=     fltarr(n1+1,nx,ny)
ts_alp_k=     fltarr(n1+1,nx,ny)
ts_Uz=        fltarr(n1+1,nx,ny)
ts_alp_k_part=fltarr(n1+1,nx,ny,nspiral)
ts_Uz_part=   fltarr(n1+1,nx,ny,nspiral)
ts_rpeak=     fltarr(n1+1,nspiral) 
ts_Bzmod=     fltarr(n1+1,nx,ny) 
ts_alp=       fltarr(n1+1,nx,ny)
ts_Bzero=     fltarr(n1+1,3,nx)
ts_Bmmax=     fltarr(n1+1,mmax,3,nx)
readf,11,r,phi,h,om,Uz,alp_k,Beq,Br_init,Bp_init
if (nvar eq 2 or nvar eq 4) then begin	;kinematic regime
  readf,12,t,Br,Bp,Bzero,Bmmax
  readf,13,ts_t,ts_rc, ts_Omega, ts_Br,ts_Bp,ts_Bzero,ts_Bmmax,ts_Uz,ts_Uz_part,ts_rpeak,ts_Bzmod
  alp_m=0.*alp_k
endif else begin
  readf,12,t,Br,Bp,alp_m,Bzero,Bmmax
  readf,13,ts_t,ts_rc,ts_Omega,ts_Br,ts_Bp,ts_alp_m,ts_alp_k,ts_alp_k_part,ts_Bzero,ts_Bmmax,ts_Uz,ts_Uz_part,ts_rpeak,ts_Bzmod
endelse
alp=alp_k +alp_m
ts_alp=ts_alp_k +ts_alp_m
for ifile=10,13 do begin
  close,ifile
endfor
cd,programdirec
end
