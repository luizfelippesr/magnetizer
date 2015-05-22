s0='CE_010'
s1='/Users/luke/fortran_pde/1D/telegraph/r_noz/'
s2='_00000001'
simulationoutputdirec=   strjoin([s1,'output/'   ,s0])
programdirec=strjoin([s1,'program/',s0])
cd,simulationoutputdirec
openr,10,strjoin(['param',s2,'.out'])
openr,11,strjoin(['init' ,s2,'.out'])
openr,12,strjoin(['run'  ,s2,'.out'])
openr,13,strjoin(['ts'   ,s2,'.out'])
xchar=1.4
t1=8.	;time 1 for calculating growth rate Gamma
t2=9.5	;time 2 for calculating growth rate Gamma
loadct,5
;read parameter values
readf,10,t0_Gyr,t0_kpcskm,h0_kpc,etat0_cm2s,n0_cm3,B0_mkG
readf,10,nvar,dt,n1,n2,dx,nxphys,nxghost,nx
readf,10,l_sol_kpc,r_l_kpc,v_sol_kms,r_v_kpc,n_sol_cm3,r_n_kpc,Uz_sol_kms,r_Uz_kpc,h_sol_kpc,r_h_kpc,Uphi_sol_kms,r_om_kpc,R_kappa
readf,10,r_in,r_disk_kpc,r_sol_kpc,r1_kpc,ctau,nn,lambda,C_alp,alpceil,Rm_inv
readf,10,etat_sol_cm2s,td_sol_Gyr,tau_sol_Gyr,etat_sol,td_sol,tau_sol
close,10
print,'UNITS'
print,'t0 (Gyr)=',t0_Gyr,'   t0 (kpc s/km)=',t0_kpcskm,'   h0 (kpc)=',h0_kpc,'   etat0 (cm2s)= ',etat0_cm2s
print,''
print,'NUMERICS:'
print,'nvar=     ',nvar     
print,'dt=       ',dt       ,'   n1=       ',n1       ,'   n2=     ',n2
print,'dx=       ',dx       ,'   nxphys=   ',nxphys   ,'   nxghost=',nxghost,'   nx=',nx
print,''
print,'INPUT PARAMETERS:'
print,'r_disk_kpc= ',r_disk_kpc ,'   R_kappa=',R_kappa
print,'l_sol_kpc=  ',l_sol_kpc  ,'   r_l_kpc=',r_l_kpc,'   v_sol_kms=   ',v_sol_kms    ,'   r_v_kpc= ',r_v_kpc
print,'n_sol_cm3=  ',n_sol_cm3  ,'   r_n_kpc=',r_n_kpc,'   Uz_sol_kms=  ',Uz_sol_kms   ,'   r_Uz_kpc=',r_Uz_kpc
print,'h_sol_kpc=  ',h_sol_kpc  ,'   r_h_kpc=',r_h_kpc,'   Uphi_sol_kms=',Uphi_sol_kms ,'   r_om_kpc=',r_om_kpc
print,''
print,'OTHER FREE PARAMETERS:'
print,'r_in= ',r_in ,'   r_disk_kpc=',r_disk_kpc,'   r_sol_kpc=',r_sol_kpc,'   r1_kpc=',r1_kpc
print,'ctau= ',ctau ,'   nn=        ',nn        ,'   lambda=   ',lambda
print,'C_alp=',C_alp,'   alpceil=   ',alpceil   ,'   Rm_inv=   ',Rm_inv
print,''
print,'IMPORTANT CALCULATED PARAMETERS'
print,'etat_sol_cm2s=   ',etat_sol_cm2s,'   td_sol_Gyr=',td_sol_Gyr,'   tau_sol_Gyr=',tau_sol_Gyr
print,'etat_sol=        ',etat_sol     ,'   td_sol=    ',td_sol    ,'   tau_sol     ',tau_sol
nx=fix(nx) & nvar=fix(nvar) & n1=fix(n1) & n2=fix(n2)
;
r=           fltarr(nx)
h_init=      fltarr(nx)
om_init=     fltarr(nx)
G_init=      fltarr(nx)
l_init=      fltarr(nx)
v_init=      fltarr(nx)
etat_init=   fltarr(nx)
tau_init=    fltarr(nx)
alp_k_init=  fltarr(nx)
Uz_init=     fltarr(nx)
Ur_init=     fltarr(nx)
alp_m_init=  fltarr(nx)
n_init=      fltarr(nx)
Beq_init=    fltarr(nx)
Br_init=     fltarr(nx)
Bp_init=     fltarr(nx)
Bzmod_init=  fltarr(nx)
;
Br=          fltarr(nx)
Bp=          fltarr(nx)
h=           fltarr(nx)
om=          fltarr(nx)
G=           fltarr(nx)
l=           fltarr(nx)
v=           fltarr(nx)
etat=        fltarr(nx)
tau=         fltarr(nx)
alp_k=       fltarr(nx)
Uz=          fltarr(nx)
Ur=          fltarr(nx)
alp_m=       fltarr(nx)
n=           fltarr(nx)
Beq=         fltarr(nx)
Br=          fltarr(nx)
Bp=          fltarr(nx)
Bzmod=       fltarr(nx)
alp=         fltarr(nx)
;
ts_t=        fltarr(n1+1) 
ts_Br=       fltarr(n1+1,nx)
ts_Bp=       fltarr(n1+1,nx)
ts_alp_m=    fltarr(n1+1,nx)
ts_Bzmod=    fltarr(n1+1,nx) 
ts_h=        fltarr(n1+1,nx)
ts_om=       fltarr(n1+1,nx)
ts_G=        fltarr(n1+1,nx)
ts_l=        fltarr(n1+1,nx)
ts_v=        fltarr(n1+1,nx)
ts_etat=     fltarr(n1+1,nx)
ts_tau=      fltarr(n1+1,nx)
ts_alp_k=    fltarr(n1+1,nx)
ts_Uz=       fltarr(n1+1,nx)
ts_Ur=       fltarr(n1+1,nx)
ts_n=        fltarr(n1+1,nx)
ts_Beq=      fltarr(n1+1,nx)
ts_rmax=     fltarr(n1+1) 
ts_delta_r=  fltarr(n1+1) 
ts_alp=      fltarr(n1+1,nx)
;read initial values
readf,11,r,h_init,om_init,Uz_init,Ur_init,l_init,v_init,etat_init,tau_init,alp_k_init,n_init,Beq_init,Br_init,Bp_init,Bzmod_init,alp_init
close,11
if (nvar eq 2 or nvar eq 4) then begin	;no dynamical quenching
  ;read final timestep values
  readf,12,   t,   Br,   Bp,            Bzmod,   h,   om,   G,   l,   v,   etat,   tau,   alp_k,   Uz,   Ur,   n,   Beq
  close,12
  ;read full time arrays
  readf,13,ts_t,ts_Br,ts_Bp,         ts_Bzmod,ts_h,ts_om,ts_G,ts_l,ts_v,ts_etat,ts_tau,ts_alp_k,ts_Uz,ts_Ur,ts_n,ts_Beq,ts_rmax,ts_delta_r,ts_alp
  close,13
  alp_m=0.*alp_k
endif else begin
  ;read final timestep values
  readf,12,   t,   Br,   Bp,   alp_m,   Bzmod,   h,   om,   G,   l,   v,   etat,   tau,   alp_k,   Uz,   Ur,   n,   Beq
  close,12
  ;read full time arrays
  readf,13,ts_t,ts_Br,ts_Bp,ts_alp_m,ts_Bzmod,ts_h,ts_om,ts_G,ts_l,ts_v,ts_etat,ts_tau,ts_alp_k,ts_Uz,ts_Ur,ts_n,ts_Beq,ts_rmax,ts_delta_r,ts_alp
  close,13
endelse
;alp=alp_k +alp_m
;ts_alp=ts_alp_k +ts_alp_m
cd,programdirec
end
