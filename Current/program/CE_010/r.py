import sys
import numpy as np
import os

s0='CE_010'
s1='/Users/luke/fortran_pde/1D/telegraph/r_noz/'
s2='_00000001'
simulationoutputdirec= s1+'output/'+s0
programdirec= s1+'program/'+ s0
print 'simulationoutputdirec=', simulationoutputdirec
print 'programdirec=', programdirec
os.chdir(simulationoutputdirec)
print 'current working directory=',os.getcwd()

paramf=open('param'+s2+'.out','r')
words=paramf.readline().split()	#reads line1 into string list
t0_Gyr=        float(words[ 0])
t0_kpcskm=     float(words[ 1])
h0_kpc=        float(words[ 2])
etat0_cm2s=    float(words[ 3])
n0_cm3=        float(words[ 4])
B0_mkG=        float(words[ 5])
words=paramf.readline().split()	#reads line2 into string list
nvar=          float(words[ 0])
dt=            float(words[ 1])
n1=            int  (words[ 2])
n2=            int  (words[ 3])
dx=            float(words[ 4])
nxphys=        int  (words[ 5])
nxghost=       int  (words[ 6])
nx=            int  (words[ 7])
words=paramf.readline().split()	#reads line3 into string list
l_sol_kpc=     float(words[ 0])
r_l_kpc=       float(words[ 1])
v_sol_kms=     float(words[ 2])
r_v_kpc=       float(words[ 3])
n_sol_cm3=     float(words[ 4])
r_n_kpc=       float(words[ 5])
Uz_sol_kms=    float(words[ 6])
r_Uz_kpc=      float(words[ 7])
h_sol_kpc=     float(words[ 8])
r_h_kpc=       float(words[ 9])
Uphi_sol_kms=  float(words[10])
r_om_kpc=      float(words[11])
R_kappa=       float(words[12])
words=paramf.readline().split()	#reads line4 into string list
r_in=          float(words[ 0])
r_disk_kpc=    float(words[ 1])
r_sol_kpc=     float(words[ 2])
r1_kpc=        float(words[ 3])
ctau=          float(words[ 4])
nn=            int  (words[ 5])
lam=           float(words[ 6])
C_alp=         float(words[ 7])
alpceil=       float(words[ 8])
Rm_inv=        float(words[ 9])
words=paramf.readline().split()	#reads line5 into string list
etat_sol_cm2s= float(words[ 0])
td_sol_Gyr=    float(words[ 1])
tau_sol_Gyr=   float(words[ 2])
etat_sol=      float(words[ 3])
td_sol=        float(words[ 4])
tau_sol=       float(words[ 5])

print 'UNITS'
print 't0_Gyr=%f, t0_kpcskm=%f, h0_kpc=%f, etat0_cm2s=%f, n0_cm3=%f, B0_mkG %f' \
       %(t0_Gyr,t0_kpcskm,h0_kpc,etat0_cm2s,n0_cm3,B0_mkG)
print ''
print 'NUMERICS:'
print 'nvar=%d, dt=%f, n1=%d, n2=%d, dx=%f, nxphys=%d, nxghost=%d, nx=%d' \
       %(nvar,dt,n1,n2,dx,nxphys,nxghost,nx)
print ''
print 'INPUT PARAMETERS:'
print 'l_sol_kpc=%f, r_l_kpc=%f, v_sol_kms=%f, r_v_kpc=%f, n_sol_cm3=%f, r_n_kpc=%f, \
       Uz_sol_kms=%f, r_Uz_kpc=%f, h_sol_kpc=%f, r_h_kpc=%f, Uphi_sol_kms=%f, r_om_kpc=%f, R_kappa=%f' \
       %(l_sol_kpc,r_l_kpc,v_sol_kms,r_v_kpc,n_sol_cm3,r_n_kpc,Uz_sol_kms,r_Uz_kpc,h_sol_kpc,r_h_kpc,Uphi_sol_kms,r_om_kpc,R_kappa)
print ''
print 'OTHER FREE PARAMETERS:'
print 'r_in=%f, r_disk_kpc=%f, r_sol_kpc=%f, r1_kpc=%f, ctau=%f, nn=%d, lam=%f, C_alp=%f, alpceil=%f, Rm_inv=%f' \
       %(r_in,r_disk_kpc,r_sol_kpc,r1_kpc,ctau,nn,lam,C_alp,alpceil,Rm_inv)
print ''
print 'IMPORTANT CALCULATED PARAMETERS'
print 'etat_sol_cm2s=%f, td_sol_Gyr=%f, tau_sol_Gyr=%f, etat_sol=%f, td_sol=%f, tau_sol=%f' \
       %(etat_sol_cm2s,td_sol_Gyr,tau_sol_Gyr,etat_sol,td_sol,tau_sol)

rows=n1+1
cols=nx

#openr,11,strjoin(['init' ,s2,'.out'])

isf=open('init'+s2+'.out','r')

def readfn(file,numel):
  words=file.readline().split()
  par=np.zeros(numel)
  for i in range(len(words)):
    par[i]=float(words[i])
  return par

r=          readfn(isf,cols)
h_init=     readfn(isf,cols)
om_init=    readfn(isf,cols)
G_init=     readfn(isf,cols)
Uz_init=    readfn(isf,cols)
Ur_init=    readfn(isf,cols) 
l_init=     readfn(isf,cols)
v_init=     readfn(isf,cols)
etat_init=  readfn(isf,cols)
tau_init=   readfn(isf,cols)
alpk_init=  readfn(isf,cols)
n_init=     readfn(isf,cols)
Beq_init=   readfn(isf,cols)
Br_init=    readfn(isf,cols)
Bp_init=    readfn(isf,cols)
Bzmod_init= readfn(isf,cols)

print 'r=          ',r
print 'h_init=     ',h_init
print 'om_init=    ',om_init  
print 'G_init=     ',G_init  
print 'Uz_init=    ',Uz_init
print 'Ur_init=    ',Ur_init
print 'l_init=     ',l_init
print 'v_init=     ',v_init
print 'etat_init=  ',etat_init 
print 'tau_init=   ',tau_init  
print 'alpk_init=  ',alpk_init
print 'n_init=     ',n_init
print 'Beq_init=   ',Beq_init
print 'Br_init=    ',Br_init
print 'Bp_init=    ',Bp_init
print 'Bzmod_init= ',Bzmod_init


tsf=open('ts'+s2+'.out','r')

ts_t=     readfn(tsf,rows)
ts_Br=    readfn(tsf,rows*cols).reshape(cols,rows).T #read as 1-D array, convert to 2-D array, transpose

ts_Bp=    readfn(tsf,rows*cols).reshape(cols,rows).T
if (nvar!=2 or nvar!=4):	#condition that dynamical quenching was used so alpha_m was calculated
  ts_alpm=readfn(tsf,rows*cols).reshape(cols,rows).T
ts_Bzmod= readfn(tsf,rows*cols).reshape(cols,rows).T
ts_h=     readfn(tsf,rows*cols).reshape(cols,rows).T
ts_Om=    readfn(tsf,rows*cols).reshape(cols,rows).T
ts_G=     readfn(tsf,rows*cols).reshape(cols,rows).T
ts_l=     readfn(tsf,rows*cols).reshape(cols,rows).T
ts_v=     readfn(tsf,rows*cols).reshape(cols,rows).T
ts_etat=  readfn(tsf,rows*cols).reshape(cols,rows).T
ts_tau=   readfn(tsf,rows*cols).reshape(cols,rows).T
ts_alpk=  readfn(tsf,rows*cols).reshape(cols,rows).T
ts_Uz=    readfn(tsf,rows*cols).reshape(cols,rows).T
ts_Ur=    readfn(tsf,rows*cols).reshape(cols,rows).T
ts_n=     readfn(tsf,rows*cols).reshape(cols,rows).T
ts_Beq=   readfn(tsf,rows*cols).reshape(cols,rows).T
#
#
print 'ts_t[  1]=', ts_t[  0]
print 'ts_t[421]=', ts_t[420]
print 'ts_t[100]=', ts_t[ 99]
#
print 'ts_Br [  1, 1]=',ts_Br [  0, 0]
print 'ts_Br [  1,37]=',ts_Br [  0,36]
print 'ts_Br [421, 1]=',ts_Br [420, 0]
print 'ts_Br [421,37]=',ts_Br [420,36]
print 'ts_Br [100,10]=',ts_Br [ 99, 9]
#                                       #               
print 'ts_Beq[  1, 1]=',ts_Beq[  0, 0]
print 'ts_Beq[  1,37]=',ts_Beq[  0,36]
print 'ts_Beq[421, 1]=',ts_Beq[420, 0]
print 'ts_Beq[421,37]=',ts_Beq[420,36]
print 'ts_Beq[100,10]=',ts_Beq[ 99, 9]

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
