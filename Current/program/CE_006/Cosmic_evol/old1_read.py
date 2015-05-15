import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator
import os

s0='CE_006'
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
temp=np.zeros(cols)

isf=open('init'+s2+'.out','r')

words=isf.readline().split()	#r
for i in range(len(words)):
  temp[i]=float(words[i])
r=temp

words=isf.readline().split()	#h_init
for i in range(len(words)):
  temp[i]=float(words[i])
h_init=temp

words=isf.readline().split()	#om_init
for i in range(len(words)):
  temp[i]=float(words[i])
om_init=temp

words=isf.readline().split()	#G_init
for i in range(len(words)):
  temp[i]=float(words[i])
G_init=temp

words=isf.readline().split()	#l_init
for i in range(len(words)):
  temp[i]=float(words[i])
l_init=temp

words=isf.readline().split()	#v_init
for i in range(len(words)):
  temp[i]=float(words[i])
v_init=temp

words=isf.readline().split()	#etat_init
for i in range(len(words)):
  temp[i]=float(words[i])
etat_init=temp

words=isf.readline().split()	#tau_init
for i in range(len(words)):
  temp[i]=float(words[i])
tau_init=temp

words=isf.readline().split()	#alpk_init
for i in range(len(words)):
  temp[i]=float(words[i])
alpk_init=temp

words=isf.readline().split()	#Uz_init
for i in range(len(words)):
  temp[i]=float(words[i])
Uz_init=temp

words=isf.readline().split()	#Ur_init
for i in range(len(words)):
  temp[i]=float(words[i])
Ur_init=temp

words=isf.readline().split()	#alpm_init
for i in range(len(words)):
  temp[i]=float(words[i])
alpm_init=temp

words=isf.readline().split()	#n_init
for i in range(len(words)):
  temp[i]=float(words[i])
n_init=temp

words=isf.readline().split()	#Beq_init
for i in range(len(words)):
  temp[i]=float(words[i])
Beq_init=temp

words=isf.readline().split()	#Br_init
for i in range(len(words)):
  temp[i]=float(words[i])
Br_init=temp

words=isf.readline().split()	#Bp_init
for i in range(len(words)):
  temp[i]=float(words[i])
Bp_init=temp

words=isf.readline().split()	#Bzmod_init
for i in range(len(words)):
  temp[i]=float(words[i])
Bzmod_init=temp

temp=np.zeros(rows*cols)

tsf=open('ts'+s2+'.out','r')

words=tsf.readline().split()	#ts_t
for i in range(len(words)):
  temp[i]=float(words[i])
ts_t=temp

print 'ts_t[  0]=', ts_t[  0]
print 'ts_t[421]=', ts_t[420]
print 'ts_t[100]=', ts_t[ 99]

words=tsf.readline().split()	#ts_Br
for i in range(len(words)):
  temp[i]=float(words[i])
ts_Br=temp.reshape(cols,rows)

print 'ts_Br[  0, 0]=', ts_Br[ 0,  0]
print 'ts_Br[ 37, 0]=', ts_Br[36,  0]
print 'ts_Br[421, 0]=', ts_Br[ 0,420]
print 'ts_Br[421,37]=', ts_Br[36,420]
print 'ts_Br[100,10]=', ts_Br[ 9, 99]

words=tsf.readline().split()	#ts_Bp
for i in range(len(words)):
  temp[i]=float(words[i])
ts_Bp=    temp.reshape(cols,rows)

if (nvar!=2 or nvar!=4):	#condition that dynamical quenching was used so alpha_m was calculated
  words=tsf.readline().split()#ts_alpm
  for i in range(len(words)):
    temp[i]=float(words[i])
  ts_alpm=temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_Bzmod
for i in range(len(words)):
    temp[i]=float(words[i])
ts_Bzmod= temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_h
for i in range(len(words)):
    temp[i]=float(words[i])
ts_h=     temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_Om
for i in range(len(words)):
    temp[i]=float(words[i])
ts_Om=    temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_G
for i in range(len(words)):
    temp[i]=float(words[i])
ts_G=     temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_l
for i in range(len(words)):
    temp[i]=float(words[i])
ts_l=     temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_v
for i in range(len(words)):
    temp[i]=float(words[i])
ts_v=     temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_etat
for i in range(len(words)):
    temp[i]=float(words[i])
ts_etat=  temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_tau
for i in range(len(words)):
    temp[i]=float(words[i])
ts_tau=   temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_alpk
for i in range(len(words)):
    temp[i]=float(words[i])
ts_alpk=  temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_Uz
for i in range(len(words)):
    temp[i]=float(words[i])
ts_Uz=   temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_Ur
for i in range(len(words)):
    temp[i]=float(words[i])
ts_Ur=   temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_n
for i in range(len(words)):
    temp[i]=float(words[i])
ts_n=    temp.reshape(cols,rows)

words=tsf.readline().split()	#ts_Beq
for i in range(len(words)):
    temp[i]=float(words[i])
ts_Beq=  temp.reshape(cols,rows)

#readf,13,ts_t,ts_Br,ts_Bp,ts_alp_m,ts_Bzmod,ts_h,ts_om,ts_G,ts_l,ts_v,ts_etat,ts_tau,ts_alp_k,ts_Uz,ts_Ur,ts_n,ts_Beq

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
