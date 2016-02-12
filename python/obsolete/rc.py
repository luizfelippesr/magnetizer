""" this code compares output from multiple simulation runs """
import sys
import numpy as np
import os

nruns=2

s0=[]
s0.append('CE_006')
s0.append('CE_007')

s1='/Users/luke/fortran_pde/1D/telegraph/r_noz/'
s2='_00000001'

simulationoutputdirec=[]
simulationoutputdirec.append(s1+'output/'+s0[0])
simulationoutputdirec.append(s1+'output/'+s0[1])

t0_Gyr=        np.zeros(nruns) 
t0_kpcskm=     np.zeros(nruns) 
h0_kpc=        np.zeros(nruns) 
etat0_cm2s=    np.zeros(nruns) 
n0_cm3=        np.zeros(nruns) 
B0_mkG=        np.zeros(nruns) 

nvar=          np.zeros(nruns) 
dt=            np.zeros(nruns) 
n1=            np.zeros(nruns) 
n2=            np.zeros(nruns) 
dx=            np.zeros(nruns) 
nxphys=        np.zeros(nruns) 
nxghost=       np.zeros(nruns) 
nx=            np.zeros(nruns) 

l_sol_kpc=     np.zeros(nruns) 
r_l_kpc=       np.zeros(nruns) 
v_sol_kms=     np.zeros(nruns) 
r_v_kpc=       np.zeros(nruns) 
n_sol_cm3=     np.zeros(nruns) 
r_n_kpc=       np.zeros(nruns) 
Uz_sol_kms=    np.zeros(nruns) 
r_Uz_kpc=      np.zeros(nruns) 
h_sol_kpc=     np.zeros(nruns) 
r_h_kpc=       np.zeros(nruns) 
Uphi_sol_kms=  np.zeros(nruns) 
r_om_kpc=      np.zeros(nruns) 
R_kappa=       np.zeros(nruns) 

r_in=          np.zeros(nruns) 
r_disk_kpc=    np.zeros(nruns) 
r_sol_kpc=     np.zeros(nruns) 
r1_kpc=        np.zeros(nruns) 
ctau=          np.zeros(nruns) 
nn=            np.zeros(nruns) 
lam=           np.zeros(nruns) 
C_alp=         np.zeros(nruns) 
alpceil=       np.zeros(nruns) 
Rm_inv=        np.zeros(nruns) 

etat_sol_cm2s= np.zeros(nruns) 
td_sol_Gyr=    np.zeros(nruns) 
tau_sol_Gyr=   np.zeros(nruns) 
etat_sol=      np.zeros(nruns) 
td_sol=        np.zeros(nruns) 
tau_sol=       np.zeros(nruns) 

simulationoutputdirec[0]= s1+'output/'+s0[0]
simulationoutputdirec[1]= s1+'output/'+s0[1]
programdirec= s1+'program/'
print 'simulationoutputdirec[0]=', simulationoutputdirec[0]
print 'simulationoutputdirec[0]=', simulationoutputdirec[1]
print 'programdirec=', programdirec

for i in range(nruns):
  os.chdir(simulationoutputdirec[i])
  print 'current working directory=',os.getcwd()
  
  paramf=open('param'+s2+'.out','r')
  words=paramf.readline().split()	#reads line1 into string list
  t0_Gyr       [i]= float(words[ 0])
  t0_kpcskm    [i]= float(words[ 1])
  h0_kpc       [i]= float(words[ 2])
  etat0_cm2s   [i]= float(words[ 3])
  n0_cm3       [i]= float(words[ 4])
  B0_mkG       [i]= float(words[ 5])
  words=paramf.readline().split()	#reads line2 into string list
  nvar         [i]= float(words[ 0])
  dt           [i]= float(words[ 1])
  n1           [i]= int  (words[ 2])
  n2           [i]= int  (words[ 3])
  dx           [i]= float(words[ 4])
  nxphys       [i]= int  (words[ 5])
  nxghost      [i]= int  (words[ 6])
  nx           [i]= int  (words[ 7])
  words=paramf.readline().split()	#reads line3 into string list
  l_sol_kpc    [i]= float(words[ 0])
  r_l_kpc      [i]= float(words[ 1])
  v_sol_kms    [i]= float(words[ 2])
  r_v_kpc      [i]= float(words[ 3])
  n_sol_cm3    [i]= float(words[ 4])
  r_n_kpc      [i]= float(words[ 5])
  Uz_sol_kms   [i]= float(words[ 6])
  r_Uz_kpc     [i]= float(words[ 7])
  h_sol_kpc    [i]= float(words[ 8])
  r_h_kpc      [i]= float(words[ 9])
  Uphi_sol_kms [i]= float(words[10])
  r_om_kpc     [i]= float(words[11])
  R_kappa      [i]= float(words[12])
  words=paramf.readline().split()	#reads line4 into string list
  r_in         [i]= float(words[ 0])
  r_disk_kpc   [i]= float(words[ 1])
  r_sol_kpc    [i]= float(words[ 2])
  r1_kpc       [i]= float(words[ 3])
  ctau         [i]= float(words[ 4])
  nn           [i]= int  (words[ 5])
  lam          [i]= float(words[ 6])
  C_alp        [i]= float(words[ 7])
  alpceil      [i]= float(words[ 8])
  Rm_inv       [i]= float(words[ 9])
  words=paramf.readline().split()	#reads line5 into string list
  etat_sol_cm2s[i]= float(words[ 0])
  td_sol_Gyr   [i]= float(words[ 1])
  tau_sol_Gyr  [i]= float(words[ 2])
  etat_sol     [i]= float(words[ 3])
  td_sol       [i]= float(words[ 4])
  tau_sol      [i]= float(words[ 5])
  
for i in range(nruns-1):
  if nx[i] != nx[i+1]:
    print ''
    print '***WARNING: grids for runs ',i,' and',i+1,' have different sizes: ',nx[i],' and ',nx[i+1]

print ''
for i in range(nruns):
  print 'RUN #',i,' UNITS'
  print '   ','t0_Gyr=%f, t0_kpcskm=%f, h0_kpc=%f, etat0_cm2s=%f, n0_cm3=%f, B0_mkG %f' \
         %(t0_Gyr[i],t0_kpcskm[i],h0_kpc[i],etat0_cm2s[i],n0_cm3[i],B0_mkG[i])
for i in range(nruns):
  print 'RUN #',i,' NUMERICS:'
  print '   ','nvar=%d, dt=%f, n1=%d, n2=%d, dx=%f, nxphys=%d, nxghost=%d, nx=%d' \
         %(nvar[i],dt[i],n1[i],n2[i],dx[i],nxphys[i],nxghost[i],nx[i])
print ''
for i in range(nruns):
  print 'RUN #',i,' INPUT PARAMETERS:'
  print '   ','l_sol_kpc=%f, r_l_kpc=%f, v_sol_kms=%f, r_v_kpc=%f, n_sol_cm3=%f, r_n_kpc=%f,' \
        %(l_sol_kpc[i],r_l_kpc[i],v_sol_kms[i],r_v_kpc[i],n_sol_cm3[i],r_n_kpc[i])
  print '   ','Uz_sol_kms=%f, r_Uz_kpc=%f, h_sol_kpc=%f, r_h_kpc=%f, Uphi_sol_kms=%f, r_om_kpc=%f, R_kappa=%f' \
        %(Uz_sol_kms[i],r_Uz_kpc[i],h_sol_kpc[i],r_h_kpc[i],Uphi_sol_kms[i],r_om_kpc[i],R_kappa[i])
print ''
for i in range(nruns):
  print 'RUN #',i,' OTHER FREE PARAMETERS:'
  print '   ','r_in=%f, r_disk_kpc=%f, r_sol_kpc=%f, r1_kpc=%f, ctau=%f, nn=%d, lam=%f, C_alp=%f, alpceil=%f, Rm_inv=%f' \
         %(r_in[i],r_disk_kpc[i],r_sol_kpc[i],r1_kpc[i],ctau[i],nn[i],lam[i],C_alp[i],alpceil[i],Rm_inv[i])
print ''
for i in range(nruns):
  print 'RUN #',i,' IMPORTANT CALCULATED PARAMETERS'
  print '   ','etat_sol_cm2s=%f, td_sol_Gyr=%f, tau_sol_Gyr=%f, etat_sol=%f, td_sol=%f, tau_sol=%f' \
         %(etat_sol_cm2s[i],td_sol_Gyr[i],tau_sol_Gyr[i],etat_sol[i],td_sol[i],tau_sol[i])
print ''
  
rows=n1[0]+1
cols=nx[0]
rows=int(rows)
cols=int(cols)

def readfn(file,numel):
  words=file.readline().split()
  par=np.zeros(numel)
  for j in range(len(words)):
    par[j]=float(words[j])
  return par
 
r         =np.zeros((nruns,cols))
h_init    =np.zeros((nruns,cols))
om_init   =np.zeros((nruns,cols))
G_init    =np.zeros((nruns,cols))
Uz_init   =np.zeros((nruns,cols))
Ur_init   =np.zeros((nruns,cols))
l_init    =np.zeros((nruns,cols))
v_init    =np.zeros((nruns,cols))
etat_init =np.zeros((nruns,cols))
tau_init  =np.zeros((nruns,cols))
alpk_init =np.zeros((nruns,cols))
n_init    =np.zeros((nruns,cols))
Beq_init  =np.zeros((nruns,cols))
Br_init   =np.zeros((nruns,cols))
Bp_init   =np.zeros((nruns,cols))
Bzmod_init=np.zeros((nruns,cols))

for i in range(nruns):
  os.chdir(simulationoutputdirec[i])
  print 'current working directory=',os.getcwd()
  
  isf=open('init'+s2+'.out','r')
  
  r         [i,:]= readfn(isf,cols)
  h_init    [i,:]= readfn(isf,cols)
  om_init   [i,:]= readfn(isf,cols)
  G_init    [i,:]= readfn(isf,cols)
  Uz_init   [i,:]= readfn(isf,cols)
  Ur_init   [i,:]= readfn(isf,cols) 
  l_init    [i,:]= readfn(isf,cols)
  v_init    [i,:]= readfn(isf,cols)
  etat_init [i,:]= readfn(isf,cols)
  tau_init  [i,:]= readfn(isf,cols)
  alpk_init [i,:]= readfn(isf,cols)
  n_init    [i,:]= readfn(isf,cols)
  Beq_init  [i,:]= readfn(isf,cols)
  Br_init   [i,:]= readfn(isf,cols)
  Bp_init   [i,:]= readfn(isf,cols)
  Bzmod_init[i,:]= readfn(isf,cols)

for i in range(nruns):
  print 'Run #',i,' r=          ',r         [i]
for i in range(nruns):
  print 'Run #',i,' h_init=     ',h_init    [i]
for i in range(nruns):
  print 'Run #',i,' om_init=    ',om_init   [i]
for i in range(nruns):
  print 'Run #',i,' G_init=     ',G_init    [i]
for i in range(nruns):
  print 'Run #',i,' Uz_init=    ',Uz_init   [i]
for i in range(nruns):
  print 'Run #',i,' Ur_init=    ',Ur_init   [i]
for i in range(nruns):
  print 'Run #',i,' l_init=     ',l_init    [i]
for i in range(nruns):
  print 'Run #',i,' v_init=     ',v_init    [i]
for i in range(nruns):
  print 'Run #',i,' etat_init=  ',etat_init [i]
for i in range(nruns):
  print 'Run #',i,' tau_init=   ',tau_init  [i] 
for i in range(nruns):
  print 'Run #',i,' alpk_init=  ',alpk_init [i]
for i in range(nruns):
  print 'Run #',i,' n_init=     ',n_init    [i]
for i in range(nruns):
  print 'Run #',i,' Beq_init=   ',Beq_init  [i]
for i in range(nruns):
  print 'Run #',i,' Br_init=    ',Br_init   [i]
for i in range(nruns):
  print 'Run #',i,' Bp_init=    ',Bp_init   [i]
for i in range(nruns):
  print 'Run #',i,' Bzmod_init= ',Bzmod_init[i]

ts_t=     np.zeros((nruns,rows))
ts_Br=    np.zeros((nruns,rows,cols))
ts_Bp=    np.zeros((nruns,rows,cols))
ts_alpm=  np.zeros((nruns,rows,cols))
ts_Bzmod= np.zeros((nruns,rows,cols))
ts_h=     np.zeros((nruns,rows,cols))
ts_Om=    np.zeros((nruns,rows,cols))
ts_G=     np.zeros((nruns,rows,cols))
ts_l=     np.zeros((nruns,rows,cols))
ts_v=     np.zeros((nruns,rows,cols))
ts_etat=  np.zeros((nruns,rows,cols))
ts_tau=   np.zeros((nruns,rows,cols))
ts_alpk=  np.zeros((nruns,rows,cols))
ts_Uz=    np.zeros((nruns,rows,cols))
ts_Ur=    np.zeros((nruns,rows,cols))
ts_n=     np.zeros((nruns,rows,cols))
ts_Beq=   np.zeros((nruns,rows,cols))

for i in range(nruns):
  os.chdir(simulationoutputdirec[i])
  print 'current working directory=',os.getcwd()
  
  tsf=open('ts'+s2+'.out','r')
  
  ts_t     [i,:]=  readfn(tsf,rows)
  ts_Br    [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T #read as 1-D array, convert to 2-D array, transpose
  ts_Bp    [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  if (nvar[i]!=2 or nvar[i]!=4):	#condition that dynamical quenching was used so alpha_m was calculated
    ts_alpm[i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Bzmod [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_h     [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Om    [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_G     [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_l     [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_v     [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_etat  [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_tau   [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_alpk  [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Uz    [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Ur    [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_n     [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Beq   [i,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T

for i in range(nruns):
  print 'Run #',i,' t:'
  print 'ts_t[  1]=', ts_t[i,  0]
  print 'ts_t[421]=', ts_t[i,420]
  print 'ts_t[100]=', ts_t[i, 99]
  print ''
  print 'Run #',i,' Br(t,r):'
  print 'ts_Br [  1, 1]=',ts_Br [i,  0, 0]
  print 'ts_Br [  1,37]=',ts_Br [i,  0,36]
  print 'ts_Br [421, 1]=',ts_Br [i,420, 0]
  print 'ts_Br [421,37]=',ts_Br [i,420,36]
  print 'ts_Br [100,10]=',ts_Br [i, 99, 9]
  print ''
  print 'Run #',i,' Beq(t,r):'
  print 'ts_Beq[  1, 1]=',ts_Beq[i,  0, 0]
  print 'ts_Beq[  1,37]=',ts_Beq[i,  0,36]
  print 'ts_Beq[421, 1]=',ts_Beq[i,420, 0]
  print 'ts_Beq[421,37]=',ts_Beq[i,420,36]
  print 'ts_Beq[100,10]=',ts_Beq[i, 99, 9]

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
