import sys
import numpy as np
import os

#this code compares output from multiple simulation runs

ngals=2 #number of galaxies

s0='CE_006' #name of run
s1='/Users/luke/fortran_pde/1D/telegraph/r_noz/'
simulationoutputdirec= s1+'output/'+s0
programdirec= s1+'program/'+ s0
print 'simulationoutputdirec=', simulationoutputdirec
print 'programdirec=', programdirec
os.chdir(simulationoutputdirec)
print 'current working directory=',os.getcwd()

t0_Gyr=        np.zeros(ngals) 
t0_kpcskm=     np.zeros(ngals) 
h0_kpc=        np.zeros(ngals) 
etat0_cm2s=    np.zeros(ngals) 
n0_cm3=        np.zeros(ngals) 
B0_mkG=        np.zeros(ngals) 

nvar=          np.zeros(ngals) 
dt=            np.zeros(ngals) 
n1=            np.zeros(ngals) 
n2=            np.zeros(ngals) 
dx=            np.zeros(ngals) 
nxphys=        np.zeros(ngals) 
nxghost=       np.zeros(ngals) 
nx=            np.zeros(ngals) 

l_sol_kpc=     np.zeros(ngals) 
r_l_kpc=       np.zeros(ngals) 
v_sol_kms=     np.zeros(ngals) 
r_v_kpc=       np.zeros(ngals) 
n_sol_cm3=     np.zeros(ngals) 
r_n_kpc=       np.zeros(ngals) 
Uz_sol_kms=    np.zeros(ngals) 
r_Uz_kpc=      np.zeros(ngals) 
h_sol_kpc=     np.zeros(ngals) 
r_h_kpc=       np.zeros(ngals) 
Uphi_sol_kms=  np.zeros(ngals) 
r_om_kpc=      np.zeros(ngals) 
R_kappa=       np.zeros(ngals) 

r_in=          np.zeros(ngals) 
r_disk_kpc=    np.zeros(ngals) 
r_sol_kpc=     np.zeros(ngals) 
r1_kpc=        np.zeros(ngals) 
ctau=          np.zeros(ngals) 
nn=            np.zeros(ngals) 
lam=           np.zeros(ngals) 
C_alp=         np.zeros(ngals) 
alpceil=       np.zeros(ngals) 
Rm_inv=        np.zeros(ngals) 

etat_sol_cm2s= np.zeros(ngals) 
td_sol_Gyr=    np.zeros(ngals) 
tau_sol_Gyr=   np.zeros(ngals) 
etat_sol=      np.zeros(ngals) 
td_sol=        np.zeros(ngals) 
tau_sol=       np.zeros(ngals) 

def s2(j):
  id='_'+'%08d' %(j+1) #id number of first galaxy to compare
  return id

os.chdir(simulationoutputdirec)
print 'current working directory=',os.getcwd()
for j in range(ngals):
  paramf=open('param'+s2(j)+'.out','r')
  words=paramf.readline().split()	#reads line1 into string list
  t0_Gyr       [j]= float(words[ 0])
  t0_kpcskm    [j]= float(words[ 1])
  h0_kpc       [j]= float(words[ 2])
  etat0_cm2s   [j]= float(words[ 3])
  n0_cm3       [j]= float(words[ 4])
  B0_mkG       [j]= float(words[ 5])
  words=paramf.readline().split()	#reads line2 into string list
  nvar         [j]= float(words[ 0])
  dt           [j]= float(words[ 1])
  n1           [j]= int  (words[ 2])
  n2           [j]= int  (words[ 3])
  dx           [j]= float(words[ 4])
  nxphys       [j]= int  (words[ 5])
  nxghost      [j]= int  (words[ 6])
  nx           [j]= int  (words[ 7])
  words=paramf.readline().split()	#reads line3 into string list
  l_sol_kpc    [j]= float(words[ 0])
  r_l_kpc      [j]= float(words[ 1])
  v_sol_kms    [j]= float(words[ 2])
  r_v_kpc      [j]= float(words[ 3])
  n_sol_cm3    [j]= float(words[ 4])
  r_n_kpc      [j]= float(words[ 5])
  Uz_sol_kms   [j]= float(words[ 6])
  r_Uz_kpc     [j]= float(words[ 7])
  h_sol_kpc    [j]= float(words[ 8])
  r_h_kpc      [j]= float(words[ 9])
  Uphi_sol_kms [j]= float(words[10])
  r_om_kpc     [j]= float(words[11])
  R_kappa      [j]= float(words[12])
  words=paramf.readline().split()	#reads line4 into string list
  r_in         [j]= float(words[ 0])
  r_disk_kpc   [j]= float(words[ 1])
  r_sol_kpc    [j]= float(words[ 2])
  r1_kpc       [j]= float(words[ 3])
  ctau         [j]= float(words[ 4])
  nn           [j]= int  (words[ 5])
  lam          [j]= float(words[ 6])
  C_alp        [j]= float(words[ 7])
  alpceil      [j]= float(words[ 8])
  Rm_inv       [j]= float(words[ 9])
  words=paramf.readline().split()	#reads line5 into string list
  etat_sol_cm2s[j]= float(words[ 0])
  td_sol_Gyr   [j]= float(words[ 1])
  tau_sol_Gyr  [j]= float(words[ 2])
  etat_sol     [j]= float(words[ 3])
  td_sol       [j]= float(words[ 4])
  tau_sol      [j]= float(words[ 5])

for j in range(ngals):
  if nx[j] != nx[0]:
    print ''
    print '***WARNING: grid for galaxy ',s2(j),' has a different size then for run 0, galaxy _00000001'

print 'UNITS'
for j in range(ngals):
  print '  Galaxy #',s2(j)
  print '   ','t0_Gyr=%f, t0_kpcskm=%f, h0_kpc=%f, etat0_cm2s=%f, n0_cm3=%f, B0_mkG %f' \
         %(t0_Gyr[j],t0_kpcskm[j],h0_kpc[j],etat0_cm2s[j],n0_cm3[j],B0_mkG[j])
print 'NUMERICS:'
for j in range(ngals):
  print '  Galaxy #',s2(j)
  print '   ','nvar=%d, dt=%f, n1=%d, n2=%d, dx=%f, nxphys=%d, nxghost=%d, nx=%d' \
         %(nvar[j],dt[j],n1[j],n2[j],dx[j],nxphys[j],nxghost[j],nx[j])
print ''
print 'INPUT PARAMETERS:'
for j in range(ngals):
  print '  Galaxy #',s2(j)
  print '   ','l_sol_kpc=%f, r_l_kpc=%f, v_sol_kms=%f, r_v_kpc=%f, n_sol_cm3=%f, r_n_kpc=%f,' \
        %(l_sol_kpc[j],r_l_kpc[j],v_sol_kms[j],r_v_kpc[j],n_sol_cm3[j],r_n_kpc[j])
  print '   ','Uz_sol_kms=%f, r_Uz_kpc=%f, h_sol_kpc=%f, r_h_kpc=%f, Uphi_sol_kms=%f, r_om_kpc=%f, R_kappa=%f' \
        %(Uz_sol_kms[j],r_Uz_kpc[j],h_sol_kpc[j],r_h_kpc[j],Uphi_sol_kms[j],r_om_kpc[j],R_kappa[j])
print ''
print 'OTHER FREE PARAMETERS:'
for j in range(ngals):
  print '  Galaxy #',s2(j)
  print '   ','r_in=%f, r_disk_kpc=%f, r_sol_kpc=%f, r1_kpc=%f, ctau=%f, nn=%d, lam=%f, C_alp=%f, alpceil=%f, Rm_inv=%f' \
         %(r_in[j],r_disk_kpc[j],r_sol_kpc[j],r1_kpc[j],ctau[j],nn[j],lam[j],C_alp[j],alpceil[j],Rm_inv[j])
print ''
print 'IMPORTANTALCULATED PARATERS'
for j in range(ngals):
  print '  Galaxy #',s2(j)
  print '   ','etat_sol_cm2s=%f, td_sol_Gyr=%f, tau_sol_Gyr=%f, etat_sol=%f, td_sol=%f, tau_sol=%f' \
         %(etat_sol_cm2s[j],td_sol_Gyr[j],tau_sol_Gyr[j],etat_sol[j],td_sol[j],tau_sol[j])
print ''
  
rows=n1[0]+1
cols=nx[0]
rows=int(rows)
cols=int(cols)
print rows,cols

def readfn(file,numel):
  words=file.readline().split()
  par=np.zeros(numel)
  for j in range(len(words)):
    par[j]=float(words[j])
  return par
 
r         =np.zeros((ngals,cols))
h_init    =np.zeros((ngals,cols))
om_init   =np.zeros((ngals,cols))
G_init    =np.zeros((ngals,cols))
Uz_init   =np.zeros((ngals,cols))
Ur_init   =np.zeros((ngals,cols))
l_init    =np.zeros((ngals,cols))
v_init    =np.zeros((ngals,cols))
etat_init =np.zeros((ngals,cols))
tau_init  =np.zeros((ngals,cols))
alpk_init =np.zeros((ngals,cols))
n_init    =np.zeros((ngals,cols))
Beq_init  =np.zeros((ngals,cols))
Br_init   =np.zeros((ngals,cols))
Bp_init   =np.zeros((ngals,cols))
Bzmod_init=np.zeros((ngals,cols))

os.chdir(simulationoutputdirec)
print 'current working directory=',os.getcwd()

for j in range(ngals):
  isf=open('init'+s2(j)+'.out','r')
  
  r         [j,:]= readfn(isf,cols)
  h_init    [j,:]= readfn(isf,cols)
  om_init   [j,:]= readfn(isf,cols)
  G_init    [j,:]= readfn(isf,cols)
  Uz_init   [j,:]= readfn(isf,cols)
  Ur_init   [j,:]= readfn(isf,cols) 
  l_init    [j,:]= readfn(isf,cols)
  v_init    [j,:]= readfn(isf,cols)
  etat_init [j,:]= readfn(isf,cols)
  tau_init  [j,:]= readfn(isf,cols)
  alpk_init [j,:]= readfn(isf,cols)
  n_init    [j,:]= readfn(isf,cols)
  Beq_init  [j,:]= readfn(isf,cols)
  Br_init   [j,:]= readfn(isf,cols)
  Bp_init   [j,:]= readfn(isf,cols)
  Bzmod_init[j,:]= readfn(isf,cols)

for j in range(ngals):
    print ' Galaxy',s2(j),' r=          ',r         [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' h_init=     ',h_init    [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' om_init=    ',om_init   [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' G_init=     ',G_init    [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' Uz_init=    ',Uz_init   [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' Ur_init=    ',Ur_init   [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' l_init=     ',l_init    [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' v_init=     ',v_init    [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' etat_init=  ',etat_init [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' tau_init=   ',tau_init  [j] 
for j in range(ngals):
    print ' Galaxy',s2(j),' alpk_init=  ',alpk_init [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' n_init=     ',n_init    [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' Beq_init=   ',Beq_init  [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' Br_init=    ',Br_init   [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' Bp_init=    ',Bp_init   [j]
for j in range(ngals):
    print ' Galaxy',s2(j),' Bzmod_init= ',Bzmod_init[j]

ts_t=     np.zeros((ngals,rows))
ts_Br=    np.zeros((ngals,rows,cols))
ts_Bp=    np.zeros((ngals,rows,cols))
ts_alpm=  np.zeros((ngals,rows,cols))
ts_Bzmod= np.zeros((ngals,rows,cols))
ts_h=     np.zeros((ngals,rows,cols))
ts_Om=    np.zeros((ngals,rows,cols))
ts_G=     np.zeros((ngals,rows,cols))
ts_l=     np.zeros((ngals,rows,cols))
ts_v=     np.zeros((ngals,rows,cols))
ts_etat=  np.zeros((ngals,rows,cols))
ts_tau=   np.zeros((ngals,rows,cols))
ts_alpk=  np.zeros((ngals,rows,cols))
ts_Uz=    np.zeros((ngals,rows,cols))
ts_Ur=    np.zeros((ngals,rows,cols))
ts_n=     np.zeros((ngals,rows,cols))
ts_Beq=   np.zeros((ngals,rows,cols))

os.chdir(simulationoutputdirec)
print 'current working directory=',os.getcwd()

for j in range(ngals):
  tsf=open('ts'+s2(j)+'.out','r')

  ts_t     [j,:]=  readfn(tsf,rows)
  ts_Br    [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T #read as 1-D array, convert to 2-D array, transpose
  ts_Bp    [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  if (nvar [j]!=2 or nvar[j]!=4):	#condition that dynamical quenching was used so alpha_m was calculated
    ts_alpm[j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Bzmod [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_h     [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Om    [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_G     [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_l     [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_v     [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_etat  [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_tau   [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_alpk  [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Uz    [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Ur    [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_n     [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
  ts_Beq   [j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T

for j in range(ngals):
  print ' Galaxy',s2(j),' t:'
  print 'ts_t[  1]=', ts_t[j,  0]
  print 'ts_t[421]=', ts_t[j,420]
  print 'ts_t[100]=', ts_t[j, 99]
  print ''
  print ' Galaxy',s2(j),' Br(t,r):'
  print 'ts_Br [  1, 1]=',ts_Br [j,  0, 0]
  print 'ts_Br [  1,37]=',ts_Br [j,  0,36]
  print 'ts_Br [421, 1]=',ts_Br [j,420, 0]
  print 'ts_Br [421,37]=',ts_Br [j,420,36]
  print 'ts_Br [100,10]=',ts_Br [j, 99, 9]
  print ''
  print ' Galaxy',s2(j),' Beq(t,r):'
  print 'ts_Beq[  1, 1]=',ts_Beq[j,  0, 0]
  print 'ts_Beq[  1,37]=',ts_Beq[j,  0,36]
  print 'ts_Beq[421, 1]=',ts_Beq[j,420, 0]
  print 'ts_Beq[421,37]=',ts_Beq[j,420,36]
  print 'ts_Beq[100,10]=',ts_Beq[j, 99, 9]
  print ''

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
