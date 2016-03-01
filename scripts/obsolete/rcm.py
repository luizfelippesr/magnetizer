""" this code compares output from multiple simulation runs """
import sys
import numpy as np
import os

nruns=1 #number of runs
ngals=np.zeros(nruns) #number of galaxies in each run
ngals[0]=2 #number of galaxies in first run
#ngals[1]=5 #number of galaxies in second run
#ngals[2]=3 #number of galaxies in second run
ngmax=int(np.max(ngals))
print 'ngmax=',ngmax

s0=[]
s0.append('CE_006') #name of first run to compare
#s0.append('CE_009') #name of second run to compare

s1='/Users/luke/fortran_pde/1D/telegraph/r_noz/'

simulationoutputdirec=[]
simulationoutputdirec.append(s1+'output/'+s0[0])
if s0[-1] != s0[0]:
  simulationoutputdirec.append(s1+'output/'+s0[1])

t0_Gyr=        np.zeros((nruns,ngmax)) 
t0_kpcskm=     np.zeros((nruns,ngmax)) 
h0_kpc=        np.zeros((nruns,ngmax)) 
etat0_cm2s=    np.zeros((nruns,ngmax)) 
n0_cm3=        np.zeros((nruns,ngmax)) 
B0_mkG=        np.zeros((nruns,ngmax)) 

nvar=          np.zeros((nruns,ngmax)) 
dt=            np.zeros((nruns,ngmax)) 
n1=            np.zeros((nruns,ngmax)) 
n2=            np.zeros((nruns,ngmax)) 
dx=            np.zeros((nruns,ngmax)) 
nxphys=        np.zeros((nruns,ngmax)) 
nxghost=       np.zeros((nruns,ngmax)) 
nx=            np.zeros((nruns,ngmax)) 

l_sol_kpc=     np.zeros((nruns,ngmax)) 
r_l_kpc=       np.zeros((nruns,ngmax)) 
v_sol_kms=     np.zeros((nruns,ngmax)) 
r_v_kpc=       np.zeros((nruns,ngmax)) 
n_sol_cm3=     np.zeros((nruns,ngmax)) 
r_n_kpc=       np.zeros((nruns,ngmax)) 
Uz_sol_kms=    np.zeros((nruns,ngmax)) 
r_Uz_kpc=      np.zeros((nruns,ngmax)) 
h_sol_kpc=     np.zeros((nruns,ngmax)) 
r_h_kpc=       np.zeros((nruns,ngmax)) 
Uphi_sol_kms=  np.zeros((nruns,ngmax)) 
r_om_kpc=      np.zeros((nruns,ngmax)) 
R_kappa=       np.zeros((nruns,ngmax)) 

r_in=          np.zeros((nruns,ngmax)) 
r_disk_kpc=    np.zeros((nruns,ngmax)) 
r_sol_kpc=     np.zeros((nruns,ngmax)) 
r1_kpc=        np.zeros((nruns,ngmax)) 
ctau=          np.zeros((nruns,ngmax)) 
nn=            np.zeros((nruns,ngmax)) 
lam=           np.zeros((nruns,ngmax)) 
C_alp=         np.zeros((nruns,ngmax)) 
alpceil=       np.zeros((nruns,ngmax)) 
Rm_inv=        np.zeros((nruns,ngmax)) 

etat_sol_cm2s= np.zeros((nruns,ngmax)) 
td_sol_Gyr=    np.zeros((nruns,ngmax)) 
tau_sol_Gyr=   np.zeros((nruns,ngmax)) 
etat_sol=      np.zeros((nruns,ngmax)) 
td_sol=        np.zeros((nruns,ngmax)) 
tau_sol=       np.zeros((nruns,ngmax)) 

simulationoutputdirec[0]= s1+'output/'+s0[0]
if s0[-1] != s0[0]:
  simulationoutputdirec[1]= s1+'output/'+s0[1]
programdirec= s1+'program/'
print 'simulationoutputdirec[0]=', simulationoutputdirec[0]
if s0[-1] != s0[0]:
  print 'simulationoutputdirec[1]=', simulationoutputdirec[1]
print 'programdirec=', programdirec

def s2(j):
  id='_'+'%08d' %(j+1) #id number of first galaxy to compare
  return id

for i in range(nruns):
  os.chdir(simulationoutputdirec[i])
  print 'current working directory=',os.getcwd()
  for j in range(int(ngals[i])):
    paramf=open('param'+s2(j)+'.out','r')
    words=paramf.readline().split()	#reads line1 into string list
    t0_Gyr       [i,j]= float(words[ 0])
    t0_kpcskm    [i,j]= float(words[ 1])
    h0_kpc       [i,j]= float(words[ 2])
    etat0_cm2s   [i,j]= float(words[ 3])
    n0_cm3       [i,j]= float(words[ 4])
    B0_mkG       [i,j]= float(words[ 5])
    words=paramf.readline().split()	#reads line2 into string list
    nvar         [i,j]= float(words[ 0])
    dt           [i,j]= float(words[ 1])
    n1           [i,j]= int  (words[ 2])
    n2           [i,j]= int  (words[ 3])
    dx           [i,j]= float(words[ 4])
    nxphys       [i,j]= int  (words[ 5])
    nxghost      [i,j]= int  (words[ 6])
    nx           [i,j]= int  (words[ 7])
    words=paramf.readline().split()	#reads line3 into string list
    l_sol_kpc    [i,j]= float(words[ 0])
    r_l_kpc      [i,j]= float(words[ 1])
    v_sol_kms    [i,j]= float(words[ 2])
    r_v_kpc      [i,j]= float(words[ 3])
    n_sol_cm3    [i,j]= float(words[ 4])
    r_n_kpc      [i,j]= float(words[ 5])
    Uz_sol_kms   [i,j]= float(words[ 6])
    r_Uz_kpc     [i,j]= float(words[ 7])
    h_sol_kpc    [i,j]= float(words[ 8])
    r_h_kpc      [i,j]= float(words[ 9])
    Uphi_sol_kms [i,j]= float(words[10])
    r_om_kpc     [i,j]= float(words[11])
    R_kappa      [i,j]= float(words[12])
    words=paramf.readline().split()	#reads line4 into string list
    r_in         [i,j]= float(words[ 0])
    r_disk_kpc   [i,j]= float(words[ 1])
    r_sol_kpc    [i,j]= float(words[ 2])
    r1_kpc       [i,j]= float(words[ 3])
    ctau         [i,j]= float(words[ 4])
    nn           [i,j]= int  (words[ 5])
    lam          [i,j]= float(words[ 6])
    C_alp        [i,j]= float(words[ 7])
    alpceil      [i,j]= float(words[ 8])
    Rm_inv       [i,j]= float(words[ 9])
    words=paramf.readline().split()	#reads line5 into string list
    etat_sol_cm2s[i,j]= float(words[ 0])
    td_sol_Gyr   [i,j]= float(words[ 1])
    tau_sol_Gyr  [i,j]= float(words[ 2])
    etat_sol     [i,j]= float(words[ 3])
    td_sol       [i,j]= float(words[ 4])
    tau_sol      [i,j]= float(words[ 5])
  
for i in range(nruns):
  for j in range(int(ngals[i])):
    if nx[i,j] != nx[0,0]:
      print ''
      print '***WARNING: grid for run ',i,' galaxy ',s2(j),' has a different size then for run 0, galaxy _00000001'

for i in range(nruns):
  print 'RUN #',i,' UNITS'
  for j in range(int(ngals[i])):
    print '  Galaxy #',s2(j)
    print '   ','t0_Gyr=%f, t0_kpcskm=%f, h0_kpc=%f, etat0_cm2s=%f, n0_cm3=%f, B0_mkG %f' \
           %(t0_Gyr[i,j],t0_kpcskm[i,j],h0_kpc[i,j],etat0_cm2s[i,j],n0_cm3[i,j],B0_mkG[i,j])
for i in range(nruns):
  print 'RUN #',i,' NUMERICS:'
  for j in range(int(ngals[i])):
    print '  Galaxy #',s2(j)
    print '   ','nvar=%d, dt=%f, n1=%d, n2=%d, dx=%f, nxphys=%d, nxghost=%d, nx=%d' \
           %(nvar[i,j],dt[i,j],n1[i,j],n2[i,j],dx[i,j],nxphys[i,j],nxghost[i,j],nx[i,j])
print ''
for i in range(nruns):
  print 'RUN #',i,' INPUT PARAMETERS:'
  for j in range(int(ngals[i])):
    print '  Galaxy #',s2(j)
    print '   ','l_sol_kpc=%f, r_l_kpc=%f, v_sol_kms=%f, r_v_kpc=%f, n_sol_cm3=%f, r_n_kpc=%f,' \
          %(l_sol_kpc[i,j],r_l_kpc[i,j],v_sol_kms[i,j],r_v_kpc[i,j],n_sol_cm3[i,j],r_n_kpc[i,j])
    print '   ','Uz_sol_kms=%f, r_Uz_kpc=%f, h_sol_kpc=%f, r_h_kpc=%f, Uphi_sol_kms=%f, r_om_kpc=%f, R_kappa=%f' \
          %(Uz_sol_kms[i,j],r_Uz_kpc[i,j],h_sol_kpc[i,j],r_h_kpc[i,j],Uphi_sol_kms[i,j],r_om_kpc[i,j],R_kappa[i,j])
print ''
for i in range(nruns):
  print 'RUN #',i,' OTHER FREE PARAMETERS:'
  for j in range(int(ngals[i])):
    print '  Galaxy #',s2(j)
    print '   ','r_in=%f, r_disk_kpc=%f, r_sol_kpc=%f, r1_kpc=%f, ctau=%f, nn=%d, lam=%f, C_alp=%f, alpceil=%f, Rm_inv=%f' \
           %(r_in[i,j],r_disk_kpc[i,j],r_sol_kpc[i,j],r1_kpc[i,j],ctau[i,j],nn[i,j],lam[i,j],C_alp[i,j],alpceil[i,j],Rm_inv[i,j])
print ''
for i in range(nruns):
  print 'RUN #',i,' IMPORTANT CALCULATED PARAMETERS'
  for j in range(int(ngals[i])):
    print '  Galaxy #',s2(j)
    print '   ','etat_sol_cm2s=%f, td_sol_Gyr=%f, tau_sol_Gyr=%f, etat_sol=%f, td_sol=%f, tau_sol=%f' \
           %(etat_sol_cm2s[i,j],td_sol_Gyr[i,j],tau_sol_Gyr[i,j],etat_sol[i,j],td_sol[i,j],tau_sol[i,j])
print ''
  
rows=n1[0,0]+1
cols=nx[0,0]
rows=int(rows)
cols=int(cols)
print rows,cols

def readfn(file,numel):
  words=file.readline().split()
  par=np.zeros(numel)
  for j in range(len(words)):
    par[j]=float(words[j])
  return par
 
r         =np.zeros((nruns,ngmax,cols))
h_init    =np.zeros((nruns,ngmax,cols))
om_init   =np.zeros((nruns,ngmax,cols))
G_init    =np.zeros((nruns,ngmax,cols))
Uz_init   =np.zeros((nruns,ngmax,cols))
Ur_init   =np.zeros((nruns,ngmax,cols))
l_init    =np.zeros((nruns,ngmax,cols))
v_init    =np.zeros((nruns,ngmax,cols))
etat_init =np.zeros((nruns,ngmax,cols))
tau_init  =np.zeros((nruns,ngmax,cols))
alpk_init =np.zeros((nruns,ngmax,cols))
n_init    =np.zeros((nruns,ngmax,cols))
Beq_init  =np.zeros((nruns,ngmax,cols))
Br_init   =np.zeros((nruns,ngmax,cols))
Bp_init   =np.zeros((nruns,ngmax,cols))
Bzmod_init=np.zeros((nruns,ngmax,cols))

for i in range(nruns):
  os.chdir(simulationoutputdirec[i])
  print 'current working directory=',os.getcwd()
  
  for j in range(int(ngals[i])):
    isf=open('init'+s2(j)+'.out','r')
    
    r         [i,j,:]= readfn(isf,cols)
    h_init    [i,j,:]= readfn(isf,cols)
    om_init   [i,j,:]= readfn(isf,cols)
    G_init    [i,j,:]= readfn(isf,cols)
    Uz_init   [i,j,:]= readfn(isf,cols)
    Ur_init   [i,j,:]= readfn(isf,cols) 
    l_init    [i,j,:]= readfn(isf,cols)
    v_init    [i,j,:]= readfn(isf,cols)
    etat_init [i,j,:]= readfn(isf,cols)
    tau_init  [i,j,:]= readfn(isf,cols)
    alpk_init [i,j,:]= readfn(isf,cols)
    n_init    [i,j,:]= readfn(isf,cols)
    Beq_init  [i,j,:]= readfn(isf,cols)
    Br_init   [i,j,:]= readfn(isf,cols)
    Bp_init   [i,j,:]= readfn(isf,cols)
    Bzmod_init[i,j,:]= readfn(isf,cols)

for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' r=          ',r         [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' h_init=     ',h_init    [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' om_init=    ',om_init   [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' G_init=     ',G_init    [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' Uz_init=    ',Uz_init   [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' Ur_init=    ',Ur_init   [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' l_init=     ',l_init    [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' v_init=     ',v_init    [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' etat_init=  ',etat_init [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' tau_init=   ',tau_init  [i] 
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' alpk_init=  ',alpk_init [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' n_init=     ',n_init    [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' Beq_init=   ',Beq_init  [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' Br_init=    ',Br_init   [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' Bp_init=    ',Bp_init   [i]
for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' Bzmod_init= ',Bzmod_init[i]

ts_t=     np.zeros((nruns,ngmax,rows))
ts_Br=    np.zeros((nruns,ngmax,rows,cols))
ts_Bp=    np.zeros((nruns,ngmax,rows,cols))
ts_alpm=  np.zeros((nruns,ngmax,rows,cols))
ts_Bzmod= np.zeros((nruns,ngmax,rows,cols))
ts_h=     np.zeros((nruns,ngmax,rows,cols))
ts_Om=    np.zeros((nruns,ngmax,rows,cols))
ts_G=     np.zeros((nruns,ngmax,rows,cols))
ts_l=     np.zeros((nruns,ngmax,rows,cols))
ts_v=     np.zeros((nruns,ngmax,rows,cols))
ts_etat=  np.zeros((nruns,ngmax,rows,cols))
ts_tau=   np.zeros((nruns,ngmax,rows,cols))
ts_alpk=  np.zeros((nruns,ngmax,rows,cols))
ts_Uz=    np.zeros((nruns,ngmax,rows,cols))
ts_Ur=    np.zeros((nruns,ngmax,rows,cols))
ts_n=     np.zeros((nruns,ngmax,rows,cols))
ts_Beq=   np.zeros((nruns,ngmax,rows,cols))

for i in range(nruns):
  os.chdir(simulationoutputdirec[i])
  print 'current working directory=',os.getcwd()
  
  for j in range(int(ngals[i])):
    tsf=open('ts'+s2(j)+'.out','r')
  
    ts_t     [i,j,:]=  readfn(tsf,rows)
    ts_Br    [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T #read as 1-D array, convert to 2-D array, transpose
    ts_Bp    [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    if (nvar[i,j]!=2 or nvar[i,j]!=4):	#condition that dynamical quenching was used so alpha_m was calculated
      ts_alpm[i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_Bzmod [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_h     [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_Om    [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_G     [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_l     [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_v     [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_etat  [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_tau   [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_alpk  [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_Uz    [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_Ur    [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_n     [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T
    ts_Beq   [i,j,:]=  readfn(tsf,rows*cols).reshape(cols,rows).T

for i in range(nruns):
  for j in range(int(ngals[i])):
    print 'Run #',i,' Galaxy',s2(j),' t:'
    print 'ts_t[  1]=', ts_t[i,j,  0]
    print 'ts_t[421]=', ts_t[i,j,420]
    print 'ts_t[100]=', ts_t[i,j, 99]
    print ''
    print 'Run #',i,' Galaxy',s2(j),' Br(t,r):'
    print 'ts_Br [  1, 1]=',ts_Br [i,j,  0, 0]
    print 'ts_Br [  1,37]=',ts_Br [i,j,  0,36]
    print 'ts_Br [421, 1]=',ts_Br [i,j,420, 0]
    print 'ts_Br [421,37]=',ts_Br [i,j,420,36]
    print 'ts_Br [100,10]=',ts_Br [i,j, 99, 9]
    print ''
    print 'Run #',i,' Galaxy',s2(j),' Beq(t,r):'
    print 'ts_Beq[  1, 1]=',ts_Beq[i,j,  0, 0]
    print 'ts_Beq[  1,37]=',ts_Beq[i,j,  0,36]
    print 'ts_Beq[421, 1]=',ts_Beq[i,j,420, 0]
    print 'ts_Beq[421,37]=',ts_Beq[i,j,420,36]
    print 'ts_Beq[100,10]=',ts_Beq[i,j, 99, 9]

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
