nt=200
ntfrac=0.8
print 'nt=%d' %nt
t1_Gyr=0.8
t2_Gyr=1.2
t1=t1_Gyr/t0_Gyr
t2=t2_Gyr/t0_Gyr
print t1_Gyr, t2_Gyr

plotdirec= s1+'plots/'+ s0
print plotdirec

r_kpc=r*r_disk_kpc
print r_kpc

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
