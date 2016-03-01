import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator




fo=open("Files/ts_00000001.out", "r")
line=fo.readline()
Thing=line.split()
a=np.zeros(len(Thing))
for i in range(len(Thing)):
    a[i]=float(Thing[i])
print len(a)


line=fo.readline()
Thing=line.split()
b=np.zeros(len(Thing))
for i in range(len(Thing)):
    b[i]=float(Thing[i])

rows=len(a)
cols=len(b)/len(a)
bb=b.reshape(cols,rows)
print np.shape(bb)
print "0,0",bb[0,0]
print "37,0",bb[36,0]
print "421,0", bb[0,420]
print "421,37", bb[36,420]
print "100, 10", bb[9,99]



fo2=open("Files/param_00000001.out", "r")
line=fo2.readline()
Thing=line.split()
var1=float(Thing[0])
var2=float(Thing[1])

vari_l1=np.zeros(len(Thing))
for i in range(len(Thing)):
    vari_l1[i]=float(Thing[i])



#Very Quick Plot
CMap='spectral'

fig=plt.figure(figsize=(10,10))
left=0.1
base=0.1
w=0.8
h=0.8
ax=fig.add_axes([left,base,w,h])
#ax.imshow(bb,cmap=CMap)

nx, ny = (421, 37)
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
xv, yv = np.meshgrid(x, y)
print np.shape(xv)

ax.pcolormesh(xv,yv,bb,cmap=CMap)


plt.savefig("First2DPlot.png", format='png')

