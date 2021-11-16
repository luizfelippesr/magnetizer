import numpy as np
import matplotlib.pyplot as plt

data  = np.loadtxt('IQU_conv.txt')
steps = data[:,3]
i = data[:,10]
q = data[:,11]
u = data[:,12]

steps =  steps*steps/10**4
plt.xscale('log')
plt.yscale('symlog')
plt.xlabel('steps/time')
plt.ylabel('I/I$_{true}$', rotation='horizontal')

plt.annotate('I ~ 1.e-7', xy=(1, -200), size=10 )
plt.annotate('Q ~ 1.e-8', xy=(1, -600), size=10 )
plt.annotate('U ~ 1.e-12', xy=(1, -2000), size=10 )

plt.plot(steps, i/i[-1], label='I')
plt.plot(steps, q/q[-1], label='Q')
plt.plot(steps, u/u[-1], label='U')
plt.legend(loc="lower right",  prop={'size': 15})

plt.savefig('convergence.png')
plt.show()
