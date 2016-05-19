import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif',size = 24)
#read data

radial=np.loadtxt("hist.data")
separation=np.linspace(0,18.1/2.0,len(radial))
#dsep = separation[1]-separation[0]
#separation += dsep/0.5
#ra`dial = radial/dsep
#rho = 1.0
#for i in range(len(radial)):
    #r = dsep * (i+0.5)
    #r = dsep
    #vol = ((i+1)**3 - (i**3))*dsep**3
    #numpart = (4.0/3.0)*np.pi*vol*rho
    #radial[i] /= (numpart)
#radial/=4*np.pi*separation
radial/=np.sum(radial[0.5*len(radial):len(radial)])/((1-0.5)*len(radial))
#radial/=radial[len(radial)-1]
print "max", max(radial), separation[np.argmax(radial)]
print "all"
#for i in range(len(radial)):
#    print radial[i], separation[i]
#radial/=np.sum(radial[len(radial)/2:len(radial)])/(0.5*len(radial))
#xminimum= min(separation)-0.05
#xmaximum = max(separation)+0.05
#yminimum= min(waterDimer)-0.005
#ymaximum = max(waterDimer)+0.005
#plot commands
fig = plt.figure(facecolor='white',figsize=(8,8))
ax = fig.add_subplot(111)
plt.xlim(0, max(separation))
plt.ylim(0, np.max(radial)+0.5)
print separation[np.argmax(radial)]
#plt.plot(radial,color='red', bins=50)
plt.plot(separation,radial,color='red')
plt.plot([separation[np.argmax(radial)], separation[np.argmax(radial)]], [0, 10], color='slategray', linestyle='--', linewidth=2)
for axis in ['top','bottom','left','right']:
	ax.spines[axis].set_linewidth(3)
#plt.xticks([0,np.pi,2*np.pi],['0',r'$\pi$',r'$2\pi$'])
ax.xaxis.set_tick_params(width=2.5,length=7)
ax.yaxis.set_tick_params(width=2.5,length=7)
plt.xlabel(r'R($\AA$)')
plt.ylabel(r'g(r)')
plt.title(r'Argon Radial Distribution function',y=1.08)
fig.tight_layout()
fig.savefig('prettyplot.png')
#plt.show()
