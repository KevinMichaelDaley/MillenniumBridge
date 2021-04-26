#top=0.96,
#bottom=0.06,
#left=0.06,
#right=0.98,
#hspace=0.12,
#wspace=0.2
figX=3.5
figY=8
LW=3.25
LW3=3.75
LW2=3
pt3=pt4=15
pt5=16
pt6=17
pt7=25
pt8=18
import numpy as np, matplotlib
from scipy.interpolate import interp1d as interp
import matplotlib.pyplot as plt
r0=np.linspace(0.25,2.8,100)
omega0=[0.6,0.9,1.13,1.2]
ratio,fp=np.meshgrid(r0,omega0)
omega=2.0*np.pi*ratio*fp
Omega=np.sqrt(9.81/1.17)
matplotlib.font_manager.findSystemFonts(fontpaths='/usr/share/fonts/TTF', fontext='ttf')
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
    "font.size": 7.0,
    "pdf.fonttype" : 42
})

T=1.0/fp
a3=1+np.exp(Omega*T/2.0)
a1=1-np.exp(Omega*T/2.0)*np.cos(omega*T/2.0)
a2=np.exp(Omega*T/2.0)*np.sin(omega*T/2.0)
m=76.9
kappa0=0
kappa1=0
A=1/(1+Omega**2/omega**2)
K1=1
K2=omega/Omega
b=0.0157
X=0.006
cp=-2*m*Omega*Omega*A*A/(omega*omega*T)*(K1*a1-K2*a2)

#sigma1=2*epsilon^2*m*g/L*A^2*X/(Omega^2*T)*(K1*a1-K2*a2)
#Ncr=sqrt(c0/(epsilon^2 sigma1))

Hv=cp*omega*X
fig=plt.figure()
fig.set_size_inches((3.4,4.5))
ax=plt.subplot(2,1,1)
ratioX=np.loadtxt('../m/Fh5.csv', delimiter=',').flatten()
XX=np.array(np.loadtxt('../m/XX5.csv', delimiter=',').flatten().tolist()*64)
#YY=np.loadtxt('../m/Fj8.csv', delimiter=',').flatten()
ip=np.arange(1,65)
YY=0.65+0.6*(ip-1)/64
YY=np.array(YY.tolist()*231).reshape([231,64]).T.flatten()
#YY=np.array(np.loadtxt('../m/YY8.csv', delimiter=',').flatten().tolist()*16).reshape([16,-1]).T.flatten()
ratioY=np.nan_to_num(np.loadtxt('../m/Fh5.csv', delimiter=',')).flatten()/2.0
print(ratioY.shape,XX.shape,YY.shape) 
#jxx=np.where(np.logical_and(XX/YY>0.6, XX/YY<2.2))[0]
#XX=XX[jxx]
#YY=YY[jxx]
#ratioY=ratioY[jxx]
#ratioX=ratioX[jxx]
cx=np.where(np.logical_and(np.logical_and(np.abs(YY-0.9)<0.06,  np.abs(ratioY)>0.000), np.abs(ratioY)<250))[0]
#cx=np.where(np.logical_and(np.logical_and(np.abs(YY-0.9)<0.04,(np.abs(ratioY)>0.000)), np.abs(ratioY-)<30))[0]
hist, bin_edges=np.histogram((XX/YY)[cx], bins=40)
H0=[]
B0=[]
XXb=[]
c=[]

for i in range(len(bin_edges)-1):
    a=bin_edges[i]
    b=bin_edges[i+1]
    hh=(-ratioY/XX/2/np.pi/X)[cx][np.where(np.logical_and((XX/YY)[cx]>=a,(XX/YY)[cx]<b))]
    if len(hh)<=0:
         continue
    if 1.03/0.9>=a and 1.03/0.9<b:
         bin_edges[i]=1.03/0.9-0.05
         bin_edges[i+1]=1.03/0.9+0.05
         c.append([0.8,0,1,1])
         i1=i
    elif 1.03/1.05>=a and 1.03/1.05<b:
         bin_edges[i]=1.03/1.05-0.05
         bin_edges[i+1]=1.03/1.05+0.05
         c.append([1,0.8,0,1])
         i2=i
    else:
         c.append([0,0,0,1])         

H1=[]
for i in range(len(bin_edges)-1):
    a=bin_edges[i]
    b=bin_edges[i+1]
    hh=(-ratioY/XX/2/np.pi/X)[cx][np.where(np.logical_and((XX/YY)[cx]>=a-0.15,(XX/YY)[cx]<b+0.15))]
    if len(hh):
     H0.append(np.mean(hh))
     H1.append(np.percentile(hh,5))
     B0.append((a+b)/2)

#ax.plot(ratio[0,:].flatten(), -76.9/113000*(Hv[0,:]).flatten()/omega[0], label='$f_p=%.1f$'%fp[0,0], c='r',linewidth=6)
#ax.plot(ratio[2,:].flatten(), -76.9/113000*(Hv[2,:]).flatten()/omega[2], label='$f_p=%.1f$'%fp[2,0], c='g', linewidth=6)
ax.plot(ratio[1,:].flatten(), (cp[1]).flatten(), label='$f_p=%.1f$'%fp[1,0], c='maroon', linewidth=LW)
	#ax.plot(ratio[2,:].flatten(), -76.9/113000*(Hv[2,:]).flatten()/omega[2], label='$f_p=%.1f$'%fp[1,0], c='black', linewidth=LW)

ax.plot(B0,H0, c='black', linewidth=LW, linestyle='dotted')
ax.plot(B0,H1, c='r', linewidth=LW, linestyle='dotted')
#ax.scatter(B0,H0,50,c=c, zorder=10000)
ax.scatter(B0[i1:i1+1], H0[i1:i1+1],50,c=c[i1:i1+1], zorder=10001)
ax.scatter(B0[i2:i2+1],H0[i2:i2+1],50,c=c[i2:i2+1], zorder=10001)

cx=np.where(np.logical_and(np.logical_and(YY>0.6,YY<1.2),  np.abs(ratioY)>0.000))[0]
A=ax.scatter(XX[cx]/(YY[cx]),-ratioY[cx]/XX[cx]/2/np.pi/X, 10, c='cornflowerblue', alpha=0.5, label='model 2')
x0,x1=sorted(np.argsort(np.abs(np.array(B0)-1.03/0.9))[:2])
x2,x3=sorted(np.argsort(np.abs(np.array(B0)-1.03/1.35))[:2])
a=((1.03/0.9)-B0[x0])/(B0[x1]-B0[x0])
b=((1.03/1.35)-B0[x2])/(B0[x3]-B0[x2])
#ax.scatter([1.03/0.9], [a*H0[x0]+(1-a)*H0[x1]], 50, c=[0.8,0,1,1], zorder=30000000)
         #ax.scatter([1.03/1.35], [b*H0[x3]+(1-b)*H0[x2]], 120, marker='d', c='violet', zorder=30000000)
         #ax.plot([1.03/0.9, 1.03/0.9], [-100,100],'--', color='orange', linewidth=2.5)
         #ax.plot([1.03/1.35, 1.03/1.35],[-100,100], '-.', color='purple', linewidth=2.5)

ax.set_xlim([0.25,2.0])
#ax.set_ylim([-0.4,0.4])
#ax.set_ylim([-12,12])
	#cb=fig.colorbar(A,ax=ax, ticks=[np.min(YY[cx])+0.01, np.max(YY[cx])*0.99])
	#cb.ax.set_yticklabels(['0.6','1.19'], fontsize=35)
ax.plot(ax.get_xlim(), [0,0], 'k--', linewidth=LW2)
ax.set_ylim([-400, 400])
ax.set_yticks([-300, 0, 300])
ax.set_yticklabels(['-300', '0', '300'], fontsize=pt6, usetex=False, fontweight='bold')

#ax.legend(fontsize=19, loc='upper right')
plt.setp(ax.get_xticklabels(), fontsize=pt7)
plt.setp(ax.get_yticklabels(), fontsize=pt7)
#plt.xlabel('$f_b/f_p$', fontsize=pt5, labelpad=10)
plt.ylabel(r'$ {\overline{\sigma}_1}$',fontsize=pt7, labelpad=8)

plt.gca().text(0,1.05,r' N $\cdot$ s/m', usetex=False, fontsize=pt6, transform=ax.transAxes)
ax2=plt.subplot(2,1,2)
#P2=np.loadtxt("../violin_macdonald_0.9.txt")
P=np.loadtxt("../../src2/violin99.txt")
Y=P[:,0]
X=P[:,1]/P[:,2]
#jx=np.where(np.logical_and(P[:,3]>=0.7, P[:,3]<=1.1))[0]
#X=X[jx]
#Y=Y[jx]
#X2=P2[:,1]/P2[:,3]
#X3=P3[:,1]/P3[:,2]
#X4=P4[:,1]/P4[:,2]
#Y2=P2[:,0]
#Y3=P3[:,0]
#Y4=P4[:,0]
#jx=np.argsort(X2)
#X2=X2[jx]
#Y2=Y2[jx]
#jx=np.argsort(X3)
#X3=X3[jx]
#Y3=Y3[jx]
#jx=np.argsort(X4)
#X4=X4[jx]
#Y4=Y4[jx]
hist, bin_edges=np.histogram(X, bins=40)
H=[]
B=[]
for i in range(len(bin_edges)-1):
    a=bin_edges[i]
    b=bin_edges[i+1]
    hh=Y[np.where(np.logical_and(X>=a,X<b))]
    if len(hh):
     H.append(hh)
     B.append((a+b)/2)
#parts1=ax2.violinplot(H[:], B[:], widths=(bin_edges[1]-bin_edges[0])*0.8,showmeans=True)
ax2.scatter(X,Y,8,color='cornflowerblue')
#ax2.scatter(X, Y, 130,'cornflowerblue', label='empirical (adaptive)')
#ax2.plot(X2,Y2,'k--', linewidth=7, label='$N_2$')
#ax2.plot(X3,Y3,'c--', linewidth=7)
#ax2.plot(X4,Y4,'m--', linewidth=7)
eps=np.sqrt(76.9/113e3)
N=np.nan_to_num(0.04*omega[1]*113e3/-cp[1])*(-cp[1,:]>0)+10000*(-cp[1,:]<0)
#N2=np.sqrt(0.04*omega[2]/(eps**4*Hv[2]/omega[2]))*(cp[2,:]>0)+10000*(cp[2,:]<0)
N3=np.nan_to_num(0.04*0.9*113e3*np.array(B0)*2*np.pi/(-np.array(H0)))*(-np.array(H0)>0)+10000*(-np.array(H0)<0)
print(N)
ax2.plot(omega[1]/(2*np.pi*0.9),N,'maroon', label='$N_1$', linewidth=LW, c='maroon')
#ax2.plot(ratio[2,:].flatten(),N2,'black', label='$N_2$', linewidth=LW, c='black')

ax2.plot(B0,N3, c='black', linewidth=LW, linestyle='dotted')
#ax2.scatter(B0,N3,50,color='black')

#ax2.plot(ax.get_xlim(), [166, 166], 'r--', linewidth=9)
#ax2.legend(fontsize=19, loc='upper right')
ax2.set_xlim(ax.get_xlim())
ax2.set_ylim([0,600])
ax2.set_xticks([0.5,1.0,1.5,2.0])
#ax.set_xticks(ax.get_xticks()[::2])
ax.set_xticks([])
ax2.set_yticks([0,200,400,600])
plt.setp(ax2.get_xticklabels(), fontsize=pt5, usetex=False, fontweight='bold')
plt.setp(ax.get_xticklabels(), fontsize=pt5, usetex=False, fontweight='bold')
plt.setp(ax2.get_yticklabels(), fontsize=pt5, usetex=False, fontweight='bold')
plt.setp(ax.get_yticklabels(), fontsize=pt5, usetex=False, fontweight='bold')
#matplotlib.mathtext.SHRINK_FACTOR = 0.7/2.0
#matplotlib.mathtext.GROW_FACTOR = 2.0 / 0.7
ax.set_ylabel(r'$\overline{\sigma}_1$', fontsize=pt7, labelpad=8)
plt.ylabel('$N_{crit}$',fontsize=pt7, labelpad=8)
plt.xlabel(r'$\left[\Omega/\bar{\omega}\right]$',fontsize=pt7, labelpad=8)
plt.tight_layout(pad=0.9)
plt.savefig('mcrobie.pdf')
plt.show()
