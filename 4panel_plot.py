import numpy as np, matplotlib, random
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
figX=2.3622
figY=4.5
LW=1.25
LW3=1.75
LW2=1
pt3=pt4=5
pt5=6
pt6=7
pt7=7
pt8=8
#matplotlib.use('Qt4Agg')

plt.rc('text', usetex=True)
matplotlib.font_manager.findSystemFonts(fontpaths='/usr/share/fonts/truetype/msttcorefonts', fontext='ttf')
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": "Arial",
    "font.size": 7.0,
    "pdf.fonttype" : 42
})
Nmax=220
Nmin=40
def sgn(x):
	return 1-2.0*(x<0)
def rgb(*args):
        return np.array((*args,255))/255.0
threshc=rgb(255,157,179)
cmb=rgb(91, 77, 183)
thresh2c=rgb(130,173,255)
VSPAN_ALPHA=0.4
cmg=rgb(129,209,82)
cmr=rgb(255,83,71)
cmk=rgb(33,33,32)
LW4=2.5
tsc=[cmg, cmr,  thresh2c, cmb, threshc]*4
def label_all(ax, thresh, thresh2, xlabel, ylabel):
    label1=ylabel.split('$')[0].strip()
    label2=ylabel.split('$')[1].strip()
    plt.rc('text',usetex=False)
    if 'order' in label1:
            ax.set_ylabel(  # position text relative to Axes
                        label1,
                         fontsize=pt7, fontname='sans-serif', rotation=90, labelpad=5
                )
    if xlabel!='':
        ax.set_xlabel( xlabel, fontsize=pt7, fontname='sans-serif', labelpad=5)
#    plt.rc('text',usetex=True)
    ax.tick_params(axis='x', which='major', labelsize=pt6)
    ax.tick_params(axis='y', which='major', labelsize=pt6)
    ax.tick_params(axis='both', which='minor', labelsize=pt6)
    ax.set_xticks([50,100,150,200])
#    ax.axvspan(2,trans,edgecolor='white', alpha=0.4, zorder=-100, color='blue', hatch='\\', linewidth=LW)
#    ax.axvspan(trans,200,edgecolor='white', alpha=0.4, zorder=-100, color='red', hatch='/', linewidth=LW)
 
def panel4(X, B, thresh, thresh2, out, ts1, ts2, ts3, N1, N2, N3, belykh=False, adapt=True, sigma1=0, sigma2=0, sigma3=0):
    AMP=0.02
    N=X[:,0]
    C=X[:,1]
    fp=X[:,5]

    fb=[]
    for i in range(2,300,8):
                zx=np.where(np.logical_and(np.diff(B[:,2]>0)>0,np.logical_and(B[1:,0]>=i,B[1:,0]<i+8)))[0]
                ts=B[zx,1]        
                period=np.mean(np.diff(ts))
                fb.append(1.0/period)
                print(fb[-1])
    if belykh:
        for i in range(2,min(300,int(np.max(N))),1):
                zx=np.where(np.floor(B[:,0])==i)[0]
                ts=B[zx,1]        
                ys=B[zx,4:]
                tlast=np.zeros(i)
                X[i-2,4]=0
                phases=[]
                for j,tj in enumerate(ts):
                            for k in range(i):
                                if j>=1:
                                     if ys[j-1,k]<0 and ys[j,k]>=0:
                                        tlast[k]=tj
                                        phases.append(2*np.pi*fp[i-2]*(tlast[k]-tlast[0]))
                pc=np.cos(np.array(phases))
                ps=np.sin(np.array(phases))
                X[i-2,4]=np.sqrt(np.sum(pc)**2+np.sum(ps)**2)/len(phases)
#                print(X[i-2,4])
                if X[i-2,4]>0.4 and X[i-3,4]<0.4 and i>3:
                        thresh2=[i]
    #fb=X[:,-4]    
    A=X[:,2]/10
    R1=X[:,3]
    R2=X[:,4]
    T=X[:,4]
    F=plt.figure(figsize=(figX,figY))
    gs=F.add_gridspec(nrows=3, ncols=1, width_ratios=[0.9], height_ratios=[0.2,0.2, 0.15])

#    thresh=209
#    thresh=156
    ix=np.where(np.logical_and(np.isfinite(R1),np.isfinite(R2)))[0]
    ax1=F.add_subplot(gs[0,0])
    ax3 = plt.axes([0,0,1,1],label='ts1')
# Manually set the position and relative size of the inset axes within ax1
    ax1.set_xlim([Nmin,Nmax])
    points= np.array([N[ix],A[ix]]).transpose().reshape((-1,1,2))
    if not belykh:
        ax1.set_ylim(-0.001, 0.025)
        mul=1
    else:
        ax1.set_ylim(-0.001, 0.025)
        mul=1
    ip1 = InsetPosition(ax1, [(77-Nmin)/(Nmax-Nmin),(0.11),0.3,0.25])
    kmax=150
    ax3.set_axes_locator(ip1)
    ax1.plot([100,100],[A[100-2],0.03*0.025/0.1], '--', color=cmk, linewidth=LW2)
    ax1.plot([128,155],[0.06*0.025/0.1,A[155-2]], '--', color=cmk, linewidth=LW2)
    ax1.plot([180,200],[0.06*0.025/0.1,A[200-2]], '--', color=cmk, linewidth=LW2)
    for ki, k in enumerate(np.random.permutation(N1)[:20]):
            ax3.plot(ts1[:kmax,1], (ts1-0.63*sgn(ts1))[:kmax, k+4]/10, alpha=0.5, linewidth=LW2, c=tsc[ki])
    ax3.plot(ts1[:kmax,1], ts1[:kmax,2]/10, color=cmk,linewidth=LW4)
    ax3.set_title('N=%i'%N1, fontsize=pt5)
    ax4 = plt.axes([0,0,1,1], label='ts2')
    ip2 = InsetPosition(ax1, [(77-Nmin)/(Nmax-Nmin),(0.1+0.32+0.16),0.3,0.25])
    ax4.set_axes_locator(ip2)
    for ki,k in enumerate(np.random.permutation(N2)[:20]):
            ax4.plot(ts2[:kmax,1], (ts2-0.63*sgn(ts2))[:kmax, 4+k]/10, alpha=0.5, linewidth=LW2, c=tsc[ki])
    ax4.plot(ts2[:kmax,1], ts2[:kmax,2]/10, color=cmk,linewidth=LW4)
    ax4.set_ylim([-AMP*2,AMP*2])
    ax4.set_title('N=%i'%N2, fontsize=pt5)
    ax5 = plt.axes([0,0,1,1], label='ts3')
    ip3 = InsetPosition(ax1, [(162-Nmin)/(Nmax-Nmin),(0.1+0.32+0.16),0.3,0.25])
    ax5.set_axes_locator(ip3)
    for ki, k in enumerate(np.random.permutation(N3)[:20]):
        ax5.plot(ts3[:kmax,1], (ts3-0.63*sgn(ts3))[:kmax, k+4]/10.0, alpha=0.5, linewidth=LW2, c=tsc[ki])
    ax5.plot(ts3[:kmax,1], ts3[:kmax,2]/10.0, color=cmk,linewidth=LW4)
    ax5.set_title('N=%i'%N3, fontsize=pt5)
    ax5.set_ylim([-0.08,0.08])
    ax3.set_xlim([ts1[0,1],ts1[:kmax,1][-1]])
    ax3.set_xticks([])
    ax3.set_ylim([-0.08,0.08])
    ax4.set_xlim([ts2[0,1],ts2[:kmax,1][-1]])
    ax4.set_xticks([])
    ax5.set_xlim([ts3[0,1],ts3[:kmax,1][-1]])
    ax4.set_ylim([-0.08,0.08]) 
    ax4.set_yticks([-0.06,0,0.06])
    ax4.set_yticklabels(['-0.06','0','0.06'], fontsize=pt3)
    ax4.tick_params(axis='both', which='major', pad=1)
    ax3.set_yticks([])
    ax5.set_yticks([])
    ax5.set_xticks([])
    jx=[0]+(np.where(C<C[0]-0.001)[0])
    ax1.plot(N[ix],A[ix], linewidth=LW, zorder=200, color=cmb)
    ax100=F.add_subplot(gs[1,0])
    
    ax2=ax100.twinx()
    ax2.set_zorder(100)
    label_all(ax1,thresh,thresh2,'', 'bridge amplitude $A_x$')
    rline=ax2.plot(N[ix],R2[ix],  linewidth=LW, color=cmb, zorder=201, label='order parameter\n (CoM)')
    ax2.set_xlim([Nmin,Nmax])

    ax100.set_ylim([-10,5])
    ax100.set_xlim([Nmin,Nmax])
    y0=ax100.get_ylim()
    yn=(0-y0[0])/(y0[1]-y0[0])
    ax2.set_ylim(-0.125,1.125/yn-0.125)    
    print(1.125/yn-0.125)

    ax2.set_yticks([0,0.25,0.5,0.75,1.0])
    ax2.set_yticklabels(['0','','0.5','','1'])
    label_all(ax2, thresh, thresh2, '', 'order parameter $r$')
    cline=ax100.plot(N, C, linewidth=LW, color=cmr, label="effective damping", zorder=300)
    ax100.set_ylim([-10,5])
    ax100.set_xlim([Nmin,Nmax])
    label_all(ax100,thresh,thresh2,'', 'effective \n damping $c_eff$')
    ax100.plot([2,Nmax], [0,0], '--', linewidth=LW2, color=cmk, alpha=0.7)
    if len(thresh2)>0:
        ax1.axvspan(thresh, thresh2[0], alpha=VSPAN_ALPHA, color=threshc,linewidth=0)
        ax100.axvspan(thresh, thresh2[0], alpha=VSPAN_ALPHA, color=threshc,linewidth=0)
        ax1.axvspan(thresh2[0], Nmax, alpha=VSPAN_ALPHA, color=thresh2c,linewidth=0)
        ax100.axvspan(thresh2[0], Nmax, alpha=VSPAN_ALPHA, color=thresh2c,linewidth=0)
    if False:
        ip9 = InsetPosition(ax100, [0.65,0.3,0.3,0.65])
        ax9 = plt.axes([0,0,1,1], label='ts4')
        ax9.set_axes_locator(ip9)
        ax9.plot(N, C, linewidth='7',zorder=200)
        ax9.set_xlim(np.min(N), np.max(N))
        ax9.set_ylim(np.min(C), np.max(C))
    ax7=plt.axes([0,0,1,1])
    ax7.plot(N,sigma1,label='$\sigma_1$',linewidth=LW,color=cmr,zorder=200)
    ax7.plot(N,sigma2,label='$\sigma_2$',linewidth=LW, color=cmb,zorder=200)
    ax7.plot(N,sigma3,label='$\sigma_3$', linewidth=LW, color=cmg,zorder=105)
    axl1=ax7.legend(fontsize=pt4,facecolor='white', framealpha=1)
    axl1.set_zorder(500)
    ax7.set_xlim([Nmin,Nmax])
    ax7.set_ylim([-5*76.9/113e3,5*76.9/113e3])
    ax6=F.add_subplot(gs[2,0])

    ax6.fill_between(N,(X[:,5]+X[:,6]), (X[:,5]-X[:,6]),color='violet', zorder=100)
    ax6.plot(N,X[:,5],'--', color='purple', label='mean foot placement\n frequency', linewidth=LW, zorder=200)
    #ax6.plot(N,(X[:,5]+X[:,6]),'b')
    #ax6.plot(N,(X[:,5]-X[:,6]),'b')
    ax6.plot(np.arange(2,300,8),np.array(fb)/2.0, color=cmk, linewidth=LW, label='bridge frequency', zorder=202)
    ax6.set_xlim([Nmin,Nmax])
    axl=ax6.legend(fontsize=pt3, facecolor='white', framealpha=1, loc='upper left')
    axl.set_zorder(500)
    axl2=ax2.legend(handles=[cline[0],rline[0]], fontsize=pt3, loc='upper left', bbox_to_anchor=[0.01,0.55], facecolor='white', framealpha=1)
    axl2.set_zorder(600)
    ax100.yaxis.offsetText.set_text('$N \cdot s/m$')
    ax100.yaxis.offsetText.set_fontsize(pt5)
    ax1.yaxis.offsetText.set_text('$m$')
    ax1.yaxis.offsetText.set_fontsize(pt5)
    ax100.yaxis.offsetText.set_visible(True)
    ax1.yaxis.offsetText.set_visible(True)
    ax6.set_ylim([0.5,2.0])
    label_all(ax7, thresh, thresh2, '', 'damping \n terms $\sigma_{1,2,3}$')
    label_all(ax6, thresh, thresh2,'number of pedestrians', 'frequency  $f_{b},f_{p}$')
#    plt.tight_layout()
    ax6.axvspan(thresh, thresh2[0], alpha=VSPAN_ALPHA, color=threshc,linewidth=0)
    ax7.axvspan(thresh, thresh2[0], alpha=VSPAN_ALPHA, color=threshc,linewidth=0)
    ax6.axvspan(thresh2[0], Nmax, alpha=VSPAN_ALPHA, color=thresh2c,linewidth=0)
    ax7.axvspan(thresh2[0], Nmax, alpha=VSPAN_ALPHA, color=thresh2c, linewidth=0)
    ax7.set_ylim([-300, 300])
    from matplotlib.ticker import ScalarFormatter
    xfmt = ScalarFormatter()
    xfmt.set_powerlimits((-4,4))
    ax7.set_yticks([-0.002,0.0,0.002,0.004,0.006])
    ax7.yaxis.set_major_formatter(xfmt)
    ax7.yaxis.get_offset_text().set_fontsize(pt5)
    ax7.yaxis.get_offset_text().set_visible(False)
    ax7.text(0,1.15,r'$\times 10^{-3}$ N $\cdot$ s/m',transform=ax7.transAxes, fontsize=pt5)
    ax7.set_visible(False)
    ax6.text(0,1.05,'1/s',transform=ax6.transAxes, fontsize=pt5)
    ax1.text(0,1.05,'m', transform=ax1.transAxes, fontsize=pt5)
    ax100.text(0,1.05,'1/s',transform=ax100.transAxes, fontsize=pt5)

    ax1.text(thresh,0.04*0.025/0.1*1.1, 'unstable (no sync)', fontsize=pt4, rotation=90, usetex=False, fontname='sans-serif')

    plt.tight_layout(pad=0.1,w_pad=0.1,h_pad=0.1)

    F.savefig(out+'.svg')

    bb1=ax1.yaxis.get_label().get_window_extent()
    bb2=ax2.yaxis.get_label().get_window_extent()
    x0=(bb1.x0+bb1.width)/F.get_dpi()
    x1=figX
    ax100.yaxis.get_label().set_visible(False)
    ax6.yaxis.get_label().set_visible(False)
        
    F.savefig(out+'_cropped.pdf', bbox_inches=matplotlib.transforms.Bbox([[x0,0],[x1,figY]]))
    return F
import sys
A=np.loadtxt('stdout.txt')
B=np.loadtxt('stderr.txt')
sigma1=A[:,-3]
sigma2=A[:,-2]
N=A[:,0]
A[:,1]=0.04*2*np.pi*1.03-np.sqrt(76.9/113e3)*(sigma1+sigma2)
tr100=B[np.where(B[:,0]==100)[0]]
tr155=B[np.where(B[:,0]==155)[0]]
tr200=B[np.where(B[:,0]==200)[0]]
panel4(A,B, 146, [165], 'our_final3',tr100, tr155,tr200, 100,155,200,adapt=False, belykh=True,sigma1=-A[:,-3]*76.9/N,sigma2=-A[:,-2]*76.9/N, sigma3=-A[:,-1]*76.9/N)
plt.show()
