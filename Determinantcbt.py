import matplotlib.pyplot as plt
import numpy as np
from numpy import e,pi,cos,sin,log10,log,arctan,sqrt,arcsin,arccos,tan
from mpl_toolkits.mplot3d import Axes3D

def Det(b,c,p):
	return -49152*pow(b,4)*(-1+pow(b,2)+pow(c,2))*pow(c,4)*pow((-3*pow(c,2)+pow(\
b,2)*(-3+64*p*pow(c,2))),-1)

# values for d in the case of having ...+sqrt(...) or ...-sqrt(...)
# to be noted that when pow acts as a sqrt on negative numbers it returns np.nan
def isdRealPlus(b,c,t):
	return pow(2,-0.5)*pow((1+-1*pow(b,2)+-1*pow(c,2)+pow((-3*pow(c,2)+pow(b,2)*(\
-3+64*p*pow(c,2))),-1)*pow((-1+pow(b,2)+pow(c,2))*(-3*pow(c,2)+pow(b,\
2)*(-3+64*p*pow(c,2)))*(3*pow(c,2)+64*p*pow(b,2)*pow(c,2)*(-1+pow(b,2)\
+pow(c,2))+-3*(pow(b,4)+-1*pow(b,2)*(1+2*pow(c,2))+pow(c,4))),0.5)),0.5)

def isdRealMinus(b,c,t):
	return pow(2,-0.5)*pow((1+-1*pow(b,2)+-1*pow(c,2)+-1*pow((-3*pow(c,2)+pow(b,\
2)*(-3+64*p*pow(c,2))),-1)*pow((-1+pow(b,2)+pow(c,2))*(-3*pow(c,2)+\
pow(b,2)*(-3+64*p*pow(c,2)))*(3*pow(c,2)+64*p*pow(b,2)*pow(c,2)*(-1+\
pow(b,2)+pow(c,2))+-3*(pow(b,4)+-1*pow(b,2)*(1+2*pow(c,2))+pow(c,4))),\
0.5)),0.5)


N = 100#1000 # resolution of the grid
b = np.linspace(0,1,N)
c = np.linspace(0,1,N)
b,c = np.meshgrid(b, c)
# value of Tr[Q^{-1}], you may want to see how it changes changing p from 0.75 up to infinity
p = 1 # <--------------------------------- main parameter to be changed

dp = isdRealPlus(b,c,p)
dm = isdRealMinus(b,c,p)
Z = Det(b,c,p)
# a main concern that I missed is that for b=0 or c=0, especially with low N you
# can see how there are points that gives DetQ = 0
mask = (b**2+c**2<=1) & (~np.isnan(dp)) & (~np.isnan(dm))# & (Z>0) # recommended putting the Z>0 to see better
Z[~mask] = np.nan

# plotting
fig = plt.figure(figsize=(10, 5))
gs = fig.add_gridspec(1, 2)
# 3D Surface plot
ax1 = fig.add_subplot(gs[0], projection='3d')
sur = ax1.plot_surface(b,c, Z, cmap='autumn')

ax1.view_init(elev=49, azim=-42)
ax1.set_xlabel('b',fontsize=16, labelpad=15)
ax1.set_ylabel('c',fontsize=16, labelpad=15)
ax1.set_zlabel(r"Tr[$Q^{-1}$]",fontsize=16, labelpad=15)


ax2 = fig.add_subplot(gs[1])
c = ax2.pcolormesh(b,c,Z, shading='auto', cmap='autumn')
fig.colorbar(c, ax=ax2)
ax2.set_aspect('equal')
ax2.set_xlabel('b',fontsize=16)
ax2.set_ylabel('c',fontsize=16)


ax1.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
plt.tight_layout()
plt.subplots_adjust(wspace=0.5)

plt.show()

