#!/usr/bin/env python
# This Python file uses the following encoding: utf-8


#This script is used to draw the stator shape in matplotlib.



import numpy as np
import matplotlib.pyplot as plt

##############################################################################
##############################################################################
class Common_tangent():
    # This class is used to calculate the common tangents for two circles.
    # 这个类是用来计算两个给定圆的公切线，及切点

    def __init__(self,x1,y1,r1,x2,y2,r2):
        self.x1  = x1
        self.y1  = y1
        self.r1  = r1
        self.x2  = x2
        self.y2  = y2
        self.r2  = r2
        self.outer_tangent()

    def outer_tangent(self):
        #这个函数的功能是返回两个给定圆的公切线的直线参数 k,b.  y = kx+b
        # x1,y1 是主圆的圆心坐标 r1 是主圆的半径
        x1  = self.x1
        y1  = self.y1
        r1  = self.r1
        x2  = self.x2
        y2  = self.y2
        r2  = self.r2


        #http://bbs.emath.ac.cn/thread-5803-1-1.html
        kA = (-x2+x1+r1-r2) * (x2-x1+r1-r2)
        kB = 2*(y1-y2) * (x1-x2)
        kC = (y2+r1-y1-r2) * (-y2+r1+y1-r2)

        k1 = (-kB + np.sqrt(kB*kB - 4.*kA*kC )) / (2.*kA)
        k2 = (-kB - np.sqrt(kB*kB - 4.*kA*kC )) / (2.*kA)

        bA = (-x2+x1+r1-r2)*(x2-x1+r1-r2)
        bB = (2*r1*r1*y2 - 2*r1*r2*y1 - 2*r1*r2*y2 + 2*r2*r2*y1 - 2*x1*x1*y2+2*x1*x2*y1+2*x1*x2*y2 - 2*x2*x2*y1)
        bC = r1*r1*x2*x2+r1*r1*y2*y2 - 2*r1*r2*x1*x2 - 2*r1*r2*y1*y2 + r2*r2*x1*x1 + r2*r2*y1*y1 - x1*x1*y2*y2+ 2*x1*x2*y1*y2 - x2*x2*y1*y1
        
        b1 = (-bB + np.sqrt(bB*bB - 4.*bA*bC)) / (2.*bA)
        b2 = (-bB - np.sqrt(bB*bB - 4.*bA*bC)) / (2.*bA)
        self.outer_tangent_k1 = k1
        self.outer_tangent_k2 = k2
        self.outer_tangent_b1 = b1
        self.outer_tangent_b2 = b2
        return 



    def outer_tangent_points_on_master(self):

        #这个函数的功能是，找出两条外公切线与主圆的两个切点
        #沿着主圆心到次圆心的方向，a1b1是在右手边的切点， a2b2是在左手边的切点。所以要好好使用这个特性
        # a1,b1 
        a1 =    (self.x1 + (self.outer_tangent_b1 + self.y1)*self.outer_tangent_k1) /(1 + self.outer_tangent_k1**2)
        b1 =    self.outer_tangent_k1*a1 - self.outer_tangent_b1

        a2 =    (self.x1 + (self.outer_tangent_b2 + self.y1)*self.outer_tangent_k2) /(1 + self.outer_tangent_k2**2)
        b2 =    self.outer_tangent_k2*a2 - self.outer_tangent_b2
        return a1, b1, a2, b2

##############################################################################
##############################################################################
def cart2polar(x,y,z):
    """ 
    convertes cartesian cooridinates to polar coordinates
    """
    r = (x**2+y**2)**0.5    
    theta = np.arctan2(y,x)
    return r,theta,z
##
##############################################################################
##############################################################################

def polar2cart(r,theta,z):
    """ 
    convertes polar cooridinates to cartesian coordinates
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta) 
    return x,y,z
##############################################################################


# Copy the output list here, add comma
properties = [0.0441168, 20, 0.001, 0.638805473827, 0.231348978961, 282.942902455, 71.4127130025, 73.8007495533, 90.187769206, 268.184362585  ]


# define the basic parameters for the figure
Origin    = (0,0)
start_theta = 0. /180.*np.pi




# rotor properties/target properties
r4       = 0.042016   # m
alpha4   = 73.300756  # degree



# input value for draw a stator 
r3          = properties[0]            # stator outlet radius
r1          = r3 + 0.2*r4              # define the rotor inlet radius
Zs          = properties[1]            # stator number [-]
delta       = 2./180*np.pi             # convert to radians
t_thickness = properties[2]            # trailing edge thickness
r_trailing  = t_thickness / 2.         # radius of trailing edge circle
alpha3b     = properties[7]/180.*np.pi # blade angle, convert to radians


#Initialize the figure
fig, ax = plt.subplots( figsize = (10,10)) # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()





# plot the rotor  inlet circle
circ_rotor  = plt.Circle( Origin, r4, color = 'k', fill = False)

# plot the stator outlet circle
circ_stator = plt.Circle( Origin, r3, color = 'k', fill = False ) 

# plot the stator inlet circle
circ_inlet  = plt.Circle( Origin, r1, color = 'k', fill = False)

# plot the trailing edge circle

Tx = r3 * np.cos(start_theta)
Ty = r3 * np.sin(start_theta)
circ_trailing = plt.Circle( (Tx, Ty), t_thickness/2., color = 'r', fill =False   )



# plot the leading edge circle
# solve Ax**2 +Bx + C = 0 to get the leading edge circle centre
A    = 1 + np.tan(alpha3b)**2 
B    = 2 * r3
C    = r3**2 - r1**2
delta_x = (-B + np.sqrt( B**2 - 4 * A * C  ) ) / (2*A)

Ly   = delta_x * np.tan(alpha3b)
Lx   = Tx + delta_x
fig, ax = plt.subplots( figsize = (10,10))
chordL = np.sqrt((Tx - Lx)**2 + (Ty - Ly)**2 ) # the chord length
r_leading = t_thickness / 2. + chordL * np.tan(delta)
circ_leading = plt.Circle( (Lx, Ly), r_leading, color = 'r', fill =False )



# plot the common tangents line:

tangents_on_traling_circle = Common_tangent(Tx, Ty, r_trailing , Lx, Ly, r_leading)
T1x,T1y, T2x,T2y = tangents_on_traling_circle.outer_tangent_points_on_master()


tangents_on_leading_circle = Common_tangent(Lx, Ly, r_leading, Tx, Ty, r_trailing )
L1x,L1y, L2x,L2y = tangents_on_leading_circle.outer_tangent_points_on_master()

ax.plot([T1x,L1x],[T1y,L1y], 'r-')
ax.plot([T2x,L2x],[T2y,L2y], 'r-')





ax.add_artist(circ_rotor)
ax.add_artist(circ_stator)
ax.add_artist(circ_inlet)
ax.add_artist(circ_trailing)
ax.add_artist(circ_leading)
#ax.add_artist(line1)
ax.set_xlim(0.03,0.06)
ax.set_ylim(-0.01,0.02)
plt.tight_layout()

plt.show()


