#!/usr/bin/env python
# This Python file uses the following encoding: utf-8


#This script is used to draw the stator shape in matplotlib.



import numpy as np
import matplotlib.pyplot as plt



##set plot format

plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['font.serif'] = 'Times New Roman'
#plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 20
plt.rcParams['axes.labelsize'] = 22
#plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['figure.titlesize'] = 18
textXLoc = 193
#fig = plt.figure(figsize=(16,18))







# Control the figure plot:
showQuater = 1  # set 1 to show the quater, otherwise show whole

#xlimit = (0, 0.04)
#ylimit = (0, 0.04)
lineWidth = 2.

# after which blade, you want to plot the velocity vector?
stator_flow_vec_pos = 0
rotor_flow_vec_pos  = 1 

# Copy the output list here, add comma
#1.05 12
properties = [0.03321045, 12, 0.0009, 0.67609781745, 0.184640264705, 298.277636921, 65.4370729862, 73.1480623984,123.99174314, 271.285083114 ,15089202.4999, 793.962816949,98.5678976607]
#1.05 16
properties = [0.03321045, 16, 0.0009, 0.677693783761, 0.182526427468, 298.930153716, 65.164848048, 67.6583831808, 125.553337193, 271.285083114, 15070006.7417, 793.789439692, 98.4692565668]
# 1.05 22
properties = [0.03321045, 22, 0.0009, 0.681305836338 ,0.178016269223 ,300.405811272, 64.5635234796, 61.5049932346, 129.027342551, 271.285083114 ,15026512.3731 ,793.395898482 ,98.245626317 ]
# 1.05 16 with loss
properties = [0.03321045, 16, 0.0009, 0.646835985316 ,0.182560673609 ,299.0196043 ,65.1278395993, 67.615724891, 125.766161727, 271.285083114 ,20808743.7497, 838.794853038, 126.583245625 ]
# define the basic parameters for the figure
Origin    = (0,0)
start_theta = 0. /180.*np.pi


# rotor properties/target properties
r4       =  0.031629   # m
r6h      =  0.009489
r6t      =  0.018531
Zr       =  13
alpha4   =  66.16126 / 180.*np.pi # rotor inlet flow angle.
C4       =  311.417414            # define the rotor inlet velocity here




# input value for draw a stator 
r3          = properties[0]            # stator outlet radius
r1          = r3 + 0.2*r4              # define the rotor inlet radius
Zs          = properties[1]            # stator number [-]
delta       = 2./180*np.pi             # convert to radians
t_thickness = properties[2]            # trailing edge thickness
r_trailing  = t_thickness / 2.         # radius of trailing edge circle
C3          = properties[5]            # stator outlet flow velocity.
alpha3      = properties[6]/180.*np.pi # stator outlet flow angle.
alpha3b     = properties[7]/180.*np.pi # blade angle, convert to radians






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
def cart2polar(x,y):
    """ 
    convertes cartesian cooridinates to polar coordinates
    """
    r = (x**2+y**2)**0.5    
    theta = np.arctan2(y,x)
    return r,theta
##
##############################################################################
##############################################################################

def polar2cart(r,theta):
    """ 
    convertes polar cooridinates to cartesian coordinates
    """
    x = r * np.cos(theta)
    y = r * np.sin(theta) 
    return x,y
##
##############################################################################
##############################################################################
def rotatePoints(theta, X_start, Y_start):
    """
    function to calculate the rotated points coordinates, with the origin as centre
    """
    R_start , T_start = cart2polar(X_start, Y_start)
    X_end , Y_end = polar2cart( R_start, T_start+theta)
    return X_end, Y_end
##############################################################################



#Initialize the figure
##################################################
fig, ax = plt.subplots( figsize = (8,8)) # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()
##################################################


# plot the rotor outlet circle
circ_rotor_out_tip  = plt.Circle( Origin, r6t, color = 'k', fill = False, lw = 2.)
circ_rotor_out_hub  = plt.Circle( Origin, r6h, color = 'k', fill = False, lw = 2.)

# Initialize the theta angle list, for proper copy the blades drawing.
theta_list = np.arange( 0.,  2*np.pi + 2*np.pi/Zs,  2*np.pi/Zs )
rotor_theta_list =  np.arange( 0.,  2*np.pi + 2*np.pi/Zr,  2*np.pi/Zr )
###
##
# plot the rotor blades at outlet
Bhx, Bhy = r6h, 0.
Btx, Bty = r6t, 0.

rotor_out_hub_list = []
rotor_out_tip_list = []
for i in range(len(rotor_theta_list)):
    rotor_out_hub_list.append( rotatePoints(rotor_theta_list[i], Bhx, Bhy) )
    rotor_out_tip_list.append( rotatePoints(rotor_theta_list[i], Btx, Bty) )

for j in range(len(rotor_theta_list)):
    ax.plot([ rotor_out_hub_list[j][0], rotor_out_tip_list[j][0] ],[ rotor_out_hub_list[j][1], rotor_out_tip_list[j][1]],'k-', lw =1.5)
##
###
# plot the rotor  inlet circle
circ_rotor  = plt.Circle( Origin, r4, color = 'k', fill = False , lw = lineWidth)
##
###
# plot the stator outlet circle
circ_stator = plt.Circle( Origin, r3, color = 'k', fill = False, lw = lineWidth ) 
##
###
# plot the stator inlet circle
circ_inlet  = plt.Circle( Origin, r1, color = 'k', fill = False, lw = lineWidth)

##
###
# plot the trailing edge circle
Tx = r3 * np.cos(start_theta)
Ty = r3 * np.sin(start_theta)


trailing_circle_centre_X = []
trailing_circle_centre_Y = []
trailing_circle_list     = []
for i in range(len(theta_list)):
    tempX, tempY = rotatePoints(theta_list[i], Tx, Ty)
    trailing_circle_centre_X.append(tempX)
    trailing_circle_centre_Y.append(tempY)


for j in range(len(theta_list)):
    trailing_circle_list.append( plt.Circle( (trailing_circle_centre_X[j], trailing_circle_centre_Y[j]),
                                 r_trailing, color = 'r', fill =False, lw = lineWidth)   )
##
###
# plot the leading edge circle
# solve Ax**2 +Bx + C = 0 to get the leading edge circle centre
A    = 1 + np.tan(alpha3b)**2 
B    = 2 * r3
C    = r3**2 - r1**2
delta_x = (-B + np.sqrt( B**2 - 4 * A * C  ) ) / (2*A)

Ly   = -delta_x * np.tan(alpha3b)
Lx   = Tx + delta_x

chordL = np.sqrt((Tx - Lx)**2 + (Ty - Ly)**2 ) # the chord length
r_leading = t_thickness / 2. + chordL * np.tan(delta)

leading_circle_centre_X = []
leading_circle_centre_Y = []
leading_circle_list     = []
for i in range(len(theta_list)):
    tempX, tempY = rotatePoints(theta_list[i], Lx, Ly)
    leading_circle_centre_X.append(tempX)
    leading_circle_centre_Y.append(tempY)
for j in range(len(theta_list)):
    leading_circle_list.append( plt.Circle( (leading_circle_centre_X[j], leading_circle_centre_Y[j]),
                                 r_leading, color = 'r', fill =False, lw = lineWidth)   )

#circ_leading = plt.Circle( (Lx, Ly), r_leading, color = 'r', fill =False )

##
###
# plot the common tangents line:
tangents_on_traling_circle = Common_tangent(Tx, Ty, r_trailing , Lx, Ly, r_leading)
T1x,T1y, T2x,T2y = tangents_on_traling_circle.outer_tangent_points_on_master()

tangents_on_leading_circle = Common_tangent(Lx, Ly, r_leading, Tx, Ty, r_trailing )
L1x,L1y, L2x,L2y = tangents_on_leading_circle.outer_tangent_points_on_master()

T1_list = []
T2_list = []
L1_list = []
L2_list = []

for i in range(len(theta_list)):
    T1_list.append(  rotatePoints(theta_list[i], T1x, T1y) )
    T2_list.append(  rotatePoints(theta_list[i], T2x, T2y) )
    L1_list.append(  rotatePoints(theta_list[i], L1x, L1y) )
    L2_list.append(  rotatePoints(theta_list[i], L2x, L2y) )


for j in range(len(theta_list)):
    ax.plot([ T1_list[j][0], L1_list[j][0] ],[ T1_list[j][1], L1_list[j][1]],'r-', lw = lineWidth)
    ax.plot([ T2_list[j][0], L2_list[j][0] ],[ T2_list[j][1], L2_list[j][1]],'r-', lw = lineWidth)

##
###
# plot the stator velocity vec arrow
# The centre for the trailing edge
Tx_vec =     trailing_circle_centre_X[stator_flow_vec_pos]
Ty_vec =     trailing_circle_centre_Y[stator_flow_vec_pos]

# then to find the line equation in a form of Ax + By + C = 0
C1x    = T1_list[stator_flow_vec_pos+1][0]
C1y    = T1_list[stator_flow_vec_pos+1][1]

C2x    = L1_list[stator_flow_vec_pos+1][0]
C2y    = L1_list[stator_flow_vec_pos+1][1]

#ax.plot([C1x,C2x],[C1y,C2y], 'b-')
# create the line from the centre to the next common tangent line. 
# 已知点的坐标（x0，y0），直线上的两点（x1，y1）、（x2，y2）；求点在直线上的垂足（x, y。
# 垂足：首先，求一系数 k： 设直线的起点和终点分别为A（x1， y1）、B（x2， y2），直线外一点为C（x0， y0），垂足为D；并设k = |AD| / |AB|。
# 则，k * AB = AD = AC + CD，又 AB * CD= 0；所以，k * AB* AB = AC *AB，故 k =AC * AB / （AB * AB）。
# 带入坐标，即得， k = ( (x0- x1) * (x2 - x1) + (y0 - y1) * (y2 - y1) )  / ( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) ) ;
# 则 x = x1 + k*(x2 - x1); y = y1 + k*(y2 - y1);
k = ( (Tx_vec- C1x) * (C2x - C1x) + (Ty_vec - C1y) * (C2y - C1y))/ ((C2x - C1x)**2 + (C2y - C1y)**2)
Dx = C1x + k*(C2x - C1x)
Dy = C1y + k*(C2y - C1y)


# print the line that from the traling edge to other side of the throat
ax.plot([Tx_vec, Dx], [Ty_vec, Dy], 'b-' ,lw = lineWidth ) 

# find the middle point fo the throat line
Mx  =  0.5* (Tx_vec + Dx)
My  =  0.5* (Ty_vec + Dy)

# the arrow length is same as the common tangent line.
# should be adjust with different velocity vector 
arrowL = C3 / 2.e4
#arrowL = np.sqrt(  (C2x - C1x)**2  + (C2y - C1y)**2   )

# the arrow vector is parallel to the common tangent line
kCommon = (C2y - C1y)/(C2x - C1x)

Vx = Mx - np.sqrt(arrowL**2/(kCommon**2 +1 ))
Vy = np.sqrt(arrowL**2*kCommon**2/(kCommon**2 +1 )) + My

# now plot the arrow of the C3 here.
ax.arrow(Mx, My, Vx-Mx, Vy-My, head_width=0.0005, head_length=0.001,width = 0.0001, length_includes_head=True,fc='g', ec='g')


# now plot the arrow of the C3R here.
k = ((Vx - Mx)*(0. - Mx) + (Vy - My)*(0. - My)) / ((0.- Mx) **2. + (0.-My)**2.)
C3Rx = Mx + k * ( 0.- Mx)
C3Ry = My + k * ( 0.- My)
ax.arrow(Mx, My, C3Rx-Mx, C3Ry-My, head_width=0.0005, head_length=0.001,width = 0.0001, length_includes_head=True,fc='r', ec='r'  )
ax.arrow(C3Rx, C3Ry, Vx - C3Rx, Vy - C3Ry, head_width=0.0005, head_length=0.001,width = 0.0001, length_includes_head=True,fc='b', ec='b' )



###
#Start plotting the rotor inlet vectors here.
rarrowL = C4/ 2.e4

# Find the start point 
rtheta = np.pi * 2. / Zr

Lx, Ly = polar2cart( r4 , rtheta*(rotor_flow_vec_pos + 0.5))


# 求C4 箭头在C4R延长线上的垂足,也就是C4Rx, C4Ry

C4RL = np.abs(rarrowL * np.cos(alpha4) )
C4Rx, C4Ry = polar2cart( r4 - C4RL, rtheta*(rotor_flow_vec_pos + 0.5))

ax.arrow( Lx, Ly, C4Rx-Lx, C4Ry-Ly, head_width=0.0005, head_length=0.001,width = 0.0001, length_includes_head=True,fc='r', ec='r'  )



C4TL = np.abs(rarrowL * np.sin(alpha4) )


C4Tx =  C4Rx - np.sqrt( C4TL**2 /(  Lx**2/Ly**2 + 1. )  ) 
C4Ty = C4Ry - (C4Tx - C4Rx)*Lx/Ly 


ax.arrow( Lx, Ly, C4Tx-Lx, C4Ty-Ly, head_width=0.0005, head_length=0.001,width = 0.0001, length_includes_head=True,fc='g', ec='g'     )

ax.arrow( C4Rx, C4Ry, C4Tx - C4Rx, C4Ty - C4Ry, head_width=0.0005, head_length=0.001,width = 0.0001, length_includes_head=True,fc='b', ec='b'   )


# add the artist into the figure. 
ax.add_artist(circ_rotor_out_tip)
ax.add_artist(circ_rotor_out_hub)
ax.add_artist(circ_rotor)
ax.add_artist(circ_stator)
ax.add_artist(circ_inlet)

for i in range(len(theta_list)):
	ax.add_artist(trailing_circle_list[i])
	ax.add_artist(leading_circle_list[i])

if showQuater == 1:
    ax.set_xlim(0., r1*1.05 ) 
    ax.set_ylim(0., r1*1.05 )
else:
    ax.set_xlim(-r1*1.05,r1*1.05 )
    ax.set_ylim(-r1*1.05,r1*1.05)
plt.xlabel('size/[m]')
plt.ylabel('size/[m]')
#plt.xticks(fontsize = 16)
#plt.yticks(fontsize = 16)
#ax.ticklabel_format(scilimits= (1,1))


ax.text( 0.65*r1, 0.98*r1, '$r_3 = $' + str("%2.4f"%r3) + 'mm', fontsize = 20  )
ax.text( 0.65*r1, 0.9*r1, '$Z_3 = $' + str("%2.0f"%Zs) , fontsize = 20  )



plt.tight_layout()
plt.savefig('Plot_turbine_105loss' + str(Zs) +'.pdf')
plt.show()


