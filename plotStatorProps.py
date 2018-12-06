#!/usr/bin/env python
# This Python file uses the following encoding: utf-8

import numpy      as         np
import scipy      as sci
import matplotlib.pyplot as plt 


#Figure flags
TrailingEdgeThickness_bladeAngle_statorNumber_MachContours   =  0
TrailingEdgeThickness_bladeAngle_statorNumber_tpcontours     =  1


showfig = 0
savefig = 1
testflag = 0
stator_is_line = 1 # use lines or use contours to to create the iso-stator number line 

alpha3b_ylimit = (45, 80)

Fig_r3_alpha3b =  0


f = open('output', 'r')





linenumber = 0
dataList = []
for line in f:
    linenumber += 1
    if line.startswith('#'):
        temp = line.replace('\n',' ').split(' ')
        if linenumber ==4 :
            r3devr4min = float(temp[-2])
        if linenumber == 5:
            r3devr4max = float(temp[-2])
        if linenumber == 6:
            r3devr4step = float(temp[-2])
        if linenumber == 7:
            Z3min = float(temp[-2])
        if linenumber == 8:
            Z3max = float(temp[-2])
        if linenumber == 9:
            Z3step = float(temp[-2])
        if linenumber == 10:            
            ttmin = float(temp[-2])        
        if linenumber == 11:
            ttmax = float(temp[-2])
        if linenumber == 12:
            ttstep = float(temp[-2])
    else:
        temp = line.replace('\n','').split(' ')
        del temp[-1]
        for i in range(len(temp)):
            temp[i] = float(temp[i])
        dataList.append(temp)
f.close()




# stator properties
ttlist = []
Z3list = []
r3list = []
M3list = []
A3devb3list = []
C3list = []
alpha3list =[]
alpha3blist = []
T3list   = []
p3list   = []
rho3list = []

# Reading the list
for i in range(len(dataList)):
    r3list.append(dataList[i][0])      #r3
    Z3list.append(dataList[i][1])      #Zs
    ttlist.append(dataList[i][2])      #tt
    M3list.append(dataList[i][3])      #M3
    A3devb3list.append(dataList[i][4])
    C3list.append(dataList[i][5])   
    alpha3list.append(dataList[i][6])
    alpha3blist.append(dataList[i][7]) 
    p3list.append(dataList[i][10]/1.e6)
    T3list.append(dataList[i][11])
    rho3list.append(dataList[i][12])




ttnum = len(set(ttlist) )   #use set to return the distinct objects of a list
Z3num = len(set(Z3list) )
r3num = len(set(r3list))
color_list = ['k','y','b','g','r','c','m','k','y','b','g','r','c','m','k','y','b','g','r','c','m']

r3devr4list = np.arange(r3devr4min, r3devr4max + r3devr4step, r3devr4step)


# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()








# 这个图是绘制 trailing edge thicknees  VS alpha3b，在不同的stator number下，同时加入了Mach number contours
# This figure 
 
for i in range(r3num):
    # 这一步是将tt和alpha3b按不同的Zs分成一条条线。
    ttZ = []
    alpha3bZ = []
    alpha3Z  = []
    M3Z   =[]
    ZZ    = []
    pZ    = []
    TZ    = []
    rhoZ  = []
    for j in range(Z3num ):
        temptt = []
        temp3b = []
        tempM3 = []
        tempZ  = [] # this is used to plot the contours, rather than separated lines.
        temp3  = []
        tempp  = []
        tempT  = []
        temprho= []
        for k in range(ttnum ):
            temptt.append(     ttlist[ i*Z3num*ttnum  + j * ttnum + k  ])
            temp3b.append(alpha3blist[ i*Z3num*ttnum  + j * ttnum + k  ])
            tempM3.append(     M3list[ i*Z3num*ttnum  + j * ttnum + k  ])
            tempZ. append(     Z3list[ i*Z3num*ttnum  + j * ttnum + k  ])
            temp3. append( alpha3list[ i*Z3num*ttnum  + j * ttnum + k  ])
            tempp. append(     p3list[ i*Z3num*ttnum  + j * ttnum + k  ])
            tempT. append(     T3list[ i*Z3num*ttnum  + j * ttnum + k  ])
            temprho.append(  rho3list[ i*Z3num*ttnum  + j * ttnum + k  ])

        ttZ.     append(temptt)
        alpha3bZ.append(temp3b)
        alpha3Z. append(temp3)
        M3Z.     append(tempM3)
        ZZ.      append(tempZ)
        pZ.      append(tempp)
        TZ.      append(tempT)
        rhoZ.    append(temprho)


# 图1是x: tt, y: alpha3b, 不同的Zs曲线。 不同的r3/r4做不同的图
# 


    if 1 == TrailingEdgeThickness_bladeAngle_statorNumber_MachContours:
 
        fig = plt.figure(figsize = (8,5))
        ax = fig.add_subplot(111)

      
        if stator_is_line == 1:
            for a in range(len(ttZ)):
                if a % 2 ==0 and (a < (len(ttZ) - 1) ):
                    ax.plot(ttZ[a],alpha3bZ[a], 'k--')
                elif a % 2 ==0 : 
                	ax.plot(ttZ[a],alpha3bZ[a], 'k--', label= '$t$ vs $\\alpha_{3b}$')
           
            # Automatic generation of the line labels 
            coor_x     = ttZ[0][int(len(ttZ[0])/2)]
            coor_yZmin = alpha3bZ[ 0][int(len(ttZ[0])/2)]  *  1.01
            coor_yZmax = alpha3bZ[-1][int(len(ttZ[0])/2)]  *  0.90
            ax.text(coor_x , coor_yZmin, '$Z_s$ = ' + str(Z3list[ 0]), fontsize = 16  )
            ax.text(coor_x , coor_yZmax, '$Z_s$ = ' + str(Z3list[-1]), fontsize = 16  )

        else:
            ZCS   = ax.contour(ttZ, alpha3bZ, ZZ,
            	    9, linewidths= 2., colors = 'k')
            # Add a legend
            ZCS.collections[0].set_label('$Stator blade number$')
            ax.clabel(ZCS,  inline =2, fontsize = 16)

        # Plot the Mach number contours
        MachCS = ax.contour(ttZ, alpha3bZ, M3Z, 
                4, linewidths= 2., colors = 'k')

        #MachCS.levels = [nf(val) for val in MachCS.levels]
        #Add a legend
        MachCS.collections[0].set_label('$M_3$ contour')
        ax.clabel(MachCS, inline=2, fontsize=16, fmt = '%2.2f')
        

        # Plot the alpha3 contours
        Alpha3CS = ax.contour(ttZ, alpha3bZ, alpha3Z,
        		8, linewidths = 2., colors = 'r')
        #Add a legend

        #Alpha3CS.levels = [nf(val) for val in Alpha3CS.levels]
        Alpha3CS.collections[0].set_label('$\\alpha_3$ contour')
        ax.clabel(Alpha3CS, inline=2, fontsize=16, fmt = '%2.1f')        

        ax.plot([0.0009, 0.0009],[alpha3b_ylimit[0],alpha3b_ylimit[1] ], 'b--',lw = 2.,  label='$t$ empirical limit' )



        ax.set_ylabel('Stator blade angle, $\\alpha_{3b}$ [$^\\circ$]', fontsize = 16 )
        ax.set_xlabel('Trailing edge thickness, $t$ [m]', fontsize = 16 )
        ax.set_ylim(alpha3b_ylimit)
        ax.grid()

        ax.set_title("$r_3$/$r_4$ ={}".format(r3devr4list[i]), fontsize = 16 )

        plt.tight_layout()
        plt.legend(loc = 'best')
            
        if showfig == 1:
            plt.show()
        
        if savefig == 1:
            plt.savefig('loss_r3devr4_is_' + str(int(r3devr4list[i] *100)) +'.pdf' )
        
        if testflag ==1:
        	1/0
 


    elif 1 == TrailingEdgeThickness_bladeAngle_statorNumber_tpcontours:
        fig = plt.figure(figsize = (8,5))
        ax = fig.add_subplot(111)
        
        if stator_is_line == 1:
            for a in range(len(ttZ)):
                if a % 2 ==0 and (a < (len(ttZ) - 1) ):
                    ax.plot(ttZ[a],alpha3bZ[a], 'k--')
                elif a % 2 ==0 : 
                    ax.plot(ttZ[a],alpha3bZ[a], 'k--', label= '$t$ vs $\\alpha_{3b}$')
           
            # Automatic generation of the line labels 
            coor_x     = ttZ[0][int(len(ttZ[0])/2)]
            coor_yZmin = alpha3bZ[ 0][int(len(ttZ[0])/2)]  *  1.01
            coor_yZmax = alpha3bZ[-1][int(len(ttZ[0])/2)]  *  0.90
            ax.text(coor_x , coor_yZmin, '$Z_s$ = ' + str(Z3list[ 0]), fontsize = 16  )
            ax.text(coor_x , coor_yZmax, '$Z_s$ = ' + str(Z3list[-1]), fontsize = 16  )

        else:
            ZCS   = ax.contour(ttZ, alpha3bZ, ZZ,
                    9, linewidths= 2., colors = 'k')
            # Add a legend
            ZCS.collections[0].set_label('$Stator blade number$')
            ax.clabel(ZCS,  inline =2, fontsize = 16)

        # Plot the pressure contours
        PCS = ax.contour(ttZ, alpha3bZ, pZ, 
                4, linewidths= 2., colors = 'k')

        #Add a legend
        PCS.collections[0].set_label('$p_3$ contour [$MPa$]')
        ax.clabel(PCS, inline=2, fontsize=16, fmt = '%2.2f')
        

        # Plot the temperature contours
        TCS = ax.contour(ttZ, alpha3bZ, TZ,
                4, linewidths = 2., colors = 'r')
        #Add a legend
        TCS.collections[0].set_label('$T_3$ contour [$k$]')
        ax.clabel(TCS, inline=2, fontsize=16, fmt = '%2.1f')  


        # Plot the density contours      

        RCS = ax.contour(ttZ, alpha3bZ, rhoZ,
                4, linewidths = 2., colors = 'g')
        #Add a legend
        RCS.collections[0].set_label('$\\rho_3$ contour [$kg m^{-3}$]')
        ax.clabel(RCS, inline=2, fontsize=16, fmt = '%2.1f')  




        ax.plot([0.0009, 0.0009],[alpha3b_ylimit[0],alpha3b_ylimit[1] ], 'b--',lw = 2.,  label='$t$ empirical limit' )



        ax.set_ylabel('Stator blade angle, $\\alpha_{3b}$ [$^\\circ$]', fontsize = 16 )
        ax.set_xlabel('Trailing edge thickness, $t$ [m]', fontsize = 16 )
        ax.set_ylim(alpha3b_ylimit)
        ax.grid()

        ax.set_title("$r_3$/$r_4$ ={}".format(r3devr4list[i]), fontsize = 16 )

        plt.tight_layout()
        plt.legend(loc = 'best')
            
        if showfig == 1:
            plt.show()
        
        if savefig == 1:
            plt.savefig('loss_TPRHO_r3devr4_is_' + str(int(r3devr4list[i] *100)) +'.pdf' )
        
        if testflag ==1:
            1/0
