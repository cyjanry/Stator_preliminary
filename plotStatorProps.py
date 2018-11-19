#!/usr/bin/env python
# This Python file uses the following encoding: utf-8

import numpy      as         np
import scipy      as sci
import matplotlib.pyplot as plt 


Fig_tt_alpha3b_Zs   =  1
Fig_r3_alpha3b_


f = open('output', 'r')
dataList = []
for line in f:
    temp = line.replace('\n','').split(' ')
    del temp[-1]
    for i in range(len(temp)):
        temp[i] = float(temp[i])
    dataList.append(temp)
f.close()

# stator properties

ttlist = []
Zslist = []
r3list = []
A3devb3list = []
C3list = []
alpha3list =[]
alpha3blist = []


# Reading the list
for i in range(len(dataList)):
    ttlist.append(dataList[i][0])     #tt
    Zslist.append(dataList[i][1])     #Zs
    r3list.append(dataList[i][3])     #r3
    A3devb3list.append(dataList[i][4])
    C3list.append(dataList[i][5])   
    alpha3list.append(dataList[i][6])
    alpha3blist.append(dataList[i][7]) 

Z3max = max(Zslist)
Z3min = max(Zslist)


Znumber = Z3max - Z3min +1

listlength = len(ttlist)
iternum = listlength / Znumber

ttZ = []
alpha3bZ = []

# Arrange the list based on Zs
# 这一步是将tt和alpha3b按不同的Zs分成一条条线。 
for i in range(Znumber -1):
    temptt = []
    temp3b = []
    for j in range(iternum -1):

        temptt.append(ttlist[j * Znumber + i])
        temp3b.append(alpha3blist[j * Znumber + i])
    ttZ.append(temptt)
    alpha3bZ.append(temp3b)


# 图1是x: tt, y: alpha3b
color_list = ['k','y','b','g','r','c','m','k','y','b','g','r','c','m','k','y','b','g','r','c','m']
fig = plt.figure(figsize = (8,5))
ax = fig.add_subplot(111)
for k in range(len(ttZ)):
    ax.plot(ttZ[k],alpha3bZ[k],color_list[k]+'-')
ax.set_ylabel('Blade angle, $\\alpha_{3b}$ [$^\\circ$]')
ax.set_ylim(78,86)
ax.set_xlabel('Trailing edge thickness, $t$ [m]')
ax.grid()
plt.tight_layout()
plt.show()
