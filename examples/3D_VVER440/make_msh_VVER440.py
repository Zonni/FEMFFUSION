#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import division
from math import sqrt 
from collections import OrderedDict

msh_file = 'VVER440.msh' 

rows = [5, 10, 13, 14, 17, 18, 19, 20, 21, 20, 21, 22, 21,
        22, 21, 20, 21, 20, 19, 18, 17, 14, 13, 10, 5 ]


pitch = 14.7
node_list = []

ri = pitch / 2 # Radio circumferencia inscrita
rc = 2/3*sqrt(3) * ri # Radio circumferencia circumscrita == lado


# Points of the hexagon
p0 = [0.0, 0.0]
p1 = [0.0, - rc]
p2 = [+ri, -rc/2]
p3 = [+ri, +rc/2]
p4 = [0.0,  rc]
p5 = [-ri, +rc/2]
p6 = [-ri, -rc/2]
hex_part = [p0, p1, p2, p3, p4, p5, p6]

elements = []
lines = []


node_list = OrderedDict([])
# Make a row
n = 0
mx_rows = max(rows)
point_num = 0
for j in range(len(rows)):
    for i in range(rows[j]):
        
        p_num = 7*['a']
        # Create new points
        for num, p in enumerate(hex_part):
            node = (round(p[0] + ((mx_rows - rows[j])/2 + i )  * pitch, 3),
                    round(p[1] + j * 3/2 * rc, 3))
            if (node not in node_list):
                node_list[node] = point_num
                p_num[num] = point_num
                point_num += 1
            else :
                p_num[num] = node_list[node]
#        # Boundary Lines 
#        if i == 0:
 #           lines.append([0, ])
  #      else if i == rows[j]-1:
        
        # New element    
        elements.append([n, p_num[0], p_num[1], p_num[2], p_num[3]])
        elements.append([n, p_num[0], p_num[3], p_num[4], p_num[5]])
        elements.append([n, p_num[0], p_num[5], p_num[6], p_num[1]])
        n += 1

#print 'node_list and ', len(node_list) 
#node_list = list(OrderedDict.fromkeys(node_list))
#print 'node_list and ', len(node_list) 


 
fp = open(msh_file, "w")
fp.write("$MeshFormat\n")
fp.write("2.2 0 8\n")
fp.write("$EndMeshFormat\n")
fp.write("$Nodes\n")
fp.write(str(len(node_list)) + "\n")
for i, node in enumerate(node_list):
    fp.write(str(i+1) + " " + str(node[0]) + " " + str(node[1]) + " 0 \n")    

fp.write("$EndNodes\n")
fp.write("$Elements\n")
fp.write(str(len(lines) + len(elements)) + "\n")
i = 0
for ln in lines:
    fp.write(str(i) + " 1 2 " + str(i) + " " +
             str(ln[0]) + " " + str(i) + " " +
             str(ln[1])  + " " + str(ln[2]) + " \n")
    i += 1

for el in elements:
    for num in range(1, len(el)):
        assert(el[num] < len(node_list))
    fp.write(str(i+1) + " 3 2 " +
             str(el[0]+1) + " " +  str(i+1) + " "+
             str(el[1]+1)  + " " + str(el[2]+1)  + " " +
             str(el[3]+1) + " " +  str(el[4]+1) + " \n")
    i += 1 


fp.write("$EndElements")
fp.close()
    


#########################
#########  PLOT  ########
#########################
x = []
y = []
for l in node_list:
    x.append(l[0])
    y.append(l[1])
    
import matplotlib.pyplot as plt
from matplotlib import rcParams
plt.close('all')
plt.style.use('default')
params = {'backend': 'pgf',
          'pgf.rcfonts': False,
          'font.family': 'serif',
          'font.size': 14,
          'axes.labelsize': 14,
          'legend.fontsize': 12,
          'xtick.labelsize': 11,
          'ytick.labelsize': 11,
#          'text.usetex': True,
          'lines.linewidth': 3,
          'lines.markersize': 5,
          'lines.markeredgewidth': 1,
          'legend.numpoints': 1, 
          'figure.autolayout': True,
          }
rcParams.update(params)


fig, ax = plt.subplots()
ax.plot(x, y,  linestyle='', marker='o')
plt.axis('equal')