import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
hfont = {'fontname':'Arial'}
plt.figure(figsize=(10, 5))#默认估计是6，5，
critical=0.462
file="./Critical_temp_charge.txt"
data = pd.read_csv(file,sep=' ', names=[ 'name', 'T','u','phi'])
data['num']=range(0,len(data['T']))
for i in range(len(data['T'])):
    print(data['name'][i],data['T'][i])
    xplus=1.6
    yplus=-0.02
    if i<1:
        color_choose="tab:red"
        xplus=0.0
        yplus=-0.04
    elif i<7:
        color_choose="tab:orange"
    elif i<13:
        color_choose="tab:green"        
    else :
        color_choose="tab:blue"        
        
    plt.annotate(data['name'][i], xy = (data['num'][i]+xplus, data['T'][i]/critical+yplus), fontsize=12 ,**hfont) # 这里xy是需要标记的坐标，xytext是对应的标签坐标
    data['num'][i]+=1
    plt.scatter(data['num'][i], data['T'][i]/critical,color=color_choose)


plt.axhline(y=1,ls="--", lw=1, c="tab:red")#添加水平直线
#plt.xticks(range(0,28),fontsize=0)
plt.xlim([-0.5,25])
plt.xticks([])
plt.ylabel("critical temp Tc",fontdict={ 'size' : 15})
plt.savefig('Tc.png')
plt.show()

             
