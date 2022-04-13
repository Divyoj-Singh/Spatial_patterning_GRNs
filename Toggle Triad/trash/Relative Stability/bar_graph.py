# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 13:00:49 2022

@author: csb
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(context='notebook', style='ticks',font='arial',font_scale=1.2)


#%% Bistable
data=pd.read_csv('Bistable.csv') 


f, ax = plt.subplots(figsize=(5,4))
x_pos=["100","010"];
mean=data['Mean'].to_list();
std=data['Std'].to_list()

ax.bar(x_pos, mean, yerr=std, align='center', ecolor='black',capsize=10)

#ax.set_aspect("equal")  
ax.set_xlabel("State",weight='bold')
ax.set_ylabel("Frequency",weight='bold')
ax.set_xticks(x_pos)
f.subplots_adjust(top=0.90, bottom=0.20, left=0.20, right=0.95, hspace=0.5,wspace=0.5)     
f.savefig("Bistable"+".png",dpi=800)
plt.close()



#%% Tristable

data=pd.read_csv('Tristable.csv') 

f, ax = plt.subplots(figsize=(5,4))
x_pos=["001","010","100"];
mean=data['Mean'].to_list();
std=data['Std'].to_list()

ax.bar(x_pos, mean, yerr=std, align='center', ecolor='black',capsize=10)

#ax.set_aspect("equal")  
ax.set_xlabel("State",weight='bold')
ax.set_ylabel("Frequency",weight='bold')
ax.set_xticks(x_pos)
f.subplots_adjust(top=0.90, bottom=0.20, left=0.20, right=0.95, hspace=0.5,wspace=0.5)     
f.savefig("Tristable"+".png",dpi=800)
plt.close()
