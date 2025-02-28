import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib.animation as animation
DATA_DIR = 
data_names = ['2msun_02_SI']#,'30msun_02_SI']

col_names_structure = ['Mr','r','Lr','P','rho',
                 'T','U','S','Cp','adia1',
                 'nablaad','mu','ne','Pe','Pr',
                 'nablarad','nabla','vc','kappa','epsnuc',
                 'epspp','epscno','eps3alpha','epsvnuc','epsv', 'epsgrav',
                 'X','-','X+','Y', 'Y+', 'Y++', 'XC', 'XN', 'XO', 'Psi']

col_names = ['Step','t','M','LogL','LogR',
                 'LogTs','LogTc','Logrhoc','logPc','psic',
                 'Xc','Yc','XCc','XNc','XOc',
                 'taudyn','tauKH','taunuc','Lpp','LCNO',
                 'L3a','LZ','Lv','MHe','MC',
                 'MO','RHe','RC','RO']

fig, ax = plt.subplots()
ax.set_xlabel(r"$r/R_\odot$ [1]")
ax.set_ylabel(r"Power per unit mass [$\mathrm{W/kg}$]")
ax.set_xscale('log')
x = np.arange(0, 2*np.pi, 0.01)
linepp, = ax.plot(x, np.sin(x), label="pp")
linecno, = ax.plot(x, np.sin(x), label="CNO")
linealpha, =  ax.plot(x, np.sin(x), label = "alpha")

def animate(i):
    data_name = data_files[i]
        
    #print(data_dir, data_name)
    file_path = os.path.join(full_path, data_name)
    
    df = pd.read_csv(file_path, sep=r'\s+', header=None)
    df.columns = col_names_structure
    
    r = df['r']

    
    linepp.set_xdata(r)
    linecno.set_xdata(r)
    linealpha.set_xdata(r)
    
    linepp.set_ydata(df['epspp'])
    linecno.set_ydata(df['epscno'])
    linealpha.set_ydata(df['eps3alpha'])
    allt = [max(df['epspp']), max(df['epscno']), max(df['eps3alpha'])]
    ax.set_ylim(0,max(allt))
    #ax.set_xlim(1e-8,1e2)
    
    ax.legend()
    ax.set_title(name + rf'\ time step {i+1}')
    return linepp, linecno, linealpha,
 


def get_structure_files(data_dir):
    all_files = os.listdir(data_dir)
    structure_files = [f for f in all_files if f.startswith("structure_")]
    
    
    return structure_files

def get_max_num(data_dir):
    all_files = os.listdir(data_dir)
    structure_files = [f for f in all_files if f.startswith("structure_")]
    time_steps = [int(f[10:15]) for f in structure_files]
    max_time_step = max(time_steps) if time_steps else 0
    
    return max_time_step

        

        
for name in data_names:
    
    full_path = DATA_DIR + name
    max_iter = get_max_num(full_path)
    data_files = get_structure_files(full_path)
    ani = animation.FuncAnimation(fig, animate, interval=0.5, blit=True, save_count=max_iter)

    ani.save(fr"test.gif")
    plt.show()