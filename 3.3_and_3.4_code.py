import os
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import cm
from matplotlib.colors import Normalize

DATA_DIR = 'Data/'
data_name = [r'2msun_02_SI\summary.txt',r'30msun_02_SI\summary.txt']


plt.style.use('matplotlib.mplstyle')

col_names = ['Step','t','M','LogL','LogR',
                 'LogTs','LogTc','Logrhoc','logPc','psic',
                 'Xc','Yc','XCc','XNc','XOc',
                 'taudyn','tauKH','taunuc','Lpp','LCNO',
                 'L3a','LZ','Lv','MHe','MC',
                 'MO','RHe','RC','RO']

col_names_structure = ['Mr','r','Lr','P','rho',
                 'T','U','S','Cp','adia1',
                 'nablaad','mu','ne','Pe','Pr',
                 'nablarad','nabla','vc','kappa','epsnuc',
                 'epspp','epscno','eps3alpha','epsvnuc','epsv', 'epsgrav',
                 'X','-','X+','Y', 'Y+', 'Y++', 'XC', 'XN', 'XO', 'Psi']

plot_title = [r'$2\mathrm{M_\odot}$',r'$30\mathrm{M_\odot}$']

#%%
#3.3 - nuclear processes

#The straight values for the luminosity-contributions

def plot_lum_values(data_files, data_dir):
    fig, axs = plt.subplots(2,2, sharex='col')
    
    for i, data_name in enumerate(data_files):
        
        df = pd.read_csv(data_dir + data_name, sep=r'\s+', header=None)
        df.columns = col_names
        df.drop(columns=['Step'], inplace=True)
        
        #Plotting the luminosity values normally 
        axs[0,i].plot(df['t'] / 1e9, df['Lpp'], label = r'$L_{\mathrm{pp}}$')
        axs[0,i].plot(df['t'] / 1e9, df['LCNO'], linestyle = '--', label = r'$L_{\mathrm{CNO}}$')
        axs[0,i].plot(df['t'] / 1e9, df['L3a'], linestyle = '-.', color = 'g', label = r'$L_{3\alpha}$')
        axs[0,i].set_ylabel(r"$L/$L$_\odot$")
        #axs[0,i].set_xlabel(r"$t \ [\mathrm{10^9 \ yr}]$")
        
        #Plotting on a log-scale for readability and comparison
        axs[1,i].semilogy(df['t'] / 1e9, df['Lpp'], label = r'$L_{\mathrm{pp}}$')
        axs[1,i].semilogy(df['t'] / 1e9, df['LCNO'], linestyle = '--', label = r'$L_{\mathrm{CNO}}$')
        axs[1,i].semilogy(df['t'] / 1e9, df['L3a'], linestyle = '-.', color = 'g', label = r'$L_{3\alpha}$')
        axs[1,i].set_ylabel(r"$L/$L$_\odot$")
        axs[1,i].set_xlabel(r"$t \ [\mathrm{10^9 \ yr}]$")
        
        
        axs[0,i].set_title(plot_title[i])
        axs[0,i].grid(which='both', alpha=0.4, visible=True)
        axs[1,i].grid(which='both', alpha=0.4, visible=True)
    
    axs[0,0].legend()
    
    plt.suptitle('Luminosity contribution over time for nuclear processes')
    plt.tight_layout()
    plt.show()
    
plot_lum_values(data_name, DATA_DIR)


#%%

#The proportion of luminosity contribution

def plot_lum_prop(data_files, data_dir, zoom = False):
    fig, axs = plt.subplots(2)
    
    if zoom == True:
        axs[0].set_xlim(1.38, 1.388)
        axs[1].set_xlim(0.005, 0.0062)
    
    for i, data_name in enumerate(data_files):
        df = pd.read_csv(data_dir + data_name, sep=r'\s+', header=None)
        df.columns = col_names
        df.drop(columns=['Step'], inplace=True)
        
        #Normalising values at each time-step
        for j in range(len(df['Lpp'])):
            #tot  = df.loc[j,'LogL']#
            tot = df['Lpp'][j] + df['LCNO'][j] + df['L3a'][j]
            df.loc[j, "Lpp"] = df.loc[j, "Lpp"]/tot#np.exp(tot)
            df.loc[j,'LCNO'] = df.loc[j,'LCNO']/tot#np.exp(tot)
            df.loc[j,'L3a'] = df.loc[j,'L3a']/tot#np.exp(tot)
            
        
        axs[i].plot(df['t'] / 1e9, df['Lpp'], label = r'$L_{\mathrm{pp}}$')
        
        axs[i].plot(df['t'] / 1e9, df['LCNO'], linestyle = '--', label = r'$L_{\mathrm{CNO}}$')
        
        axs[i].plot(df['t'] / 1e9, df['L3a'], linestyle = '-.', label = r'$L_{3\alpha}$')
        
        
        axs[i].set_ylabel(r"Fraction of total $L$")
        axs[1].set_xlabel(r"$t \ [\mathrm{10^9 \ yr}]$")
        
        
    
        axs[i].set_title(plot_title[i])
        axs[i].grid(which='both', alpha=0.4, visible=True)
    fig.legend([r'$L_{\mathrm{pp}}$',r'$L_{\mathrm{CNO}}$',r'$L_{3\alpha}$'], bbox_to_anchor=(0.75, 0.05), ncol=3)
    plt.suptitle('Proportion of the luminosity contribution')
    
    plt.tight_layout()
    plt.show()

plot_lum_prop(data_name, DATA_DIR)



#%%

#Data for power per unit mass

#time_steps_2 = [10, 380, 550, 830]
time_steps_2 = [10, 380, 580, 730, 920, 950]
time_steps_30 = [10, 140, 319, 524]

data_names = ['2msun_02_SI','30msun_02_SI']

def get_structure_files(data_dir, time_steps):
    df = pd.read_csv(data_dir + '/summary.txt', sep=r'\s+', header=None)
    df.columns = col_names
    all_files = os.listdir(os.path.join(data_dir))
    structure_files = [f for f in all_files if f.startswith("structure_")]
    
    # Load first file to extract time column
    df_ref = pd.read_csv(os.path.join(data_dir, structure_files[0]), sep=r'\s+', header=None)
    df_ref.columns = col_names_structure
    
    # Find indices closest to given time fractions
    #time_frac_idx = [np.argmin(np.abs(df['t'] - (t * np.max(df['t'])))) for t in time_frac]
    
    formatted_indices = [str(idx).zfill(5) for idx in time_steps]
    filtered_files = [f for f in structure_files if f[10:15] in formatted_indices]
    
    t = [df['t'].values[idx] for idx in time_steps]
    
    return sorted(filtered_files), t

def get_max_num(data_dir):
    all_files = os.listdir(data_dir)
    structure_files = [f for f in all_files if f.startswith("structure_")]
    time_steps = [int(f[10:15]) for f in structure_files]
    max_time_step = max(time_steps) if time_steps else 0
    
    return max_time_step

def plot_rad_contr(data_files, data_dir, max_iter, name, time):
    fig, axs = plt.subplots(len(data_files), sharex= True)
    
    cmap = cm.copper
    norm = Normalize(vmin=0, vmax=1)
    for i,data_name in enumerate(data_files):
        
        file_path = os.path.join(data_dir, data_name)
        
        df = pd.read_csv(file_path, sep=r'\s+', header=None)
        df.columns = col_names_structure
        
        r = df['r']
        
        axs[i].semilogx(r, df['epspp'], linestyle = '-', label = 'pp-chain')
        axs[2].set_ylabel(r"Power per unit mass [$\mathrm{Wkg^{-1}}$]")
        
        axs[i].semilogx(r, df['epscno'], linestyle = '--', label = 'CNO-cycle')
        
        
        axs[i].semilogx(r, df['eps3alpha'], linestyle = '-.', label = r'$3\alpha$-process')
        axs[len(data_files)-1].set_xlabel(r"$r$/R$_\odot$ [1]")
        #plt.colorbar(label =r'$t/t_\mathrm{end}$', aspect = 0.8)
        
        axs[i].annotate( rf'{time[i]/1e6:.2f} $\cdot$ $10^6$ yr', xy=(.71, .5), xycoords='axes fraction',bbox=dict(facecolor='none', edgecolor='grey'))
    axs[1].legend(bbox_to_anchor=(0.3, 0.3), fontsize = 4)
    plt.suptitle(name + rf'star')
    plt.subplots_adjust(hspace=0)#, wspace=0)
    plt.show()

        
for i, name in enumerate(data_names):
    steps = [time_steps_2,time_steps_30]
    full_path = DATA_DIR + name
    max_iter = get_max_num(full_path)
    data_files = get_structure_files(full_path, steps[i])[0]
    time = get_structure_files(full_path, steps[i])[1]
    plot_rad_contr(data_files, full_path, max_iter, plot_title[i], time)
    
    

#%%
#3.4 - Neutrinos

#Separate subplots for each star

def plot_neutrino(data_files, data_dir):
    fig, axs = plt.subplots(1,2)
    
    for i, data_name in enumerate(data_files):
        
        df = pd.read_csv(data_dir + data_name, sep=r'\s+', header=None)
        df.columns = col_names
        df.drop(columns=['Step'], inplace=True)
        
        
        axs[i].semilogy(df['t'] / 1e9, df['Lv'], label = plot_title[i])

        axs[i].set_ylabel(r"$L_{\nu}/$L$_\odot$")
        axs[i].set_xlabel(r"$t \ [\mathrm{10^9 \ yr}]$")
        
        axs[i].set_title(rf"{plot_title[i]}")
        axs[i].grid(which='both', alpha=0.4, visible=True)
        
    #fig.suptitle(r'Luminosity of neutrino losses, $L_{\mathrm{v}}$')
    plt.tight_layout()
    plt.show()

plot_neutrino(data_name, DATA_DIR)
#%%

#A single plot for both, with normalised time-axis

def plot_neutrino_together(data_files, data_dir):
    for i, data_name in enumerate(data_files):
        df = pd.read_csv(data_dir + data_name, sep=r'\s+', header=None)
        df.columns = col_names
        df.drop(columns=['Step'], inplace=True)
        
        
        plt.semilogy(df['t'] / max(df['t']), df['Lv'], label = plot_title[i])
        
        plt.ylabel(r"$L_{\nu}/$L$_\odot$")
        plt.xlabel(r"Normalised time")
        
        #plt.title(r'Luminosity of neutrino losses, $L_{\mathrm{v}}$')
        plt.grid(which='both', alpha=0.4, visible=True)
        plt.legend()
    plt.show()

plot_neutrino_together(data_name, DATA_DIR)
