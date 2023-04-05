from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import brewer2mpl
import netCDF4 as nc
import numpy as np
import os

# matplotlib.use("TkAgg")

def force_aspect(ax,aspect=1):
    """Force plot aspects."""

    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)


def load_fortran_data():
    """Load fortran bfm56-pom1d data."""

    path = os.getcwd() + '/plots/model_data/bfm56_pom1d.nc'
    variables = nc.Dataset(path)
    variables = variables.variables

    # Extract fields of interest
    chlorophyll = variables['Chla'][:]
    oxygen = variables['O2o'][:]
    nitrate = variables['N3n'][:]
    phosphate = variables['N1p'][:]
    pon = variables['P1n'][:] + variables['P2n'][:] + variables['P3n'][:] + variables['P4n'][:] + variables['Z3n'][:] + variables['Z4n'][:] + variables['Z5n'][:] + variables['Z6n'][:] \
        + variables['R1n'][:] + variables['R6n'][:]
    production = (variables['ruPTc'][:] - variables['resPP'][:] - variables['resZT'][:])/12
    dic = (variables['DIC'][:])*(variables['ERHO'][:])*(12/1000)

    # Write as array
    chlorophyll = np.asarray(chlorophyll)
    oxygen = np.asarray(oxygen)
    nitrate = np.asarray(nitrate)
    phosphate = np.asarray(phosphate)
    pon = np.asarray(pon)
    production = np.asarray(production)
    dic = np.asarray(dic)

    # Load data into matrix
    data_fortran = np.zeros((7,150,3600))
    data_fortran[0,:,:] = chlorophyll.transpose()
    data_fortran[1,:,:] = oxygen.transpose()
    data_fortran[2,:,:] = nitrate.transpose()
    data_fortran[3,:,:] = phosphate.transpose()
    data_fortran[4,:,:] = pon.transpose()
    data_fortran[5,:,:] = production.transpose()
    data_fortran[6,:,:] = dic.transpose()

    # Calculate monthly averages for year 2 of simulation
    avg_data = np.zeros((7,150,12))
    for spec in range(0,7):
        for year in range(1,2):
            for month in range(0,12):
                for day in range(0,30):
                    avg_data[spec,:,month] = avg_data[spec,:,month] + data_fortran[spec,:,(day + (month*30) + (year*360))]
    avg_data = avg_data/30

    return avg_data


def load_python_data(model_name):
    """Load data from model of interest. Input model name as string."""

    path = os.getcwd() + '/plots/model_data/' + model_name + '.npz'
    model = np.load(path,allow_pickle=True)

    concentration = model['conc']
    chlorophyll = model['chl']
    production = model['npp']
    
    # Extract second year of data
    concentration = concentration[:,:,12:24]
    chlorophyll = chlorophyll[:,12:24]
    production = production[:,12:24]

    # Load data into matrix
    data = np.zeros((7,150,12))
    data[0,:,:] = chlorophyll   # Chlorophyll-a
    data[1,:,:] = concentration[:,0,:]  # Oxygen
    data[2,:,:] = concentration[:,2,:]  # Nitrate
    data[3,:,:] = concentration[:,1,:]  # Phosphate
    data[4,:,:] = concentration[:,11,:] + concentration[:,16,:] + concentration[:,20,:] + concentration[:,24,:] + concentration[:,28,:] \
                    + concentration[:,31,:] + concentration[:,34,:] + concentration[:,37,:] + concentration[:,40,:] + concentration[:,45,:] # Particulate Organic Nitrogen
    data[5,:,:] = production    # Net Primary Production
    data[6,:,:] = concentration[:,48,:] # Dissolved Inorganic Carbon

    return data


def rmse(check,comp):
    """Calculate the root mean square error between concentration matrices.
    RMSE of 'comp' with respect to 'check'."""

    er = check-comp
    se = np.power(er,2)
    mse = np.mean(se,axis=1)
    rmse = np.power(mse,0.5)

    return rmse


def standard_deviation(matrix):
    """Calculate standard deviation of concentration fields from input matrix.
    Returns a matrix containing the sandard deviation for each constituent."""

    mean = np.mean(matrix,axis=1)
    dif = np.zeros_like(matrix)
    for i in range(0,7):
        dif[i,:,:] = matrix[i,:,:] - mean[i,:]
    sq = np.power(dif,2)
    # sum = np.sum(sq,axis=1)
    sum = np.zeros(7)
    for i in range(0,7):
        sum[i] = np.sum(sq[i,:,:])/(150*12)
    std = np.power(sum,0.5)

    return std


def nrmse(rmse,std):
    """Calculate normalized root mean square error.
    Root mean square error normalized by the standard deviation of the check matrix."""

    nrmse = np.zeros_like(rmse)
    for i in range(0,7):
        nrmse[i,:] = rmse[i,:]/std[i]
    
    return nrmse


def plot_error(nrmse,model_name,chl_removed):
    """Create plot of normalized root mean square error."""

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Plot Style
    plt.rc('font', family='serif', size=14)
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.rc('axes', labelsize=20, linewidth=2)
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Legend Default
    plt.rc('legend', framealpha=1.0, facecolor='white', frameon=True, edgecolor='black')
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Plot Colors
    bmap = brewer2mpl.get_map('Paired', 'qualitative', 10)
    colors = bmap.mpl_colors
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Plot Details
    x = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5]
    legend = ['Chlorophyll-a','Oxygen','Nitrate','Phosphate','PON','NPP','DIC']
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Crop Data
    if chl_removed:
        nrmse = nrmse[0:6,:]
        legend = legend[0:6]
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Plot
    plt.figure(figsize=[12,8])
    for i in range(0,nrmse.shape[0]):
        # if i != 5:
        plt.scatter(x,nrmse[i],s=50)
    plt.xticks(x,['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.xlabel('Month')
    plt.xlim([0,12])
    plt.ylabel('NRMSE ($\%$)')
    plt.legend(legend,bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax = plt.gca()
    ax.set_ylim(bottom=0)
    ax.grid(linestyle='--')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_tick_params(which='major', size=7, width=2, direction='out')
    ax.yaxis.set_tick_params(which='major', size=7, width=2, direction='out')
    
    plt.tight_layout()
    fig_name = 'plots/figures/' + model_name + '_error.jpg'
    plt.savefig(fig_name)


def plot_fields(check,comp,model_name):
    """Create plots of concentration fields for the two input models (check and comp)."""
    
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Plot Style
    plt.rc('font', family='serif', size=20)
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    plt.rc('axes', labelsize=20, linewidth=1)
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Legend Default
    plt.rc('legend', framealpha=1.0, facecolor='white', frameon=True, edgecolor='black')
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Plot Colors
    bmap = brewer2mpl.get_map('Paired', 'qualitative', 10)
    colors = bmap.mpl_colors
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Titles
    title_check = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phoshate','(e) PON','(f) NPP','(g) DIC']
    title_comp = ['(h) Chl-a','(i) Oxygen','(j) Nitrate','(k) Phoshate','(l) PON','(m) NPP','(n) DIC']
    title = ['(a) Chl-a','(b) Oxygen','(c) Nitrate','(d) Phoshate','(e) Chl-a','(f) Oxygen','(g) Nitrate','(h) Phoshate','(i) PON','(j) NPP','(k) DIC','(l) PON','(m) NPP','(n) DIC']
    # ---------------------------------------------------------------------------------------------------------------------------------
    # Colorbar Limits
    clow   = [0,180,0,0,0.1,0,30]
    chigh  = [0.225,235,2.5,0.075,0.405,2.0,170]

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Field Plots - Style 1   
    fig,axes = plt.subplots(7,2,figsize=[11,30])
    for i in range(0,7):
        # Check
        plt.subplot(7,2,(2*i)+1)
        plt.imshow(check[i,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
        ax = plt.gca()
        plt.title(title_check[i])
        plt.xlabel('Month',fontsize=14)
        plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
        plt.ylabel('Depth (m)',fontsize=14)
        plt.yticks([0,50,100,150])
        plt.clim(clow[i],chigh[i])
        force_aspect(ax,aspect=1)

        # Comp
        plt.subplot(7,2,(2*i)+2)
        plt.imshow(comp[i,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
        ax = plt.gca()
        plt.title(title_comp[i])
        plt.xlabel('Month',fontsize=14)
        plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
        plt.yticks([0,50,100,150],[])
        plt.colorbar(orientation='vertical')
        plt.clim(clow[i],chigh[i])
        force_aspect(ax,aspect=1)


    plt.tight_layout()
    fig.subplots_adjust(wspace=-0.5,hspace=0.3)

    fig_name = 'plots/figures/' + model_name + '_fields_1.jpg'
    plt.savefig(fig_name)

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Field Plots - Style 2   
    fig,axes = plt.subplots(4,4,figsize=[16,15])
    for i in range(0,7):
        if i < 4:
            plt.subplot(4,4,i+1)
        else:
            plt.subplot(4,4,i+5)
        plt.imshow(check[i,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
        ax = plt.gca()
        plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
        plt.xlabel('Month',fontsize=14)
        if i%4 == 0:
            plt.yticks([0,50,100,150])
            plt.ylabel('Depth (m)',fontsize=14)
        else:
            plt.yticks([0,50,100,150],[])
        if i < 4:
            plt.title(title[i],fontsize=20)
        else:
            plt.title(title[i+4],fontsize=20)
        plt.clim(clow[i],chigh[i])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)

    for i in range(7,14):
        if i < 11:
            plt.subplot(4,4,i-2)
        else:
            plt.subplot(4,4,i+2)
        plt.imshow(comp[i-7,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
        ax = plt.gca()
        plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
        plt.xlabel('Month',fontsize=14)
        plt.yticks([0,50,100,150])
        if i%4 == 3:
            plt.yticks([0,50,100,150])
            plt.ylabel('Depth (m)',fontsize=14)
        else:
            plt.yticks([0,50,100,150],[])
        if i < 11:
            plt.title(title[i-3],fontsize=20)
        else:
            plt.title(title[i],fontsize=20)
        plt.clim(clow[i-7],chigh[i-7]) 
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)   

    fig.delaxes(axes[2,3])
    fig.delaxes(axes[3,3])

    plt.tight_layout(h_pad=0.75, w_pad=0.75)

    fig_name = 'plots/figures/' + model_name + '_fields_2.jpg'
    plt.savefig(fig_name)

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Field Plots - Style 3    
    fig,axes = plt.subplots(4,4,figsize=[16,15])
    for i in range(0,7):
        plt.subplot(4,4,i+1)
        plt.imshow(check[i,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
        ax = plt.gca()
        plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
        plt.xlabel('Month',fontsize=14)
        if i%4 == 0:
            plt.yticks([0,50,100,150])
            plt.ylabel('Depth (m)',fontsize=14)
        else:
            plt.yticks([0,50,100,150],[])
        plt.title(title_check[i],fontsize=20)
        plt.clim(clow[i],chigh[i])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)

    for i in range(7,14):
        plt.subplot(4,4,i+2)
        plt.imshow(comp[i-7,:,:],extent=[0,12,150,0],aspect='auto',cmap='jet')
        ax = plt.gca()
        plt.xticks([0.5,2.5,4.5,6.5,8.5,10.5], ['J','M','M','J','S','N'])
        plt.xlabel('Month',fontsize=14)
        plt.yticks([0,50,100,150])
        if i%4 == 3:
            plt.yticks([0,50,100,150])
            plt.ylabel('Depth (m)',fontsize=14)
        else:
            plt.yticks([0,50,100,150],[])
        plt.title(title_comp[i-7],fontsize=20)
        plt.clim(clow[i-7],chigh[i-7]) 
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)   

    fig.delaxes(axes[1,3])
    fig.delaxes(axes[3,3])

    plt.tight_layout(h_pad=0.75, w_pad=0.75)

    fig_name = 'plots/figures/' + model_name + '_fields_3.jpg'
    plt.savefig(fig_name)

