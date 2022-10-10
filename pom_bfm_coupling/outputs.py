import numpy as np
# import netCDF4 as nc
from scipy.io import loadmat
from matplotlib import pyplot as plt
from pom.constants import current_path, num_boxes

def plot_outputs(d3ave, chl_ave, npp_ave):
    
    chlorophylla = chl_ave.monthly_ave[:,0:12]
    oxygen = d3ave.monthly_ave[:,0,0:12]
    phosphate = d3ave.monthly_ave[:,1,0:12]
    nitrate = d3ave.monthly_ave[:,2,0:12]
    parrticulate_organic_nitrogen = d3ave.monthly_ave[:,8,0:12] + d3ave.monthly_ave[:,16,0:12] + d3ave.monthly_ave[:,20,0:12] + d3ave.monthly_ave[:,24,0:12] + \
        d3ave.monthly_ave[:,28,0:12] + d3ave.monthly_ave[:,31,0:12] + d3ave.monthly_ave[:,34,0:12] + d3ave.monthly_ave[:,37,0:12] + d3ave.monthly_ave[:,40,0:12] + d3ave.monthly_ave[:,45,0:12]
    net_primary_production = npp_ave.monthly_ave[:,0:12]
    
    # oxygen = d3ave.daily_ave[:,0,0:11]
    # phosphate = d3ave.daily_ave[:,1,0:11]
    # nitrate = d3ave.daily_ave[:,2,0:11]
    # parrticulate_organic_nitrogen = d3ave.daily_ave[:,8,0:11] + d3ave.daily_ave[:,16,0:11] + d3ave.daily_ave[:,20,0:11] + d3ave.daily_ave[:,24,0:11] + \
    #     d3ave.daily_ave[:,28,0:11] + d3ave.daily_ave[:,31,0:11] + d3ave.daily_ave[:,34,0:11] + d3ave.daily_ave[:,37,0:11] + d3ave.daily_ave[:,40,0:11] + d3ave.daily_ave[:,45,0:11]

    x_labels = ['J','M','M','J','S','N']
    
    fig1, ax1 = plt.subplots()
    chla_plot = plt.imshow(chlorophylla,extent=[0,12,150,0],aspect='auto',cmap='bone')
    fig1.colorbar(chla_plot)
    # plt.clim(180,235)
    # plt.xlabel('Time (Days)')
    plt.xlabel('Month')
    plt.ylabel('Depth (m)')
    plt.title('Oxygen (mmol O2/m3)')
    ax1.set_xticks([0.5,2.5,4.5,6.5,8.5,10.5])
    ax1.set_xticklabels(x_labels)
    plt.show()


    fig2, ax2 = plt.subplots()
    oxygen_plot = plt.imshow(oxygen,extent=[0,12,150,0],aspect='auto',cmap='bone')
    fig2.colorbar(oxygen_plot)
    plt.clim(180,235)
    # plt.xlabel('Time (Days)')
    plt.xlabel('Month')
    plt.ylabel('Depth (m)')
    plt.title('Oxygen (mmol O2/m3)')
    ax2.set_xticks([0.5,2.5,4.5,6.5,8.5,10.5])
    ax2.set_xticklabels(x_labels)
    plt.show()

    fig3, ax3 = plt.subplots()
    phosphate_plot = plt.imshow(phosphate,extent=[0,12,150,0],aspect='auto',cmap='bone')
    fig2.colorbar(phosphate_plot)
    plt.clim(0,0.075)
    # plt.xlabel('Time (Days)')
    plt.xlabel('Month')
    plt.ylabel('Depth (m)')
    plt.title('Phosphate (mmol P/m3)')
    ax3.set_xticks([0.5,2.5,4.5,6.5,8.5,10.5])
    ax3.set_xticklabels(x_labels)
    plt.show()

    fig4, ax4 = plt.subplots()
    nitrate_plot = plt.imshow(nitrate,extent=[0,12,150,0],aspect='auto',cmap='bone')
    fig4.colorbar(nitrate_plot)
    plt.clim(0,2.5)
    # plt.xlabel('Time (Days)')
    plt.xlabel('Month')
    plt.ylabel('Depth (m)')
    plt.title('Nitrate (mmol N/m3)')
    ax4.set_xticks([0.5,2.5,4.5,6.5,8.5,10.5])
    ax4.set_xticklabels(x_labels)
    plt.show()

    fig5, ax5 = plt.subplots()
    PON_plot = plt.imshow(parrticulate_organic_nitrogen,extent=[0,12,150,0],aspect='auto',cmap='bone')
    fig5.colorbar(PON_plot)
    # plt.clim(0,0.075)
    # plt.xlabel('Time (Days)')
    plt.xlabel('Month')
    plt.ylabel('Depth (m)')
    plt.title('Particulate Organic Nitrate (mg N/m3)')
    ax5.set_xticks([0.5,2.5,4.5,6.5,8.5,10.5])
    ax5.set_xticklabels(x_labels)
    plt.show()

    fig6, ax6 = plt.subplots()
    NPP_plot = plt.imshow(net_primary_production,extent=[0,12,150,0],aspect='auto',cmap='bone')
    fig6.colorbar(NPP_plot)
    # plt.clim(0,0.075)
    # plt.xlabel('Time (Days)')
    plt.xlabel('Month')
    plt.ylabel('Depth (m)')
    plt.title('Particulate Organic Nitrate (mg N/m3)')
    ax6.set_xticks([0.5,2.5,4.5,6.5,8.5,10.5])
    ax6.set_xticklabels(x_labels)
    plt.show()

# def get_observational_data():

#     data_obs = np.zeros((6,150,12))

#     # Chlorophyll-a
#     chl = loadmat(current_path + '/inputs/Data_ML/Chla_1yr_climatology.mat')
#     chl = chl['Var']
#     data_obs[0,:,:] = chl[0:150,:]

#     # Oxygen
#     oxy = loadmat(current_path + '/inputs/Data_ML/Oxy_1yr_climatology.mat')
#     oxy = oxy['Var']
#     data_obs[1,:,:] = oxy[0:150,:]

#     # Nitrate
#     nit = loadmat(current_path + '/inputs/Data_ML/Nitrate_1yr_climatology.mat')
#     nit = nit['Var']
#     data_obs[2,:,:] = nit[0:150,:]

#     # Phosphate
#     phos = loadmat(current_path + '/inputs/Data_ML/Phos_1yr_climatology.mat')
#     phos = phos['Var']
#     data_obs[3,:,:] = phos[0:150,:]

#     # Particulate Organic Nitrogen
#     pon = loadmat(current_path + '/inputs/Data_ML/PON_1yr_climatology.mat')
#     pon = pon['Var']
#     data_obs[4,:,:] = pon[0:150,:]

#     # Net Primary Production
#     npp = loadmat(current_path + '/inputs/Data_ML/NPP_1yr_climatology.mat')
#     npp = npp['Var']
#     data_obs[5,:,:] = npp[0:150,:]

#     return data_obs

# def get_fortran_data():

#     data17 = nc.Dataset(current_path + '/inputs/bfm17_pom1d.nc')
#     data17 = data17.variables

#     data_fortran17 = np.zeros((6,150,3600))

#     # Chlorophyll a
#     chl17 = data17['Chla'][:]
#     chl17 = np.asarray(chl17)

#     # Oxygen
#     oxy17 = data17['O2o'][:]
#     oxy17 = np.asarray(oxy17)

#     # Nitrate
#     nit17 = data17['N3n'][:]
#     nit17 = np.asarray(nit17)

#     # Phosphate
#     phos17 = data17['N1p'][:]
#     phos17 = np.asarray(phos17)

#     # Particulate Organic Nitrogen
#     pon17 = data17['P2n'][:] + data17['Z5n'][:] + data17['R6n'][:]
#     pon17 = np.asarray(pon17)

#     # Net Primary Production
#     npp17 = data17['ruPTc'][:] + data17['resPP'][:] + data17['resZT'][:]
#     npp17 = np.asarray(npp17)

#     data_fortran17[0,:,:] = chl17.transpose()
#     data_fortran17[1,:,:] = oxy17.transpose()
#     data_fortran17[2,:,:] = nit17.transpose()
#     data_fortran17[3,:,:] = phos17.transpose()
#     data_fortran17[4,:,:] = pon17.transpose()
#     data_fortran17[5,:,:] = npp17.transpose()

#     data56 = nc.Dataset(current_path + '/inputs/bfm56_pom1d.nc')
#     data56 = data56.variables

#     data_fortran56 = np.zeros((6,150,3600))

#     # Chlorophyll a
#     chl56 = data56['Chla'][:]
#     chl56 = np.asarray(chl56)

#     # Oxygen
#     oxy56 = data56['O2o'][:]
#     oxy56 = np.asarray(oxy56)

#     # Nitrate
#     nit56 = data56['N3n'][:]
#     nit56 = np.asarray(nit56)

#     # Phosphate
#     phos56 = data56['N1p'][:]
#     phos56 = np.asarray(phos56)

#     # Particulate Organic Nitrogen
#     pon56 = data56['P1n'][:] + data56['P2n'][:] + data56['P3n'][:] + data56['P4n'][:] + data56['Z3n'][:] + data56['Z4n'][:] + data56['Z5n'][:] + data56['Z6n'][:] \
#             + data56['R1n'][:] + data56['R6n'][:]
#     pon56 = np.asarray(pon56)

#     # Net Primary Production
#     npp56 = data56['ruPTc'][:] + data56['resPP'][:] + data56['resZT'][:]
#     npp56 = np.asarray(npp56)

#     data_fortran56[0,:,:] = chl56.transpose()
#     data_fortran56[1,:,:] = oxy56.transpose()
#     data_fortran56[2,:,:] = nit56.transpose()
#     data_fortran56[3,:,:] = phos56.transpose()
#     data_fortran56[4,:,:] = pon56.transpose()
#     data_fortran56[5,:,:] = npp56.transpose()

#     avg_data_fortran17 = np. zeros((6,num_boxes,12))
#     avg_data_fortran56 = np. zeros((6,num_boxes,12))
#     for n in range(0,6):
#         for k in range(0,1):
#             for j in range(0,12):
#                 for m in range(0,30):
#                     avg_data_fortran17[n,:,j] = avg_data_fortran17[n,:,j] + data_fortran17[n,:, m + (j*30) + (k*360)]
#                     avg_data_fortran56[n,:,j] = avg_data_fortran56[n,:,j] + data_fortran56[n,:, m + (j*30) + (k*360)]

#     avg_data_fortran17 = avg_data_fortran17/30; 
#     avg_data_fortran56 = avg_data_fortran56/30; 

#     return avg_data_fortran17, avg_data_fortran56


def normalized_rmsd(matrix1,matrix2):
    """ Calculates the normalized root mean square difference of matrix2 to matrix1. Normalized by the mean of matrix1.
        return: 6x12 matrix with the normalized RMSD for the monthly averages of each constituent
    """
    constituents, gridpoints, months = np.shape(matrix1)
    nrmsd = np.zeros((constituents,months))

    for month in range(0,12):
        for i in range(0,6):
            # nrmsd[i,month] = np.sqrt((np.sum(np.power((matrix1[i,:,month]-matrix2[i,:,month]),2)))/np.size(matrix1[i,:,month]))/np.mean(matrix1[i,:,month])
            nrmsd[i,month] = np.sqrt((np.sum(np.power((matrix1[i,:,month]-matrix2[i,:,month]),2)))/gridpoints)/np.mean(matrix1[i,:,month])
    
    return nrmsd