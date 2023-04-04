import numpy as np
from functions import *

# ---------------------------------------------------------------------------------------------------------------------------------
# Data
model_name='pyPOM1D-BFM21'
full_model='pyPOM1D-BFM50'
pyPOM21 = load_python_data(model_name)
pyPOM50 = load_python_data(full_model)

# ---------------------------------------------------------------------------------------------------------------------------------
# Error Information
rmse_pyPOM21 = rmse(pyPOM50,pyPOM21)
std_pyPOM50 = standard_deviation(pyPOM50)
nrmse_pyPOM21 = nrmse(rmse_pyPOM21,std_pyPOM50)

# ---------------------------------------------------------------------------------------------------------------------------------
# Plots
plot_fields(pyPOM50,pyPOM21,model_name) # pyPOM50 set as 'check' for plotting to remain consistent with placement in other plots
plot_error(nrmse_pyPOM21,model_name,chl_removed=True)

# ---------------------------------------------------------------------------------------------------------------------------------
# Print Max NRMSE
species = ['Chl-a','Oxygen','Nitrate','Phosphate','PON','NPP','DIC']
maxes = np.max(nrmse_pyPOM21,axis=1)

print('Max NRMSE - pyPOM21 vs pyPOM50')
for i in range(0,7):
    print(species[i],' - ',maxes[i])


# ---------------------------------------------------------------------------------------------------------------------------------
# Data
model_name='pyPOM1D-BFM21'
comp_model='pyPOM1D-BFM35'
pyPOM21 = load_python_data(model_name)
pyPOM35 = load_python_data(comp_model)

# ---------------------------------------------------------------------------------------------------------------------------------
# Error Information
rmse_pyPOM21 = rmse(pyPOM35,pyPOM21)
std_pyPOM35 = standard_deviation(pyPOM35)
nrmse_pyPOM21_35 = nrmse(rmse_pyPOM21,std_pyPOM35)

# ---------------------------------------------------------------------------------------------------------------------------------
# Plots
plot_fields(pyPOM35,pyPOM21,model_name='pyPOM1D-BFM21_35') # pyPOM50 set as 'check' for plotting to remain consistent with placement in other plots
plot_error(nrmse_pyPOM21_35,model_name='pyPOM1D-BFM21_35',chl_removed=True)

# ---------------------------------------------------------------------------------------------------------------------------------
# Print Max NRMSE
species = ['Chl-a','Oxygen','Nitrate','Phosphate','PON','NPP','DIC']
maxes = np.max(nrmse_pyPOM21_35,axis=1)

print('Max NRMSE - pyPOM21 vs pyPOM35')
for i in range(0,7):
    print(species[i],' - ',maxes[i])