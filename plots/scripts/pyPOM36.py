import numpy as np
from functions import *

# ---------------------------------------------------------------------------------------------------------------------------------
# Data
model_name='pyPOM1D-BFM36'
full_model='pyPOM1D-BFM50'
pyPOM36 = load_python_data(model_name)
pyPOM50 = load_python_data(full_model)

# ---------------------------------------------------------------------------------------------------------------------------------
# Error Information
rmse_pyPOM36 = rmse(pyPOM50,pyPOM36)
std_pyPOM50 = standard_deviation(pyPOM50)
nrmse_pyPOM36 = nrmse(rmse_pyPOM36,std_pyPOM50)

# ---------------------------------------------------------------------------------------------------------------------------------
# Plots
plot_fields(pyPOM50,pyPOM36,model_name) # pyPOM50 set as 'check' for plotting to remain consistent with placement in other plots
plot_error(nrmse_pyPOM36,model_name,chl_removed=False)

# ---------------------------------------------------------------------------------------------------------------------------------
# Print Max NRMSE
species = ['Chl-a','Oxygen','Nitrate','Phosphate','PON','NPP','DIC']
maxes = np.max(nrmse_pyPOM36,axis=1)
avgs = np.average(nrmse_pyPOM36,axis=1)

print('Max NRMSE - pyPOM36 vs pyPOM50')
for i in range(0,7):
    print(species[i],' - ',maxes[i])
print()

print('Average NRMSE - pyPOM36 vs pyPOM50')
for i in range(0,7):
    print(species[i],' - ',avgs[i])
print()
