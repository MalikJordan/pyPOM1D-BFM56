import numpy as np
from functions import *

# ---------------------------------------------------------------------------------------------------------------------------------
# Data
model_name='pyPOM1D-BFM50'
pyPOM50 = load_python_data(model_name)
bfm50 = load_fortran_data()

# ---------------------------------------------------------------------------------------------------------------------------------
# Error Information
rmse_pyPOM50 = rmse(bfm50,pyPOM50)
std_bfm50 = standard_deviation(bfm50)
nrmse_pyPOM50 = nrmse(rmse_pyPOM50,std_bfm50)

# ---------------------------------------------------------------------------------------------------------------------------------
# Plots
plot_fields(pyPOM50,bfm50,model_name) # pyPOM50 set as 'check' for plotting to remain consistent with placement in other plots
plot_error(nrmse_pyPOM50,model_name,chl_removed=False)

# ---------------------------------------------------------------------------------------------------------------------------------
# Print Max NRMSE
species = ['Chl-a','Oxygen','Nitrate','Phosphate','PON','NPP','DIC']
maxes = np.max(nrmse_pyPOM50,axis=1)
avgs = np.average(nrmse_pyPOM50,axis=1)

print('Max NRMSE - pyPOM50 vs bfm50')
for i in range(0,7):
    print(species[i],' - ',maxes[i])
print()

print('Average NRMSE - pyPOM5- vs bfm50')
for i in range(0,7):
    print(species[i],' - ',avgs[i])
print()
