from cppdefs import *
from include import POM_only
import numpy as np
import json
from pom.forcing import forcing_manager
from inputs import params_POMBFM
from pom.calculations import calculate_vertical_density_profile, create_vertical_coordinate_system
from pom.initialize_variables import get_temperature_and_salinity_initial_coditions
from pom.create_profiles import create_kinetic_energy_profile, create_vertical_diffusivity_profile, \
    calculate_vertical_temperature_and_salinity_profiles, calculate_vertical_zonal_velocity_profile, calculate_vertical_meridional_velocity_profile
from pom.data_classes import DiffusionCoefficients, ForcingManagerCounters, LeapFrogTimeLevels, MonthlyForcingData, Stresses, TemperatureSalinityData, VelocityData
from pom.constants import current_path, earth_angular_velocity, DAYI, water_specific_heat_times_density, vertical_layers, seconds_per_day, twice_the_timestep
from pom_bfm_coupling.initialize_variables import initialize_bfm_in_pom
from pom_bfm_coupling.data_classes import BfmPhysicalVariableData, D3stateAverageData, ChlAverageData, DicAverageData, NppAverageData, OutputData
from pom_bfm_coupling.coupling import pom_to_bfm, pom_bfm_1d, calculate_vertical_extinction, calculate_light_distribution
from matplotlib import pyplot as plt
import copy
import time
from reduction.modified_DRGEP import modified_DRGEP
from reduction.other_functions import solution_time
np.set_printoptions(precision=20)


with open(current_path + "/bfm/bfm56/bfm56_parameters.json") as read_parameters:
        parameters = json.load(read_parameters)
co2_flux_parameters = parameters["co2_flux_parameters"]

# GENERAL INITIALIZATION
length_scale = np.ones(vertical_layers)
length_scale[0] = 0.
length_scale[vertical_layers-1] = 0.
diffusion = DiffusionCoefficients()
kinetic_energy = LeapFrogTimeLevels(1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers))
kinetic_energy_times_length = LeapFrogTimeLevels(1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers),1.E-07 * np.ones(vertical_layers))
velocity = VelocityData()
wind_stress = Stresses()
bottom_stress = Stresses()
temperature = TemperatureSalinityData()
salinity = TemperatureSalinityData()

# DEFINE VERTICAL COORDINATE SYSTEM
vertical_grid = create_vertical_coordinate_system(params_POMBFM.kl1, params_POMBFM.kl2)
vertical_grid.length_scale = length_scale
# CORIOLIS PARAMETER
coriolis_parameter = 2. * earth_angular_velocity * np.sin(params_POMBFM.alat * 2. * np.pi / 360.)        # COR

# ITERATIONS NEEDED TO CARRY OUT AN "IDAYS" SIMULATION
iterations_needed = params_POMBFM.idays * seconds_per_day / params_POMBFM.dti                            # iend
# iterations_needed = 30 * seconds_per_day / params_POMBFM.dti                            # iend

# READ  T&S INITIAL CONDITIONS (IHOTST=0) OR RESTART FILE (IHOTST=1)
if params_POMBFM.ihotst == 0:
    t0 = 0.
    temperature.current, temperature.backward, salinity.current, salinity.backward = get_temperature_and_salinity_initial_coditions()
    # depth = vertical_grid.vertical_spacing * params_POMBFM.h
    # vertical_density_profile = calculate_vertical_density_profile(temperature.current, salinity.current, depth)
    vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_grid)
elif params_POMBFM.ihotst == 1:
    # get_rst()
    pass


REDUCE_BFM = False

if not POM_only:
    # INITIALIZATION OF BFM
    d3state, d3stateb = initialize_bfm_in_pom(vertical_grid)
    d3ave = D3stateAverageData()
    chl_ave = ChlAverageData()
    dic_ave = DicAverageData()
    npp_ave = NppAverageData()
    bfm_phys_vars = BfmPhysicalVariableData(vertical_layers-1)
    bfm_phys_vars.pH = co2_flux_parameters["ph_initial"]
    bfm_phys_vars.depth = vertical_grid.vertical_spacing[0:vertical_layers-1] * params_POMBFM.h
    outputs = OutputData()

    if REDUCE_BFM:       
        # INTEGRATION TIME OF FULL MODEL
        # integration_time_full = solution_time(d3state[0,:])

        # REDUCE MODEL
        # reduction, error_limit = modified_DRGEP(d3state)
        reduction, error_limit = modified_DRGEP(d3state,bfm_phys_vars)

        error_data = reduction['error_data']
        if error_data[-1] > error_limit:
            model = -2
        else:
            model = -1

        species_removed_data = reduction['species_removed_data']
        species_removed = species_removed_data[model]

        multiplier = np.ones(d3state.shape[1])
        d3state_reduced = copy.copy(d3state)
        for index in species_removed.values():
            d3state_reduced[:,index] = 0.0
            multiplier[index] = 0.0

        # UPDATE CONCENTRATION MATRICES
        d3state = copy.copy(d3state_reduced)
        d3stateb = copy.copy(d3state_reduced)

        # INTEGRATION TIME OF REDUCED MODEL
        # integration_time_reduced = solution_time(d3state[0,:])

        # PRINT 0D INTEGRATION TIMES
        num_species = d3state.shape[1] - len(species_removed)
        # print('Number of species in reduced model: ', num_species)
        # print('Full Model integration time: ', integration_time_full, 's')
        # print('Reduced Model integration time: ', integration_time_reduced, 's')

    else:
        multiplier = np.ones(d3state.shape[1])
        species_removed = []
else:
    bfm_phys_vars = BfmPhysicalVariableData()


# species_removed = [2,4,5,6,10,11,12,14,19,20,21,23,24,25,27,28,29,30,31,32,34,36,37,38,42,43,47,48,49]
# for i in range(0,len(species_removed)):
#     index = species_removed[i]
#     multiplier[index] = 0.0
#     d3state[:,index] = 0.0
#     d3stateb[:,index] = 0.0
# num_species = d3state.shape[1] - len(species_removed)
# print('num_species = ', num_species)
species_removed = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,20,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49] # BFM1
for i in range(0,len(species_removed)):
    index = species_removed[i]
    multiplier[index] = 0.0
    d3state[:,index] = 0.0
    d3stateb[:,index] = 0.0
num_species = d3state.shape[1] - len(species_removed)
print('num_species = ', num_species)

# bfm_phys_vars.pH = co2_flux_parameters["ph_initial"]
# BEGIN THE TIME MARCH
counters = ForcingManagerCounters()
month1_data = MonthlyForcingData()
month2_data = MonthlyForcingData()

# RECORD START TIME FOR COMPUTATION
start_time = time.process_time()
for i in range(0, int(iterations_needed)+1):

    t = t0 + (params_POMBFM.dti * i * DAYI)

    # TURBULENCE CLOSURE
    kinetic_energy.forward[:] = kinetic_energy.backward[:]
    kinetic_energy_times_length.forward[:] = kinetic_energy_times_length.backward[:]

    kinetic_energy, kinetic_energy_times_length, diffusion, vertical_grid = create_kinetic_energy_profile(vertical_grid, diffusion, temperature, salinity, vertical_density_profile, velocity,
                                                                                                          kinetic_energy, kinetic_energy_times_length, wind_stress, bottom_stress)
    # DEFINE ALL FORCINGS
    temperature.forward, temperature.interpolated, salinity.forward, salinity.interpolated, \
        shortwave_radiation, temperature.surface_flux, wind_stress, bfm_phys_vars.wgen, bfm_phys_vars.weddy, \
        month1_data, month2_data, counters, nutrients, inorganic_suspended_matter = forcing_manager(i,counters,month1_data,month2_data)

    if not POM_only:
        # ZERO FORCING DATA FOR REMOVED SPECIES
        if multiplier[0] == 0:  # Oxygen
            nutrients.O2bott = 0
        if multiplier[1] == 0:  # Phosphate
            nutrients.PO4bott = 0
            nutrients.PO4surf = 0
        if multiplier[2] == 0:  # Nitrate
            nutrients.NO3bott = 0
            nutrients.NO3surf = 0
        if multiplier[3] == 0:  # Ammonium
            nutrients.NH4surf = 0
            nutrients.PONbott_grad = 0
        if multiplier[5] == 0:  # Silicate
            nutrients.SIO4surf = 0

    # T&S COMPUTATION
    if params_POMBFM.idiagn == 0:
        # PROGNOSTIC MODE
        # T&S FULLY COMPUTED BY MODEL
        temperature.surface_value = temperature.forward[0]
        salinity.surface_value = salinity.forward[0]

        if params_POMBFM.trt != 0:
            for j in range(0, vertical_layers):
                if (-vertical_grid.vertical_coordinates_staggered[j] * params_POMBFM.h) >= params_POMBFM.upperh:
                    temperature.lateral_advection[j] = (temperature.interpolated[j] - temperature.current[j]) / (params_POMBFM.trt * seconds_per_day)

        if params_POMBFM.srt != 0:
            for j in range(0, vertical_layers):
                if (-vertical_grid.vertical_coordinates_staggered[j] * params_POMBFM.h) >= params_POMBFM.upperh:
                    salinity.lateral_advection[j] = (salinity.interpolated[j] - salinity.current[j]) / (params_POMBFM.srt * seconds_per_day)

        # COMPUTE SURFACE SALINITY FLUX
        salinity.surface_flux = -(salinity.surface_value - salinity.current[0]) * params_POMBFM.srt / seconds_per_day

        # COMPUTE TEMPREATURE
        temperature.forward[:] = temperature.backward[:] + (temperature.lateral_advection[:] * twice_the_timestep)
        calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, temperature, shortwave_radiation, params_POMBFM.nbct, params_POMBFM.umol)

        # CALCULATE SALINITY
        salinity.forward[:] = salinity.backward[:] + (salinity.lateral_advection[:] * twice_the_timestep)
        calculate_vertical_temperature_and_salinity_profiles(vertical_grid, diffusion, salinity, shortwave_radiation, params_POMBFM.nbcs, params_POMBFM.umol)

        # MIXING THE TIMESTEP (ASSELIN)
        temperature.current[:] = temperature.current[:] + 0.5 * params_POMBFM.smoth * (temperature.forward[:] + temperature.backward[:] - 2. * temperature.current[:])
        salinity.current[:] = salinity.current[:] + 0.5 * params_POMBFM.smoth * (salinity.forward[:] + salinity.backward[:] - 2. * salinity.current[:])

    velocity.zonal_forward[:] = velocity.zonal_backward[:] + twice_the_timestep * coriolis_parameter * velocity.meridional_current[:]
    velocity, bottom_stress = calculate_vertical_zonal_velocity_profile(vertical_grid, wind_stress, bottom_stress, diffusion, velocity)

    velocity.meridional_forward[:] = velocity.meridional_backward[:] - twice_the_timestep * coriolis_parameter * velocity.zonal_current[:]
    velocity, bottom_stress = calculate_vertical_meridional_velocity_profile(vertical_grid, wind_stress, bottom_stress, diffusion, velocity)

    # MIX TIME STEL (ASSELIN FILTER)
    kinetic_energy.current[:] = kinetic_energy.current[:] + 0.5 * params_POMBFM.smoth * (kinetic_energy.forward[:] + kinetic_energy.backward[:] - 2. * kinetic_energy.current[:])
    kinetic_energy_times_length.current[:] = kinetic_energy_times_length.current[:] + 0.5 * params_POMBFM.smoth * (kinetic_energy_times_length.forward[:] + kinetic_energy_times_length.backward[:] - 2. * kinetic_energy_times_length.current[:])

    velocity.zonal_current[:] = velocity.zonal_current[:] + 0.5 * params_POMBFM.smoth * (velocity.zonal_forward[:] + velocity.zonal_backward[:] - 2. * velocity.zonal_current[:])
    velocity.meridional_current[:] = velocity.meridional_current[:] + 0.5 * params_POMBFM.smoth * (velocity.meridional_forward[:] + velocity.meridional_backward[:] - 2. * velocity.meridional_current[:])

    # RESTORE TIME SEQUENCE
    kinetic_energy.backward[:] = kinetic_energy.current[:]
    kinetic_energy.current[:] = kinetic_energy.forward[:]
    kinetic_energy_times_length.backward[:] = kinetic_energy_times_length.current[:]
    kinetic_energy_times_length.current[:] = kinetic_energy_times_length.forward[:]

    velocity.zonal_backward[:] = velocity.zonal_current[:]
    velocity.zonal_current[:] = velocity.zonal_forward[:]
    velocity.meridional_backward[:] = velocity.meridional_current[:]
    velocity.meridional_current[:] = velocity.meridional_forward[:]

    temperature.backward[:] = temperature.current[:]
    temperature.current[:] = temperature.forward[:]
    salinity.backward[:] = salinity.current[:]
    salinity.current[:] = salinity.forward[:]

    # UPDATE DENSITY
    depth = vertical_grid.vertical_spacing * params_POMBFM.h
    # vertical_density_profile = calculate_vertical_density_profile(temperature.current,salinity.current,depth)
    vertical_density_profile = calculate_vertical_density_profile(temperature,salinity,vertical_grid)

    if not POM_only:
        bfm_phys_vars = pom_to_bfm(bfm_phys_vars, vertical_grid, temperature, salinity, inorganic_suspended_matter, shortwave_radiation, vertical_density_profile, wind_stress)
        bfm_phys_vars.vertical_extinction = calculate_vertical_extinction(bfm_phys_vars, d3state)
        bfm_phys_vars.irradiance = calculate_light_distribution(bfm_phys_vars)
        
        d3state, d3stateb, d3ave, chl_ave, dic_ave, npp_ave = pom_bfm_1d(i, vertical_grid, t, diffusion, nutrients, bfm_phys_vars, d3state, d3stateb, d3ave, chl_ave, dic_ave, npp_ave, multiplier)    

if POM_only:
    pom_data = [temperature,salinity,bottom_stress,velocity,kinetic_energy,kinetic_energy_times_length,vertical_density_profile]
    np.savez("pom_data",pyPOM1D=pom_data)
# Write Outputs
if not POM_only:
    # data_python = np.zeros((7,150,25))

    # data_python[0,:,:] = chl_ave.monthly_ave[:,:]   # Chl-a
    # data_python[1,:,:] = d3ave.monthly_ave[:,0,:]   # Oxygen
    # data_python[2,:,:] = d3ave.monthly_ave[:,2,:]   # Nitrate
    # data_python[3,:,:] = d3ave.monthly_ave[:,1,:]   # Phosphate
    # data_python[4,:,:] = d3ave.monthly_ave[:,11,:] + d3ave.monthly_ave[:,16,:] + d3ave.monthly_ave[:,20,:] + d3ave.monthly_ave[:,24,:] + \
    #     d3ave.monthly_ave[:,28,:] + d3ave.monthly_ave[:,31,:] + d3ave.monthly_ave[:,34,:] + d3ave.monthly_ave[:,37,:] + d3ave.monthly_ave[:,40,:] + d3ave.monthly_ave[:,45,:]   # PON
    # data_python[5,:,:] = npp_ave.monthly_ave[:,:]   # NPP
    # data_python[6,:,:] = dic_ave.monthly_ave[:,:]

    # np.savez("python_data.npz",bfm56=data_python)

    np.savez("pyPOM1D-BFM.npz",conc=d3ave.monthly_ave,chl=chl_ave.monthly_ave,dic=dic_ave.monthly_ave,npp=npp_ave.monthly_ave,species_removed=species_removed,integration_time=(time.process_time()-start_time))

# RECORD END TIME FOR COMPUTATION
end_time = time.process_time()
print('Time elapsed: ', end_time-start_time)

print('Main done')
