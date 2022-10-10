import numpy as np
from bfm.variable_info import set_variable_info_bfm

variable_abbrev, variable_names, variable_units, index, pelagic_index, seaice_index, benthic_index = set_variable_info_bfm()

def d3state_average_field(d3ave,d3state,case):
    if case == 'Initialize':
        d3ave.single_day_ave[:,:] = 0
        # d3ave.count = d3ave.count + 1

    elif case == 'Accumulate':
        d3ave.single_day_ave[:,:] = d3ave.single_day_ave[:,:] + d3state[:,:]
        # d3ave.count = d3ave.count + 1

    elif case == 'Mean':
        d3ave.single_day_ave[:,:] = d3ave.single_day_ave[:,:] + d3state[:,:]
        d3ave.single_day_ave[:,:] = d3ave.single_day_ave[:,:]/d3ave.count
        d3ave.daily_ave[:,:,d3ave.day] = d3ave.single_day_ave[:,:]
        d3ave.monthly_ave[:,:,d3ave.month] = d3ave.monthly_ave[:,:,d3ave.month] + d3ave.daily_ave[:,:,d3ave.day]
        if (d3ave.day != 0) & ((d3ave.day+1) % 30 == 0):
            # should actually be when (d3ave.day+1) % 30 == 0 because of python '0 first' indexing
            d3ave.monthly_ave[:,:,d3ave.month] = d3ave.monthly_ave[:,:,d3ave.month]/30
            d3ave.month += 1

        d3ave.day = d3ave.day + 1
        d3ave.single_day_ave[:,:] = 0
    
    elif case == 'Reset':
        # Output writing and print functions will go here
        d3ave.count = 0
        d3ave.day = 0
        d3ave.month = 0
        d3ave.single_day_ave[:,:] = 0
        d3ave.daily_ave[:,:,:] = 0
        d3ave.monthly_ave[:,:,:] = 0
        
    return d3ave 


def chl_average_field(chl_ave,chlorophylla,case):
    if case == 'Initialize':
        chl_ave.single_day_ave[:] = 0
        # d3ave.count = d3ave.count + 1

    elif case == 'Accumulate':
        chl_ave.single_day_ave[:] = chl_ave.single_day_ave[:] + chlorophylla[:]
        # d3ave.count = d3ave.count + 1

    elif case == 'Mean':
        chl_ave.single_day_ave[:] = chl_ave.single_day_ave[:] + chlorophylla[:]
        chl_ave.single_day_ave[:] = chl_ave.single_day_ave[:]/chl_ave.count
        chl_ave.daily_ave[:,chl_ave.day] = chl_ave.single_day_ave[:]
        chl_ave.monthly_ave[:,chl_ave.month] = chl_ave.monthly_ave[:,chl_ave.month] + chl_ave.daily_ave[:,chl_ave.day]
        if (chl_ave.day != 0) & ((chl_ave.day+1) % 30 == 0):
            # should actually be when (d3ave.day+1) % 30 == 0 because of python '0 first' indexing
            chl_ave.monthly_ave[:,chl_ave.month] = chl_ave.monthly_ave[:,chl_ave.month]/30
            chl_ave.month += 1

        chl_ave.day = chl_ave.day + 1
        chl_ave.single_day_ave[:] = 0
    
    elif case == 'Reset':
        # Output writing and print functions will go here
        chl_ave.count = 0
        chl_ave.day = 0
        chl_ave.month = 0
        chl_ave.single_day_ave[:] = 0
        chl_ave.daily_ave[:] = 0
        chl_ave.monthly_ave[:,:] = 0
        
    return chl_ave 


def npp_average_field(npp_ave,net_primary_production,case):
    if case == 'Initialize':
        npp_ave.single_day_ave[:] = 0
        # d3ave.count = d3ave.count + 1

    elif case == 'Accumulate':
        npp_ave.single_day_ave[:] = npp_ave.single_day_ave[:] + net_primary_production[:]
        # d3ave.count = d3ave.count + 1

    elif case == 'Mean':
        npp_ave.single_day_ave[:] = npp_ave.single_day_ave[:] + net_primary_production[:]
        npp_ave.single_day_ave[:] = npp_ave.single_day_ave[:]/npp_ave.count
        npp_ave.daily_ave[:,npp_ave.day] = npp_ave.single_day_ave[:]
        npp_ave.monthly_ave[:,npp_ave.month] = npp_ave.monthly_ave[:,npp_ave.month] + npp_ave.daily_ave[:,npp_ave.day]
        if (npp_ave.day != 0) & ((npp_ave.day+1) % 30 == 0):
            # should actually be when (d3ave.day+1) % 30 == 0 because of python '0 first' indexing
            npp_ave.monthly_ave[:,npp_ave.month] = npp_ave.monthly_ave[:,npp_ave.month]/30
            npp_ave.month += 1

        npp_ave.day = npp_ave.day + 1
        npp_ave.single_day_ave[:] = 0
    
    elif case == 'Reset':
        # Output writing and print functions will go here
        npp_ave.count = 0
        npp_ave.day = 0
        npp_ave.month = 0
        npp_ave.single_day_ave[:] = 0
        npp_ave.daily_ave[:] = 0
        npp_ave.monthly_ave[:,:] = 0
        
    return npp_ave 