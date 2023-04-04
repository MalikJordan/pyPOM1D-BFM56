import numpy
import copy
import sys
import json
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
from bfm.bfm56.BFM56_rate_eqns import reduced_bfm56_rate_eqns


def calc_error_chl_1(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak chlorophyll concentration during a spring bloom
    This is for the sum of p1l + p2l + p3l + p4l
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of chl-a concentrations (p1l + p2l + p3l + p4l)
    chl_sum_conc_full = solution_full_model.y[13] + solution_full_model.y[18] + solution_full_model.y[22] + solution_full_model.y[26]
    chl_sum_conc_reduced = solution_reduced_model.y[13] + solution_reduced_model.y[18] + solution_reduced_model.y[22] + solution_reduced_model.y[26]

    # Set time range for searching for max value
    # January to March in year 6
    t_min = 86400*365*6
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    spring_bloom_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    spring_bloom_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))
    
    p1l = 13    # P1l index
    p2l = 18    # P2l index
    p3l = 22    # P3l index
    p4l = 26    # P4l index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            spring_bloom_full[i,(j_full-index_t_min_full)] = solution_full_model.y[p1l+(i*50),j_full] + solution_full_model.y[p2l+(i*50),j_full] \
                                                               + solution_full_model.y[p3l+(i*50),j_full] + solution_full_model.y[p4l+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            spring_bloom_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[p1l+(i*50),j_reduced] + solution_reduced_model.y[p2l+(i*50),j_reduced] \
                                                                        + solution_reduced_model.y[p3l+(i*50),j_reduced] + solution_reduced_model.y[p4l+(i*50),j_reduced]
    
    # slice list of chl concentrations
    # spring_bloom_full = chl_sum_conc_full[index_t_min_full:index_t_max_full]
    # spring_bloom_reduced = chl_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find peak chl-a concentration in spring bloom
    spring_bloom_peak_full = numpy.max(spring_bloom_full)
    spring_bloom_peak_reduced = numpy.max(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_peak_full - spring_bloom_peak_reduced)/spring_bloom_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if spring_bloom_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_chl_2(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak chlorophyll concentration during a spring bloom
    This is for the sum of p2l + p3l + p4l
    """

    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of chl-a concentrations (p1l + p2l + p3l + p4l)
    chl_sum_conc_full = solution_full_model.y[18] + solution_full_model.y[22] + solution_full_model.y[26]
    chl_sum_conc_reduced =  solution_reduced_model.y[18] + solution_reduced_model.y[22] + solution_reduced_model.y[26]

    # Set time range for searching for max value
    # January to March in year 6
    t_min = 86400*365*6
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # slice list of chl concentrations
    spring_bloom_full = chl_sum_conc_full[index_t_min_full:index_t_max_full]
    spring_bloom_reduced = chl_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find peak chl-a concentration in spring bloom
    spring_bloom_peak_full = max(spring_bloom_full)
    spring_bloom_peak_reduced = max(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_peak_full - spring_bloom_peak_reduced)/spring_bloom_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if spring_bloom_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_chl_3(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak chlorophyll concentration during a spring bloom
    This is for the sum of p3l and p4l
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of chl-a concentrations (p1l + p2l + p3l + p4l)
    chl_sum_conc_full = solution_full_model.y[22] + solution_full_model.y[26]
    chl_sum_conc_reduced = solution_reduced_model.y[22] + solution_reduced_model.y[26]

    # Set time range for searching for max value
    # January to March in year 6
    t_min = 86400*365*6
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # slice list of chl concentrations
    spring_bloom_full = chl_sum_conc_full[index_t_min_full:index_t_max_full]
    spring_bloom_reduced = chl_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find peak chl-a concentration in spring bloom
    spring_bloom_peak_full = max(spring_bloom_full)
    spring_bloom_peak_reduced = max(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_peak_full - spring_bloom_peak_reduced)/spring_bloom_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if spring_bloom_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_chl_4(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average chlorophyll concentration during a spring bloom
    This is for the sum of p1l + p2l + p3l + p4l
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of chl-a concentrations (p1l + p2l + p3l + p4l)
    chl_sum_conc_full = solution_full_model.y[13] + solution_full_model.y[18] + solution_full_model.y[22] + solution_full_model.y[26]
    chl_sum_conc_reduced = solution_reduced_model.y[13] + solution_reduced_model.y[18] + solution_reduced_model.y[22] + solution_reduced_model.y[26]

    # Set time range for searching for max value
    # January to March in year 6
    t_min = 86400*365*6
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # slice list of chl concentrations
    spring_bloom_full = chl_sum_conc_full[index_t_min_full:index_t_max_full]
    spring_bloom_reduced = chl_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find average chl-a concentration in spring bloom
    spring_bloom_avg_full = numpy.mean(spring_bloom_full)
    spring_bloom_avg_reduced = numpy.mean(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_avg_full - spring_bloom_avg_reduced)/spring_bloom_avg_full)
    
    # If peak is zero in reduced model, set error to zero
    if spring_bloom_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_chl_5(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average chlorophyll concentration during a spring bloom
    This is for the sum of p2l + p3l + p4l
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of chl-a concentrations (p1l + p2l + p3l + p4l)
    chl_sum_conc_full = solution_full_model.y[18] + solution_full_model.y[22] + solution_full_model.y[26]
    chl_sum_conc_reduced =  solution_reduced_model.y[18] + solution_reduced_model.y[22] + solution_reduced_model.y[26]

    # Set time range for searching for max value
    # January to March in year 6
    t_min = 86400*365*6
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # slice list of chl concentrations
    spring_bloom_full = chl_sum_conc_full[index_t_min_full:index_t_max_full]
    spring_bloom_reduced = chl_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find average chl-a concentration in spring bloom
    spring_bloom_avg_full = numpy.mean(spring_bloom_full)
    spring_bloom_avg_reduced = numpy.mean(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_avg_full - spring_bloom_avg_reduced)/spring_bloom_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if spring_bloom_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_chl_6(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average chlorophyll concentration during a spring bloom
    This is for the sum of p3l and p4l
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of chl-a concentrations (p1l + p2l + p3l + p4l)
    chl_sum_conc_full = solution_full_model.y[22] + solution_full_model.y[26]
    chl_sum_conc_reduced = solution_reduced_model.y[22] + solution_reduced_model.y[26]

    # Set time range for searching for max value
    # January to March in year 6
    t_min = 86400*365*6
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # slice list of chl concentrations
    spring_bloom_full = chl_sum_conc_full[index_t_min_full:index_t_max_full]
    spring_bloom_reduced = chl_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find average chl-a concentration in spring bloom
    spring_bloom_avg_full = numpy.mean(spring_bloom_full)
    spring_bloom_avg_reduced = numpy.mean(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_avg_full - spring_bloom_avg_reduced)/spring_bloom_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if spring_bloom_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_dic_7(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average dissolved inorganic carbon (DIC) concentration
    during the month of january in the sixth year of the ten year simulation
    DIC is 'O3c' and its index is 48
    """
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*6
    t_max = t_min + 86400*31

    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break

    # dic_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # dic_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 48  # O3c index
    # for i in range(0,num_boxes):
    #     dic_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     dic_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of DIC concentration
    # dic_full = dic_full[:,index_t_min_full:index_t_max_full]
    # dic_reduced = dic_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    dic_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    dic_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 48  # O3c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            dic_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            dic_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # Find average DIC concentration
    dic_avg_full = numpy.mean(dic_full)
    dic_avg_reduced = numpy.mean(dic_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(dic_avg_full - dic_avg_reduced)/dic_avg_full)
    if dic_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_dic_8(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in annual average dissolved inorganic carbon (DIC) concentration
    during the sixth year of the ten year simulation
    DIC is 'O3c' and its index is 48
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*6
    t_max = 86400*365*7
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # dic_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # dic_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 48  # O3c index
    # for i in range(0,num_boxes):
    #     dic_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     dic_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of DIC concentration
    # dic_full = dic_full[:,index_t_min_full:index_t_max_full]
    # dic_reduced = dic_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    dic_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    dic_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 48  # O3c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            dic_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            dic_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # Find peak chl-a concentration
    dic_avg_full = numpy.mean(dic_full)
    dic_avg_reduced = numpy.mean(dic_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(dic_avg_full - dic_avg_reduced)/dic_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if dic_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_dic_9(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in annual average dissolved inorganic carbon (DIC) concentration
    during eighth year of the ten year simulation
    DIC is 'O3c' and its index is 48
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8
    t_max = 86400*365*9
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # dic_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # dic_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 48  # O3c index
    # for i in range(0,num_boxes):
    #     dic_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     dic_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of DIC concentration
    # dic_full = dic_full[:,index_t_min_full:index_t_max_full]
    # dic_reduced = dic_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    dic_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    dic_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 48  # O3c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            dic_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            dic_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # Find average DIC concentration
    dic_avg_full = numpy.mean(dic_full)
    dic_avg_reduced = numpy.mean(dic_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(dic_avg_full - dic_avg_reduced)/dic_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if dic_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_dic_10(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak dissolved inorganic carbon (DIC) concentration
    during the eighth year of the ten year simulation
    DIC is 'O3c' and its index is 48
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8.5
    t_max = 86400*365*9.5
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # dic_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # dic_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 48  # O3c index
    # for i in range(0,num_boxes):
    #     dic_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     dic_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of DIC concentration
    # dic_full = dic_full[:,index_t_min_full:index_t_max_full]
    # dic_reduced = dic_reduced[:,index_t_min_reduced:index_t_max_reduced]

    dic_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    dic_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 48  # O3c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            dic_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            dic_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
    # Find peak DIC concentration    
    # if statement in case reduction automation removes all spcies
    if len(species_removed) < 50:
        dic_peak_full = numpy.max(dic_full)
        dic_peak_reduced = numpy.max(dic_reduced)
        
        # Compute the error
        error = numpy.abs(100*numpy.abs(dic_peak_full - dic_peak_reduced)/dic_peak_full)
        
    # if all species are removed then error = nan
    else:
        error = numpy.nan
    
    return error, solution_reduced_model


def calc_error_dic_11(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in time of peak dissolved inorganic carbon (DIC) concentration
    during the eighth year of the ten year simulation
    DIC is 'O3c' and its index is 48
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8.5
    t_max = 86400*365*9.5
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # dic_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # dic_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 48  # O3c index
    # for i in range(0,num_boxes):
    #     dic_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     dic_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of DIC concentration
    # dic_full = dic_full[:,index_t_min_full:index_t_max_full]
    # dic_reduced = dic_reduced[:,index_t_min_reduced:index_t_max_reduced]

    dic_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    dic_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 48  # O3c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            dic_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            dic_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
    # slice list of time
    time_full = solution_full_model.t[index_t_min_full:index_t_max_full]
    time_reduced = solution_reduced_model.t[index_t_min_reduced:index_t_max_reduced]
    
    # if statement in case reduction automation removes all spcies
    if len(species_removed) < 50:
        # Find the index for the peak DIC concentration
        # dic_peak_index_full = numpy.argmax(dic_full)
        # dic_peak_index_reduced = numpy.argmax(dic_reduced)
        dic_peak_index_full = numpy.argmax(dic_full,axis=1)
        dic_peak_index_reduced = numpy.argmax(dic_reduced,axis=1)
        
        # Find the time of the peak DIC
        dic_peak_time_full = time_full[dic_peak_index_full]
        dic_peak_time_reduced = time_reduced[dic_peak_index_reduced]
        
        # Compute the error
        err = numpy.zeros(num_boxes)
        for i in range(0,num_boxes):
            err[i] = numpy.abs(100*numpy.abs(dic_peak_time_full[i] - dic_peak_time_reduced[i])/dic_peak_time_full[i])
            
            # If time of peak is zero in reduced model, set error to zero
            if dic_peak_time_reduced[i] == 0.0:
                err[i] = 100

        error = numpy.max(err)
        # error = numpy.abs(100*numpy.abs(dic_peak_time_full - dic_peak_time_reduced)/dic_peak_time_full)
        
        # If time of peak is zero in reduced model, set error to zero
        # if dic_peak_time_reduced == 0.0:
        #     error = 100
            
    # if all species are removed then error = nan
    else:
        error = numpy.nan
    
    return error, solution_reduced_model


def calc_error_r6n_12(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in the peak particulate organic phosphate concentration 
    during the eighth year of the ten year simulation.
    DIC is 'R6n' and its index is 45
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8
    t_max = 86400*365*9
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    r6n_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    r6n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 45  # R6n index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            r6n_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            r6n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # slice list of r6n concentration
    # r6n_full = solution_full_model.y[45][index_t_min_full:index_t_max_full]
    # r6n_reduced = solution_reduced_model.y[45][index_t_min_reduced:index_t_max_reduced]
    
    # Find peak r6n concentration
    r6n_peak_full = numpy.max(r6n_full)
    r6n_peak_reduced = numpy.max(r6n_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(r6n_peak_full - r6n_peak_reduced)/r6n_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if r6n_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_r6n_13(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in the average particulate organic phosphate concentration 
    during the eighth year of the ten year simulation.
    DIC is 'R6n' and its index is 45
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8
    t_max = 86400*365*9
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    r6n_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    r6n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 45  # R6n index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            r6n_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            r6n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
    # slice list of r6n concentration
    # r6n_full = solution_full_model.y[45][index_t_min_full:index_t_max_full]
    # r6n_reduced = solution_reduced_model.y[45][index_t_min_reduced:index_t_max_reduced]
    
    # Find avergae r6n concentration
    r6n_avg_full = numpy.mean(r6n_full)
    r6n_avg_reduced = numpy.mean(r6n_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(r6n_avg_full - r6n_avg_reduced)/r6n_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if r6n_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_r6n_14(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in time of peak particulate organic nitrogen (PON) concentration
    during the eighth year of the ten year simulation
    PON is 'R6n' and its index is 45
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value in full model
    t_min = 86400*365*8
    t_max = 86400*365*9
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    r6n_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    r6n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 45  # R6n index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            r6n_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            r6n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # slice list of time
    time_full = solution_full_model.t[index_t_min_full:index_t_max_full]
    time_reduced = solution_reduced_model.t[index_t_min_reduced:index_t_max_reduced]
    
    # Find the index for the peak o2o concentration
    r6n_peak_index_full = numpy.argmax(r6n_full,axis=1)
    r6n_peak_index_reduced = numpy.argmax(r6n_reduced,axis=1)
    
    # Find the time of the peak o2o
    r6n_peak_time_full = time_full[r6n_peak_index_full]
    r6n_peak_time_reduced = time_reduced[r6n_peak_index_reduced]

    err = numpy.zeros(num_boxes)
    for i in range(0,num_boxes):
        err[i] = numpy.abs(100*numpy.abs(r6n_peak_time_full[i] - r6n_peak_time_reduced[i])/r6n_peak_time_full[i])

        if r6n_peak_time_reduced[i] == 0.0:
            err[i] = 100

    error = numpy.max(err)   
    return error, solution_reduced_model


# def calc_error_r6n_14(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
#     """ calculates error in time of peak particulate organic nitrogen (PON) concentration
#     during the eighth year of the ten year simulation
#     PON is 'R6n' and its index is 45
#     """
    
#     # Integrate reduced model
#     # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
#     conc_reduced = copy.copy(c0)
#     num_boxes = conc_reduced.shape[0]
#     multiplier = numpy.ones(conc_reduced.shape[1])
#     for index in species_removed.values():
#         conc_reduced[:,index] = 0.0
#         multiplier[index] = 0.0
#     solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
#     # Reshape concentration matrix
#     conc_reduced = conc_reduced.reshape((num_boxes,50))

#     # Set time range for searching for max value in full model
#     t_min_full = 86400*365*8
#     t_max_full = 86400*365*9
    
#     # Find index associated for t_min and t_max for slicing the list
#     for index,time in enumerate(solution_full_model.t):
#         if time >= t_min_full:
#             index_t_min_full = index
#             break
#     for index,time in enumerate(solution_full_model.t):
#         if time > t_max_full:
#             index_t_max_full = index
#             break
    
#     r6n_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
#     r6n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

#     index = 45  # R6n index
#     for i in range(0,num_boxes):
#         for j_full in range(index_t_min_full,index_t_max_full):
#             r6n_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
#         for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
#             r6n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
#     # Slice list of PON concentration and time for full model
#     # r6n_full = solution_full_model.y[45][index_t_min_full:index_t_max_full]
#     time_full = solution_full_model.t[index_t_min_full:index_t_max_full]
    
#     # Find the index for the peak PON concentration
#     # r6n_peak_index_full = numpy.argmax(r6n_full)
#     r6n_peak_index_full = numpy.argmax(r6n_full,axis=1)
    
#     # Find the time of the peak PON concentration
#     r6n_peak_time_full = time_full[r6n_peak_index_full]
    
#     # Set bounds for searching for reduced model peak to be +/- 6 months around r6n_peak_time_full
#     t_min_reduced = r6n_peak_time_full - (86400*365/2)
#     t_max_reduced = r6n_peak_time_full + (86400*365/2)
    
#     # Find index associated for t_min and t_max for slicing the list
#     for index,time in enumerate(solution_reduced_model.t):
#         if time >= t_min_reduced:
#             index_t_min_reduced = index
#             break
#     for index,time in enumerate(solution_reduced_model.t):
#         if time > t_max_reduced:
#             index_t_max_reduced = index
#             break
    
#     r6n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))
#     index = 45  # R6n index
#     for i in range(0,num_boxes):
#         for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
#             r6n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
#     # Slice list of PON concentration and time for reduced model
#     # r6n_reduced = solution_reduced_model.y[45][index_t_min_reduced:index_t_max_reduced]
#     time_reduced = solution_reduced_model.t[index_t_min_reduced:index_t_max_reduced]
    
#     # Find the indices for the peak PON concentration
#     r6n_peak_indices_reduced = argrelextrema(r6n_reduced, numpy.greater, order=10)
#     # if len(species_removed) == 10:
#     #     sys.exit(r6n_peak_indices_reduced, type(r6n_peak_indices_reduced), len(r6n_peak_indices_reduced))
    
    
#     # Find the time of the peak PON concentration for each index in r6n_peak_indices_reduced
#     r6n_peak_times_reduced = numpy.zeros(len(r6n_peak_indices_reduced[0]))
#     for i, peak_index in enumerate(r6n_peak_indices_reduced[0]):
#         # sys.exit(r6n_peak_indices_reduced)
#         r6n_peak_times_reduced[i] = time_reduced[peak_index]
    
#     # Compute the error
#     error = numpy.min(numpy.abs(r6n_peak_time_full - r6n_peak_times_reduced)/(365*86400))*100
    
#     return error, solution_reduced_model


def calc_error_o2o_15(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average oxygen concentration
    during the month of january in the 8th year of the ten year simulation
    DIC is 'O2o and its index is 0
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8
    t_max = t_min + 86400*31
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    

    # o2o_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # o2o_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 0
    # for i in range(0,num_boxes):
    #     o2o_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     o2o_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of o2o concentration
    # o2o_full = o2o_full[:,index_t_min_full:index_t_max_full]
    # o2o_reduced = o2o_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    o2o_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    o2o_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 0
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            o2o_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            o2o_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
    # Find average o2o concentration
    o2o_avg_full = numpy.mean(o2o_full)
    o2o_avg_reduced = numpy.mean(o2o_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(o2o_avg_full - o2o_avg_reduced)/o2o_avg_full)

    # If avg is zero in reduced model, set error to zero
    if o2o_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_o2o_16(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in the average oxygen concentration 
    during the eighth year of the ten year simulation.
    DIC is 'O2o' and its index is 0
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8
    t_max = 86400*365*9
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # o2o_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # o2o_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 0
    # for i in range(0,num_boxes):
    #     o2o_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     o2o_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of o2o concentration
    # o2o_full = o2o_full[:,index_t_min_full:index_t_max_full]
    # o2o_reduced = o2o_reduced[:,index_t_min_reduced:index_t_max_reduced]

    o2o_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    o2o_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 0
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            o2o_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            o2o_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]  

    # Find avergae o2o concentration
    o2o_avg_full = numpy.mean(o2o_full)
    o2o_avg_reduced = numpy.mean(o2o_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(o2o_avg_full - o2o_avg_reduced)/o2o_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if o2o_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_o2o_17(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in the peak oxygen concentration 
    during the eighth year of the ten year simulation.
    DIC is 'O2o' and its index is 0
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8.5
    t_max = 86400*365*9.5
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # o2o_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # o2o_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 0
    # for i in range(0,num_boxes):
    #     o2o_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     o2o_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]

    # # slice list of o2o concentration
    # o2o_full = o2o_full[:,index_t_min_full:index_t_max_full]
    # o2o_reduced = o2o_reduced[:,index_t_min_reduced:index_t_max_reduced]
   
    o2o_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    o2o_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 0
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            o2o_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            o2o_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]


    # Find peak o2o concentration
    o2o_peak_full = numpy.max(o2o_full)
    o2o_peak_reduced = numpy.max(o2o_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(o2o_peak_full - o2o_peak_reduced)/o2o_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if o2o_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_o2o_18(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in time of peak oxygen concentration
    during the eighth year of the ten year simulation
    oxygen is 'O2o' and its index is 0
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8.5
    t_max = 86400*365*9.5
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # o2o_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # o2o_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    o2o_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    o2o_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    # index = 0
    # for i in range(0,num_boxes):
    #     o2o_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     o2o_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]
    
    index = 0
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            o2o_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            o2o_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # slice list of o2o concentration
    # o2o_full = o2o_full[:,index_t_min_full:index_t_max_full]
    # o2o_reduced = o2o_reduced[:,index_t_min_reduced:index_t_max_reduced]
   
    # slice list of time
    time_full = solution_full_model.t[index_t_min_full:index_t_max_full]
    time_reduced = solution_reduced_model.t[index_t_min_reduced:index_t_max_reduced]
    
    # Find the index for the peak o2o concentration
    # o2o_peak_index_full = numpy.argmax(o2o_full)
    # o2o_peak_index_reduced = numpy.argmax(o2o_reduced)
    
    o2o_peak_index_full = numpy.argmax(o2o_full,axis=1)
    o2o_peak_index_reduced = numpy.argmax(o2o_reduced,axis=1)
    
    # Find the time of the peak o2o
    o2o_peak_time_full = time_full[o2o_peak_index_full]
    o2o_peak_time_reduced = time_reduced[o2o_peak_index_reduced]

    # o2o_peak_time_full = numpy.zeros_like(o2o_peak_index_full)
    # o2o_peak_time_reduced = numpy.zeros_like(o2o_peak_index_reduced)
    err = numpy.zeros(num_boxes)
    for i in range(0,num_boxes):
        # o2o_peak_time_full[i] = time_full[o2o_peak_index_full[i]]
        # o2o_peak_time_reduced[i] = time_full[o2o_peak_index_reduced[i]]
        err[i] = numpy.abs(100*numpy.abs(o2o_peak_time_full[i] - o2o_peak_time_reduced[i])/o2o_peak_time_full[i])

        if o2o_peak_time_reduced[i] == 0.0:
            err[i] = 100

    error = numpy.max(err)
    # # Compute the error
    # error = numpy.abs(100*numpy.abs(o2o_peak_time_full - o2o_peak_time_reduced)/o2o_peak_time_full)
    
    # # If time of peak is zero in reduced model, set error to zero
    # if o2o_peak_time_reduced == 0.0:
    #     error = 100
    
    return error, solution_reduced_model


def calc_error_pc_19(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average phytoplankton carbon concentration during a spring bloom
    This is for the sum of p1c + p2c + p3c + p4c
    This is for the spring bloom (Jan - March) in the 8th year of a 10 year simulation
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of phytoplankton carbon concentrations (p1l + p2l + p3l + p4l)
    # pc_sum_conc_full = solution_full_model.y[10] + solution_full_model.y[15] + solution_full_model.y[19] + solution_full_model.y[23]
    # pc_sum_conc_reduced = solution_reduced_model.y[10] + solution_reduced_model.y[15] + solution_reduced_model.y[19] + solution_reduced_model.y[23]

    # Set time range for searching for max value
    # January to March in year 8
    t_min = 86400*365*8
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break

    spring_bloom_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    spring_bloom_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    p1c = 10    # P1c index
    p2c = 15    # P2c index
    p3c = 19    # P3c index
    p4c = 23    # P4c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            spring_bloom_full[i,(j_full-index_t_min_full)] = solution_full_model.y[p1c+(i*50),j_full] + solution_full_model.y[p2c+(i*50),j_full] \
                                                               + solution_full_model.y[p3c+(i*50),j_full] + solution_full_model.y[p4c+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            spring_bloom_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[p1c+(i*50),j_reduced] + solution_reduced_model.y[p2c+(i*50),j_reduced] \
                                                                        + solution_reduced_model.y[p3c+(i*50),j_reduced] + solution_reduced_model.y[p4c+(i*50),j_reduced]
    
    # slice list of phytoplankton carbon concentrations
    # spring_bloom_full = pc_sum_conc_full[index_t_min_full:index_t_max_full]
    # spring_bloom_reduced = pc_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    # Find average phytoplankton carbon concentration in spring bloom
    spring_bloom_avg_full = numpy.mean(spring_bloom_full)
    spring_bloom_avg_reduced = numpy.mean(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_avg_full - spring_bloom_avg_reduced)/spring_bloom_avg_full)
    
    # If average is zero in reduced model, set error to zero
    if spring_bloom_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_pc_20(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak phytoplankton carbon concentration during a spring bloom
    This is for the sum of p1c + p2c + p3c + p4c
    This is for the spring bloom (Jan - March) in the 8th year of a 10 year simulation
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of phytoplankton carbon concentrations (p1l + p2l + p3l + p4l)
    # pc_sum_conc_full = solution_full_model.y[10] + solution_full_model.y[15] + solution_full_model.y[19] + solution_full_model.y[23]
    # pc_sum_conc_reduced = solution_reduced_model.y[10] + solution_reduced_model.y[15] + solution_reduced_model.y[19] + solution_reduced_model.y[23]

    # Set time range for searching for max value
    # January to March in year 8
    t_min = 86400*365*8
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # Slice list of phytoplankton carbon concentrations
    # spring_bloom_full = pc_sum_conc_full[index_t_min_full:index_t_max_full]
    # spring_bloom_reduced = pc_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    
    spring_bloom_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    spring_bloom_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    p1c = 10    # P1c index
    p2c = 15    # P2c index
    p3c = 19    # P3c index
    p4c = 23    # P4c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            spring_bloom_full[i,(j_full-index_t_min_full)] = solution_full_model.y[p1c+(i*50),j_full] + solution_full_model.y[p2c+(i*50),j_full] \
                                                               + solution_full_model.y[p3c+(i*50),j_full] + solution_full_model.y[p4c+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            spring_bloom_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[p1c+(i*50),j_reduced] + solution_reduced_model.y[p2c+(i*50),j_reduced] \
                                                                        + solution_reduced_model.y[p3c+(i*50),j_reduced] + solution_reduced_model.y[p4c+(i*50),j_reduced]
    
    # Find peak phytoplankton carbon concentration in spring bloom
    spring_bloom_peak_full = numpy.max(spring_bloom_full)
    spring_bloom_peak_reduced = numpy.max(spring_bloom_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(spring_bloom_peak_full - spring_bloom_peak_reduced)/spring_bloom_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if spring_bloom_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_pc_21(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in time of peak phytoplankton carbon concentration during a spring bloom
    This is for the sum of p1c + p2c + p3c + p4c
    This is for the spring bloom (Jan - March) in the 8th year of a 10 year simulation
    """
    
    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get sum of phytoplankton carbon concentrations (p1l + p2l + p3l + p4l)
    pc_sum_conc_full = solution_full_model.y[10] + solution_full_model.y[15] + solution_full_model.y[19] + solution_full_model.y[23]
    pc_sum_conc_reduced = solution_reduced_model.y[10] + solution_reduced_model.y[15] + solution_reduced_model.y[19] + solution_reduced_model.y[23]

    # Set time range for searching for max value
    # January to March in year 8
    t_min = 86400*365*8
    t_max = t_min + 86400*90
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    spring_bloom_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    spring_bloom_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    p1c = 10    # P1c index
    p2c = 15    # P2c index
    p3c = 19    # P3c index
    p4c = 23    # P4c index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            spring_bloom_full[i,(j_full-index_t_min_full)] = solution_full_model.y[p1c+(i*50),j_full] + solution_full_model.y[p2c+(i*50),j_full] \
                                                               + solution_full_model.y[p3c+(i*50),j_full] + solution_full_model.y[p4c+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            spring_bloom_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[p1c+(i*50),j_reduced] + solution_reduced_model.y[p2c+(i*50),j_reduced] \
                                                                        + solution_reduced_model.y[p3c+(i*50),j_reduced] + solution_reduced_model.y[p4c+(i*50),j_reduced]
    
    # Slice list of phytoplankton carbon concentrations and time data
    # spring_bloom_full = pc_sum_conc_full[index_t_min_full:index_t_max_full]
    # spring_bloom_reduced = pc_sum_conc_reduced[index_t_min_reduced:index_t_max_reduced]
    time_full = solution_full_model.t[index_t_min_full:index_t_max_full]
    time_reduced = solution_reduced_model.t[index_t_min_reduced:index_t_max_reduced]
    
    # Find the index for the peak phytoplankton carbon concentration
    # pc_peak_index_full = numpy.argmax(spring_bloom_full)
    # pc_peak_index_reduced = numpy.argmax(spring_bloom_reduced)
    pc_peak_index_full = numpy.argmax(spring_bloom_full,axis=1)
    pc_peak_index_reduced = numpy.argmax(spring_bloom_reduced,axis=1)
    
    # Find the time of the peak phytoplankton carbon
    pc_peak_time_full = time_full[pc_peak_index_full]
    pc_peak_time_reduced = time_reduced[pc_peak_index_reduced]
    
    # pc_peak_full = spring_bloom_full[:,pc_peak_index_full]
    # pc_peak_reduced = spring_bloom_reduced[:,pc_peak_index_reduced]
    
    pc_peak_full = numpy.zeros(num_boxes)
    pc_peak_reduced = numpy.zeros(num_boxes)
    for i in range(0,num_boxes):
        pc_peak_full[i] = spring_bloom_full[i,pc_peak_index_full[i]]
        pc_peak_reduced[i] = spring_bloom_reduced[i,pc_peak_index_reduced[i]]
    
    err1 = numpy.zeros(num_boxes)
    err2 = numpy.zeros(num_boxes)
    for i in range(0,num_boxes):
        # o2o_peak_time_full[i] = time_full[o2o_peak_index_full[i]]
        # o2o_peak_time_reduced[i] = time_full[o2o_peak_index_reduced[i]]
        err1[i] = numpy.abs(100*numpy.abs(pc_peak_time_full[i] - pc_peak_time_reduced[i])/pc_peak_time_full[i])
        err2[i] = numpy.abs(100*numpy.abs(pc_peak_full[i] - pc_peak_reduced[i])/pc_peak_full[i])
        if pc_peak_time_reduced[i] == 0.0:
            err1[i] = 100
        if pc_peak_reduced[i] == 0.0:
            err2[i] =100

    error_time_of_peak = numpy.max(err1)
    error_peak = numpy.max(err2)
    
    
    
    
    # Compute the error in the time of peak
    # error_time_of_peak = 100*numpy.abs(pc_peak_time_full - pc_peak_time_reduced)/pc_peak_time_full
    
    # Cheack magnitude of peak
    # pc_peak_full = spring_bloom_full[pc_peak_index_full]
    # pc_peak_reduced = spring_bloom_reduced[pc_peak_index_reduced]
    # error_peak = 100*numpy.abs(pc_peak_full - pc_peak_reduced)/pc_peak_full
    
    # Compare error in time of peak and peak magnitude
    # if peak magnitude has > than 80% then output this error
    if error_peak > 80:
        error = error_peak
    else:
        error = error_time_of_peak
    
    return error, solution_reduced_model
    
    
def calc_error_n1p_22(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average phoshate concentration in the eighth year of the simulation
    """

    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Get n1p concentration
    # n1p_full = solution_full_model.y[1]
    # n1p_reduced = solution_reduced_model.y[1]

    # Set time range for searching for max value
    t_min = 86400*365*8     # beginning of year 8
    t_max = 86400*365*9     # beginning of year 9

    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    n1p_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    n1p_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 1   # N1p index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            n1p_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            n1p_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]   

    # n1p_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # n1p_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 1   # N1p index
    # for i in range(0,num_boxes):
    #     n1p_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     n1p_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]
    
    # # slice list of phosphate concentrations
    # n1p_full = n1p_full[:,index_t_min_full:index_t_max_full]
    # n1p_reduced = n1p_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    # Find average chl-a concentration in spring bloom
    n1p_avg_full = numpy.mean(n1p_full)
    n1p_avg_reduced = numpy.mean(n1p_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(n1p_avg_full - n1p_avg_reduced)/n1p_avg_full)
    
    # If peak is zero in reduced model, set error to zero
    if n1p_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_n1p_23(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak phoshate concentration in the eighth year of the simulation
    """

    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8.5
    t_max = 86400*365*9.5
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break

    n1p_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    n1p_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 1   # N1p index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            n1p_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            n1p_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]    
    
    # n1p_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # n1p_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 1   # N1p index
    # for i in range(0,num_boxes):
    #     n1p_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     n1p_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]
    
    # # slice list of phosphate concentrations
    # n1p_full = n1p_full[:,index_t_min_full:index_t_max_full]
    # n1p_reduced = n1p_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    # Find peak o2o concentration
    n1p_peak_full = numpy.max(n1p_full)
    n1p_peak_reduced = numpy.max(n1p_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(n1p_peak_full - n1p_peak_reduced)/n1p_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if n1p_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model
    
    
def calc_error_n3n_24(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in average nitrate concentration in the eighth year of the simulation
    """

    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    n3n_reduced = copy.copy(c0)
    num_boxes = n3n_reduced.shape[0]
    multiplier = numpy.ones(n3n_reduced.shape[1])
    for index in species_removed.values():
        n3n_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, n3n_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    n3n_reduced = n3n_reduced.reshape((num_boxes,50))

    # Get sum of n1p concentrations
    # nit_sum_conc_full = solution_full_model.y[2]
    # nit_sum_conc_reduced = solution_reduced_model.y[2]

    # Set time range for searching for max value
    t_min = 86400*365*8     # beginning of year 8
    t_max = 86400*365*9     # beginning of year 9

    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # n3n_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # n3n_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    # index = 2   # N3n index
    # for i in range(0,num_boxes):
    #     n3n_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     n3n_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]
    
    # # slice list of n3n concentrations
    # n3n_full = n3n_full[:,index_t_min_full:index_t_max_full]
    # n3n_reduced = n3n_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    n3n_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    n3n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 2   # N3n index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            n3n_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            n3n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]

    # Find average n3n concentration in spring bloom
    n3n_avg_full = numpy.mean(n3n_full)
    n3n_avg_reduced = numpy.mean(n3n_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(n3n_avg_full - n3n_avg_reduced)/n3n_avg_full)
    
    # If peak is zero in reduced model, set error to zero
    if n3n_avg_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


def calc_error_n3n_25(species_removed, solution_full_model, t_span, c0, bfm_phys_vars):
    """ calculates error in peak nitrate concentration in the eighth year of the simulation
    """

    # Integrate reduced model
    # Rempove indicated species in concentration list and create a multiplier for rate eqns to "remove species"
    conc_reduced = copy.copy(c0)
    num_boxes = conc_reduced.shape[0]
    multiplier = numpy.ones(conc_reduced.shape[1])
    for index in species_removed.values():
        conc_reduced[:,index] = 0.0
        multiplier[index] = 0.0
    solution_reduced_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, multiplier, True, num_boxes, True), t_span, conc_reduced.ravel(), method='RK23')
    
    # Reshape concentration matrix
    conc_reduced = conc_reduced.reshape((num_boxes,50))

    # Set time range for searching for max value
    t_min = 86400*365*8.5
    t_max = 86400*365*9.5
    
    # Find index associated for t_min and t_max for slicing the list
    for index,time in enumerate(solution_full_model.t):
        if time >= t_min:
            index_t_min_full = index
            break
    for index,time in enumerate(solution_full_model.t):
        if time > t_max:
            index_t_max_full = index
            break
        
    for index,time in enumerate(solution_reduced_model.t):
        if time >= t_min:
            index_t_min_reduced = index
            break
    for index,time in enumerate(solution_reduced_model.t):
        if time > t_max:
            index_t_max_reduced = index
            break
    
    # n3n_full = numpy.zeros((num_boxes,solution_full_model.y.shape[1]))
    # n3n_reduced = numpy.zeros((num_boxes,solution_reduced_model.y.shape[1]))

    n3n_full = numpy.zeros((num_boxes,(index_t_max_full-index_t_min_full)))
    n3n_reduced = numpy.zeros((num_boxes,(index_t_max_reduced-index_t_min_reduced)))

    index = 2   # N3n index
    for i in range(0,num_boxes):
        for j_full in range(index_t_min_full,index_t_max_full):
            n3n_full[i,(j_full-index_t_min_full)] = solution_full_model.y[index+(i*50),j_full]
        for j_reduced in range(index_t_min_reduced,index_t_max_reduced):
            n3n_reduced[i,(j_reduced-index_t_min_reduced)] = solution_reduced_model.y[index+(i*50),j_reduced]
    
    # index = 2   # N3n index
    # for i in range(0,num_boxes):
    #     n3n_full[i,:] = solution_full_model.y[index+(i*50),:]
    #     n3n_reduced[i,:] = solution_reduced_model.y[index+(i*50),:]
    
    # slice list of n3n concentrations
    # n3n_full = n3n_full[:,index_t_min_full:index_t_max_full]
    # n3n_reduced = n3n_reduced[:,index_t_min_reduced:index_t_max_reduced]
    
    # Find peak n3n concentration
    n3n_peak_full = numpy.max(n3n_full)
    n3n_peak_reduced = numpy.max(n3n_reduced)
    
    # Compute the error
    error = numpy.abs(100*numpy.abs(n3n_peak_full - n3n_peak_reduced)/n3n_peak_full)
    
    # If peak is zero in reduced model, set error to zero
    if n3n_peak_reduced == 0.0:
        error = 100
    
    return error, solution_reduced_model


