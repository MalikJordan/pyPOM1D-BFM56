import numpy
import logging
import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import brewer2mpl

from scipy.integrate import solve_ivp
from reduction.PythonBFM.BFM50_rate_eqns import bfm50_rate_eqns
from reduction.pyMARS_DRGEP_functions import get_importance_coeffs
import reduction.error_fcns as error_fcns
from reduction.other_functions import extract_dict, fill_dict
from bfm.bfm56.BFM56_rate_eqns import reduced_bfm56_rate_eqns

from pom.constants import vertical_layers
from pom_bfm_coupling.data_classes import BfmPhysicalVariableData

# Names of species in the system
species_names = ['O2o', 'N1p', 'N3n', 'N4n', 'O4n', 'N5s', 'N6r', 'B1c', 'B1n', 'B1p', 
                 'P1c', 'P1n', 'P1p', 'P1l', 'P1s', 'P2c', 'P2n', 'P2p', 'P2l',
                 'P3c', 'P3n', 'P3p', 'P3l', 'P4c', 'P4n', 'P4p', 'P4l',
                 'Z3c', 'Z3n', 'Z3p', 'Z4c', 'Z4n', 'Z4p', 'Z5c', 'Z5n', 'Z5p',
                 'Z6c', 'Z6n', 'Z6p', 'R1c', 'R1n', 'R1p', 'R2c', 'R3c', 'R6c', 
                 'R6n', 'R6p', 'R6s', 'O3c', 'O3h']


def calc_modified_DRGEP_dic(rate_eqn_fcn, conc, t, bfm_phys_vars):
    """ Calculates the percent difference between the new rate eqn and the 
    original. The new rate is based on turning one species 'off' """

    # c_original = list(conc)
    # c_new = list(conc)
    c_original = copy.copy(conc)
    c_new = copy.copy(conc)
    number_of_species = conc.shape[1]
    num_boxes = conc.shape[0]
    mult = numpy.ones(number_of_species)
    # calculate original rate values
    # dc_dt_og = rate_eqn_fcn(t, conc)
    # dc_dt_og = rate_eqn_fcn(0, bfm_phys_vars, t, conc, mult)
    # dc_dt_og = rate_eqn_fcn(t, c_new, bfm_phys_vars, mult, num_boxes, reducing=False)
    dc_dt_og = rate_eqn_fcn(t, c_new, bfm_phys_vars, mult, False, num_boxes, True)

    percent_error_matrix = numpy.zeros([number_of_species, number_of_species, num_boxes])

    # for i in range(0,conc.shape[0]):
    #     c_original = conc[i,:]
    #     c_new = conc[i,:]
    for j in range(number_of_species):
        c_new[:,j] = 0.0
        
        # dc_dt_new = rate_eqn_fcn(t, c_new)
        # dc_dt_new = rate_eqn_fcn(0, bfm_phys_vars, t, c_new, mult)
        # dc_dt_new = rate_eqn_fcn(t, c_new, bfm_phys_vars, mult, num_boxes, reducing=False)
        dc_dt_new = rate_eqn_fcn(t, c_new, bfm_phys_vars, mult, False, num_boxes, True)
        new = dc_dt_new
        og = dc_dt_og
        # percent_error = calc_percent_error(dc_dt_new, dc_dt_og)
        percent_error = calc_percent_error(new, og, number_of_species, num_boxes)
        percent_error_matrix[:,j,:] = percent_error
        c_new[:,j] = c_original[:,j]

            
    percent_error = numpy.amax(percent_error_matrix,axis=2)

    # Find the maximum value along each row (used for normalization)
    # row_max = numpy.amax(percent_error_matrix, axis=1)
    row_max = numpy.amax(percent_error, axis=1)

    # Normalize percent_error_matrix by row max which is new_dic_matrix
    dic_matrix = numpy.zeros([number_of_species, number_of_species])
    for i in range(number_of_species):
        if row_max[i] == 0:
            dic_matrix[i,:] = 0.0
        else:
            # dic_matrix[i,:] = percent_error_matrix[i,:]/row_max[i]
            dic_matrix[i,:] = percent_error[i,:]/row_max[i]

    # Set diagonals to zero to avoid self_directing graph edges
    numpy.fill_diagonal(dic_matrix, 0.0) 

    return dic_matrix


# def calc_percent_error(new_array, old_array):
#     """ Calculates the percent error between two matricies.
#     This is used for the calculating the New Method's direct interaction coeffs
#     """
#     percent_error = []
#     if len(new_array)==len(old_array):
#         for i in range(len(new_array)):
#             if new_array[i] == old_array[i]:
#                 percent_error.append(0.0)
#             else:
#                 percent_error.append(abs(new_array[i] - old_array[i])/abs(old_array[i])*100)

#     return percent_error

def calc_percent_error(new_matrix,old_matrix,num_species,num_boxes):
    """ Calculates the percent error between two matricies.
    This is used for the calculating the New Method's direct interaction coeffs
    """
    # percent_error = numpy.zeros((num_species,num_species,depth))
    percent_error = numpy.zeros_like(new_matrix)
    if new_matrix.shape == old_matrix.shape:
        for i in range(0,num_species):
            for k in range(0,num_boxes):
                if new_matrix[i,k] != old_matrix[i,k]:
                    percent_error[i,k] = 100*(abs(new_matrix[i,k] - old_matrix[i,k])/abs(old_matrix[i,k]))
                if numpy.isnan(percent_error[i,k]):
                    percent_error[i,k] = 0


    return percent_error


def reduce_modified_DRGEP(error_fcn, last_error, species_names, species_safe, threshold, overall_interaction_coeffs, last_num_species, solution_full_model, solution_reduced_model, t_span, c0, bfm_phys_vars):
# def reduce_modified_DRGEP(error_fcn, species_names, species_safe, threshold, overall_interaction_coeffs, solution_full_model, t_span, c0):
    """ calculates the number of species and error for a given threshold value
    """

    # Find species to remove using cutoff threshold
    species_removed = {}
    for species in overall_interaction_coeffs:
        if overall_interaction_coeffs[species] < threshold and species not in species_safe:
            species_removed[species] = numpy.nan
    # Assign index to dictionary values for the species_removed
    for index, species in enumerate(species_names):
        if species in species_removed:
            species_removed[species] = index

    # Count how many species are remaining
    num_species = len(species_names) - len(species_removed)

    if num_species == last_num_species:
        error = last_error
        solution_reduced_model = solution_reduced_model
    else:
        # For the reduced model, call function to get error
        error, solution_reduced_model = error_fcn(species_removed, solution_full_model, t_span, c0, bfm_phys_vars)
        # error, solution_reduced_model = error_fcn(species_removed, solution_full_model, t_span, c0)

    return error, num_species, species_removed, solution_reduced_model


def run_modified_DRGEP(overall_interaction_coeffs, error_limit, error_fcn, species_names, species_targets, species_safe, rate_eqn_fcn, t_span, c0, bfm_phys_vars):
# def run_modified_DRGEP(overall_interaction_coeffs, error_limit, error_fcn, species_names, species_targets, species_safe, rate_eqn_fcn, t_span, c0):
    
    """ Iterates through different threshold values to find a reduced model that meets error criteria
    """

    assert species_targets, 'Need to specify at least one target species.'

    # begin reduction iterations
    logging.info('Beginning reduction loop')
    logging.info(45 * '-')
    logging.info('Threshold | Number of species | Max error (%)')

    # make lists to store data
    threshold_data = []
    num_species_data = []
    error_data = []
    species_removed_data = []
    solution_reduced_models = []
    output_data = {}

    # start with detailed (starting) model
    # solution_full_model = solve_ivp(rate_eqn_fcn, t_span, c0, method='RK23')
    num_species = c0.shape[1]
    num_boxes = c0.shape[0]
    mult = numpy.ones(num_species)
    # solution_full_model = solve_ivp(rate_eqn_fcn, t_span, c0.ravel(), method='RK23',args=(bfm_phys_vars,mult,True,num_boxes,True))    
    solution_full_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, mult, True, num_boxes, True), t_span, c0.ravel(), method='RK23')#,args=(bfm_phys_vars,multiplier,True,num_boxes,True))

    first = True
    error_current = 0.0
    # threshold = 1e-5 #1e-5 #4e-4 #1e-6 #0.001
    threshold = 5e-5
    threshold_increment = 0.0001 #0.01
    threshold_multiplier = 2
    while error_current <= error_limit:
        if first: solution_reduced_model = solution_full_model
        error_current, num_species, species_removed, solution_reduced_model = reduce_modified_DRGEP(error_fcn, error_current, species_names, species_safe, threshold, overall_interaction_coeffs, num_species, solution_full_model, solution_reduced_model, t_span, c0, bfm_phys_vars)
        # error_current, num_species, species_removed, solution_reduced_model = reduce_modified_DRGEP(error_fcn, species_names, species_safe, threshold, overall_interaction_coeffs, solution_full_model,  t_span, c0)

        # reduce threshold if past error limit on first iteration
        if first and error_current > error_limit:
            error_current = 0.0
            threshold /= 10
            threshold_increment /= 10
            if threshold <= 1e-10:
                raise SystemExit(
                    'Threshold value dropped below 1e-10 without producing viable reduced model'
                    )
            logging.info('Threshold value too high, reducing by factor of 10')
            continue

        logging.info(f'{threshold:^9.2e} | {num_species:^17} | {error_current:^.5f}')

        # store data
        threshold_data.append(threshold)
        num_species_data.append(num_species)
        error_data.append(error_current)
        species_removed_data.append(species_removed)
        solution_reduced_models.append(solution_reduced_model)

        threshold += threshold_increment
        # threshold *= threshold_multiplier
        first = False

        # Stop reduction process if num species reaches one
        if num_species <= 1:
            break

        # Stop interating if threshold exceeds 1
        # if threshold >= 1:
        if threshold >= 3e-2:
            break

    if error_current > error_limit:
        threshold -= (2 * threshold_increment)
        error_current, num_species, species_removed, solution_reduced_model = reduce_modified_DRGEP(error_fcn, error_current, species_names, species_safe, threshold, overall_interaction_coeffs, num_species, solution_full_model, solution_reduced_model, t_span, c0, bfm_phys_vars)
        
    logging.info(45 * '-')
    logging.info('New method reduction complete.')

    # Store all output data to a dictionary
    output_data['threshold_data'] = threshold_data
    output_data['num_species_data'] = num_species_data
    output_data['error_data'] = error_data
    output_data['species_removed_data'] = species_removed_data
    output_data['solution_full_model'] = solution_full_model
    output_data['solution_reduced_models'] = solution_reduced_models

    return output_data


def modified_DRGEP(conc,bfm_phys_vars):

    # Information for reduction
    # User needs to update these values
    # scenario = 'chl_2'
    # species_targets = ['P2l', 'P3l', 'P4l']
    # species_safe = ['N4n', 'P1l', 'P2l', 'P3l', 'P4l']

    scenario = 'o2o_17'
    species_targets = ['O2o']
    # species_safe = []
    species_safe = ['N4n','P1l','P2l','P3l','P4l']
    error_limit = 5
    error_fcn = error_fcns.calc_error_o2o_17
    # rate_eqn_fcn = bfm50_rate_eqns
    rate_eqn_fcn = reduced_bfm56_rate_eqns
    # rate_eqn_fcn_1d = bfm50_rate_eqns

    # Time span for integration
    t_span = [0, 86400*365*10]
    # t_span = [0, 86400*365]

    # Time at which the DIC values are obtained
    t = 86400*365

    # Log input data to file
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename='output.log', level=logging.INFO)
    logging.info('Scenario: {}'.format(scenario))
    logging.info('Error function: {}'.format(error_fcn))
    logging.info('Error limit: {}'.format(error_limit))
    logging.info('Target species: {}'.format(species_targets))
    logging.info('Retained species: {}'.format(species_safe))

    # Get direct interaction coefficients
    # dc_dt_og_0 = rate_eqn_fcn(t, conc)
    # dic_matrix = calc_modified_DRGEP_dic(rate_eqn_fcn, conc, t)

    # interaction_coeffs = numpy.ones_like(conc)
    # for i in range(0,interaction_coeffs.shape[0]):
    #     dic_matrix = calc_modified_DRGEP_dic(rate_eqn_fcn, conc[i,:], t)
    #     coeffs = get_importance_coeffs(species_names, species_targets, [dic_matrix])
        
    #     interaction_coeffs[i,:] = extract_dict(coeffs)
    #     # for j in coeffs.items():
    #     #     interaction_coeffs[i,j] = coeffs[j].value()

    dic_matrix = calc_modified_DRGEP_dic(rate_eqn_fcn, conc, t, bfm_phys_vars)
    overall_interaction_coeffs = get_importance_coeffs(species_names, species_targets, [dic_matrix])

    # # Get overall interaction coefficients
    # max_coeffs = interaction_coeffs.max(axis=0)
    # min_coeffs = interaction_coeffs.min(axis=0)

    # overall_interaction_coeffs = fill_dict(coeffs,max_coeffs,species_names)
    # # overall_interaction_coeffs = get_importance_coeffs(species_names, species_targets, [dic_matrix])

    # Make dictionary of target species and their index values
    target_species = {}
    for index, species in enumerate(species_names):
        if species in species_targets:
            target_species[species] = index

    # Slice concentration matrix
    # spacing = 15
    spacing = 25
    num_boxes = int(numpy.floor(conc.shape[0]/spacing)+1)
    conc_sliced = numpy.zeros((num_boxes,conc.shape[1]))
    phys_vars_sliced = BfmPhysicalVariableData(num_boxes)
    phys_vars_sliced.pH = bfm_phys_vars.pH
    for i in range(0,num_boxes-1):
        # if i > 0: slice = spacing*i-1
        # else: slice=0
        slice = spacing*i
        conc_sliced[i,:] = conc[slice,:]
        phys_vars_sliced.depth[i] = bfm_phys_vars.depth[slice]
    conc_sliced[-1] = conc[-1,:]
    phys_vars_sliced.depth[-1] = bfm_phys_vars.depth[-1]
    mult = numpy.ones(50)
    # solution_full_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, mult, True, 150, True), t_span, conc.ravel(), method='RK23')#,args=(bfm_phys_vars,multiplier,True,num_boxes,True))

    # Run new method
    # reduction_data = run_modified_DRGEP(overall_interaction_coeffs, error_limit, error_fcn, species_names, species_targets, species_safe, rate_eqn_fcn_1d, t_span, conc[0,:])
    reduction_data = run_modified_DRGEP(overall_interaction_coeffs, error_limit, error_fcn, species_names, species_targets, species_safe, rate_eqn_fcn, t_span, conc_sliced, phys_vars_sliced)

    return reduction_data, error_limit
