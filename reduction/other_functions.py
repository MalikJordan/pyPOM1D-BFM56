import numpy
from scipy.integrate import solve_ivp
from reduction.PythonBFM.BFM50_rate_eqns import bfm50_rate_eqns
from bfm.bfm56.BFM56_rate_eqns import reduced_bfm56_rate_eqns
import time

def solution_time(conc, bfm_phys_vars, mult,):
    """ calculates the computation time for 0D integration
    """

    # Time span for integration
    t_span = [0, 86400*365*10]
    
    # Record integration start time
    start = time.process_time()

    # Integrate model
    # solution = solve_ivp(bfm50_rate_eqns, t_span, conc, method='RK23')
    solution_full_model = solve_ivp(lambda time, conc: reduced_bfm56_rate_eqns(time, conc, bfm_phys_vars, mult, True, 150, True), t_span, conc.ravel(), method='RK23')#,args=(bfm_phys_vars,multiplier,True,num_boxes,True))


    # Record integration end time
    end = time.process_time()

    # Time of integration
    integration_time = end - start

    return integration_time


def extract_dict(dictionary):

    data = list(dictionary.values())
    values = numpy.array(data)

    return values


def fill_dict(dictionary,array,species_names):

    for i, sp in enumerate(species_names):
        dictionary[sp] = array[i]
    
    return dictionary
