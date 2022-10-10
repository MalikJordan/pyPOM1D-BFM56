import numpy as np
from pom.constants import num_boxes

def settling_dynamics(d3state, detritus_sedimentation, phyto_sedimentation, pel_bottom, group, constant_parameters, phyto_parameters):

    if group == 1: # P1: Diatoms terms
        pc = d3state[10,:]
        pn = d3state[11,:]
        pp = d3state[12,:]
        pl = d3state[13,:]
        ps = d3state[14,:]
    elif group == 2: # P2: Flagellates terms
        pc = d3state[15,:]
        pn = d3state[16,:]
        pp = d3state[17,:]
        pl = d3state[18,:]
    elif group == 3: # P3: PicoPhytoplankton terms
        pc = d3state[19,:]
        pn = d3state[20,:]
        pp = d3state[21,:]
        pl = d3state[22,:]
    elif group == 4: # P4: Large Phytoplankton terms
        pc = d3state[23,:]
        pn = d3state[24,:]
        pp = d3state[25,:]
        pl = d3state[26,:]
    
    # Set the phytoplankton sinking velocity at the water sediment interface (from Settling.F90)      
    if (phyto_parameters["p_res"] > 0):    
        if (phyto_sedimentation[num_boxes-1,group] > phyto_parameters["p_burvel_PI"]):
            sedi = phyto_parameters["p_burvel_PI"]*np.ones(num_boxes)
        else:
            sedi = phyto_sedimentation[num_boxes-1,group]*np.ones(num_boxes)    # num_boxes-1 for bottom sedimentation rate

    ruQIc = np.zeros(num_boxes)
    ruQ1c = np.zeros(num_boxes)
    ruQ6c = np.zeros(num_boxes)
    for i in range(0,num_boxes):
        if (sedi[i] > 0.):
            # Phytoplankton carbon
            ruQIc[i] = sedi[i]*pc[i]
            ruQ1c[i] = constant_parameters["p_pe_R1c"]*ruQIc[i]
            ruQ6c[i] = ruQIc[i] - ruQ1c[i]
            pel_bottom[i,group] = -ruQIc[i]
            