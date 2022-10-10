from inputs import params_POMBFM

current_path = '/Users/malikjordan/Documents/GitHub/pyPOM1D-BFM56'
seconds_per_day = 86400.
earth_angular_velocity = 7.29E-5                                                        # OMEGA
vertical_layers = 151
num_boxes = vertical_layers - 1
DAYI = 1. / seconds_per_day
water_specific_heat_times_density = 4.187E6
twice_the_timestep = 2. * params_POMBFM.dti
