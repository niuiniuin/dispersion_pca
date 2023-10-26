import bilby
import numpy as np
import scipy.integrate as integrate
from matplotlib import pyplot as plt

from astropy.cosmology import Planck18
import astropy.cosmology as cosmos
from astropy import units as u
from astropy import constants as const
import lal


alpha_00 = 0.0
alpha_05 = 0.5
alpha_10 = 1.0
alpha_15 = 1.5
alpha_25 = 2.5
alpha_30 = 3.0
alpha_35 = 3.5
alpha_40 = 4.0

mass_factor = 1e3
# mass_factor = 1e4
# conver_fac = (const.c.to('km/s')/Planck18.H0).to_value(u.Mpc)
conver_fac = (const.c.to('km/s')/Planck18.H0).to_value(u.Gpc)
OmegaM = Planck18.Om0
OmegaL = Planck18.Ode0

# (1+redshift)**(1-alpha_00)
# can be cancled

@np.vectorize
def D_alpha(redshift):
    '''
    returned value is in the unit of Gpc (note the unit used in conver_fac)
    '''
    D_00 = (1+redshift)**(1-alpha_00) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_00-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_05 = (1+redshift)**(1-alpha_05) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_05-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_10 = (1+redshift)**(1-alpha_10) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_10-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_15 = (1+redshift)**(1-alpha_15) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_15-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_25 = (1+redshift)**(1-alpha_25) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_25-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_30 = (1+redshift)**(1-alpha_30) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_30-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_35 = (1+redshift)**(1-alpha_35) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_35-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_40 = (1+redshift)**(1-alpha_40) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_40-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    return D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40  # in the unit of Gpc (note the unit used in conver_fac)


def delta_phase_nonGR_dispersion_alpha_all(f_array, total_mass, chirp_mass, redshift, D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40,
                                           dphi_Aalpha_00, dphi_Aalpha_05, dphi_Aalpha_10, dphi_Aalpha_15, dphi_Aalpha_25, dphi_Aalpha_30, dphi_Aalpha_35, dphi_Aalpha_40):
    
    M_sec = total_mass * lal.MTSUN_SI
    Mc_sec = chirp_mass * lal.MTSUN_SI

    phase_nonGR_00 = dphi_Aalpha_00 * (M_sec*f_array)**(alpha_00-1) * np.pi*(1+redshift)**(alpha_00-1)/(alpha_00-1) * (total_mass/mass_factor)**(1-alpha_00) * D_00
    phase_nonGR_05 = dphi_Aalpha_05 * (M_sec*f_array)**(alpha_05-1) * np.pi*(1+redshift)**(alpha_05-1)/(alpha_05-1) * (total_mass/mass_factor)**(1-alpha_05) * D_05
    phase_nonGR_10 = dphi_Aalpha_10 * np.pi * D_10 * np.log(np.pi*Mc_sec*f_array)
    phase_nonGR_15 = dphi_Aalpha_15 * (M_sec*f_array)**(alpha_15-1) * np.pi*(1+redshift)**(alpha_15-1)/(alpha_15-1) * (total_mass/mass_factor)**(1-alpha_15) * D_15
    phase_nonGR_25 = dphi_Aalpha_25 * (M_sec*f_array)**(alpha_25-1) * np.pi*(1+redshift)**(alpha_25-1)/(alpha_25-1) * (total_mass/mass_factor)**(1-alpha_25) * D_25
    phase_nonGR_30 = dphi_Aalpha_30 * (M_sec*f_array)**(alpha_30-1) * np.pi*(1+redshift)**(alpha_30-1)/(alpha_30-1) * (total_mass/mass_factor)**(1-alpha_30) * D_30
    phase_nonGR_35 = dphi_Aalpha_35 * (M_sec*f_array)**(alpha_35-1) * np.pi*(1+redshift)**(alpha_35-1)/(alpha_35-1) * (total_mass/mass_factor)**(1-alpha_35) * D_35
    phase_nonGR_40 = dphi_Aalpha_40 * (M_sec*f_array)**(alpha_40-1) * np.pi*(1+redshift)**(alpha_40-1)/(alpha_40-1) * (total_mass/mass_factor)**(1-alpha_40) * D_40
    
    phase_nonGR_tol = phase_nonGR_00 + phase_nonGR_05 + phase_nonGR_10 + phase_nonGR_15 + phase_nonGR_25 + phase_nonGR_30 + phase_nonGR_35 + phase_nonGR_40

    return phase_nonGR_tol

z_interp = np.linspace(0.01, 1.5, 10000)
dL_interp = Planck18.luminosity_distance(z_interp).value
D_00_interp, D_05_interp, D_10_interp, D_15_interp, D_25_interp, D_30_interp, D_35_interp, D_40_interp = D_alpha(z_interp)
def interp_z_D_alpha(dL_in):
    z = np.interp(dL_in, dL_interp, z_interp)

    D_00 = np.interp(dL_in, dL_interp, D_00_interp) 
    D_05 = np.interp(dL_in, dL_interp, D_05_interp) 
    D_10 = np.interp(dL_in, dL_interp, D_10_interp) 
    D_15 = np.interp(dL_in, dL_interp, D_15_interp) 
    D_25 = np.interp(dL_in, dL_interp, D_25_interp) 
    D_30 = np.interp(dL_in, dL_interp, D_30_interp) 
    D_35 = np.interp(dL_in, dL_interp, D_35_interp) 
    D_40 = np.interp(dL_in, dL_interp, D_40_interp)

    return z, D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40

# test plot
rng = np.random.default_rng()
z_test_in = rng.uniform(0.01, 1.5, 500)
dL_test = Planck18.luminosity_distance(z_test_in).value
z_test_out, D_00_test, D_05_test, D_10_test, D_15_test, D_25_test, D_30_test, D_35_test, D_40_test = interp_z_D_alpha(dL_test)

fig, ax = plt.subplots()
ax.plot(z_interp, D_00_interp, color='C1', label='D_00_interp')
ax.plot(z_interp, D_05_interp, color='C2', label='D_05_interp')
ax.plot(z_interp, D_10_interp, color='C3', label='D_10_interp')
ax.plot(z_interp, D_15_interp, color='C4', label='D_15_interp')
ax.plot(z_interp, D_25_interp, color='C5', label='D_25_interp')
ax.plot(z_interp, D_30_interp, color='C6', label='D_30_interp')
ax.plot(z_interp, D_35_interp, color='C7', label='D_35_interp')
ax.plot(z_interp, D_40_interp, color='C8', label='D_40_interp')
ax.scatter(z_test_out, D_00_test, s=12, marker='x', color='C1', label='D_00_test')
ax.scatter(z_test_out, D_05_test, s=12, marker='x', color='C2', label='D_05_test')
ax.scatter(z_test_out, D_10_test, s=12, marker='x', color='C3', label='D_10_test')
ax.scatter(z_test_out, D_15_test, s=12, marker='x', color='C4', label='D_15_test')
ax.scatter(z_test_out, D_25_test, s=12, marker='x', color='C5', label='D_25_test')
ax.scatter(z_test_out, D_30_test, s=12, marker='x', color='C6', label='D_30_test')
ax.scatter(z_test_out, D_35_test, s=12, marker='x', color='C7', label='D_35_test')
ax.scatter(z_test_out, D_40_test, s=12, marker='x', color='C8', label='D_40_test')
ax.set_xlabel('redshift')
ax.set_ylabel('Gpc')
ax.legend()
fig.savefig('interp_Dalpha_test.png')

fig, ax = plt.subplots()
ax.plot(z_interp, dL_interp, color='C0', label='dL_interp')
ax.scatter(z_test_out, dL_test, s=12, marker='x', color='C0', label='dL_test')
ax.set_xlabel('redshift')
ax.set_ylabel('Mpc')
ax.legend()
fig.savefig('interp_DL_test.png')

################################################################################
def waveform_alpha_all(frequency_array, 
                       mass_1, mass_2, 
                       luminosity_distance, theta_jn, phase, 
                       a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, 
                       A_alpha_00, A_alpha_05, A_alpha_10, A_alpha_15, A_alpha_25, A_alpha_30, A_alpha_35, A_alpha_40, 
                       **kwargs):

    total_mass = mass_1 + mass_2
    eta = mass_1 * mass_2 / (total_mass**2)
    chirp_mass = total_mass * eta**(3/5)
    z, D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40 = interp_z_D_alpha(luminosity_distance)

    waveform_GR = bilby.gw.source.lal_binary_black_hole(
        frequency_array=frequency_array, 
        mass_1=mass_1, mass_2=mass_2, 
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase, 
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12, a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl, 
        **kwargs)
    h_plus  = waveform_GR['plus']
    h_cross = waveform_GR['cross']

    minimum_frequency = kwargs['minimum_frequency']
    maximum_frequency = frequency_array[-1]
    frequency_bounds = ((frequency_array >= minimum_frequency) * 
                        (frequency_array <= maximum_frequency))
    
    f_array = frequency_array[frequency_bounds]
    
    nonGR_tol = delta_phase_nonGR_dispersion_alpha_all(f_array, total_mass, chirp_mass, 
                                                       z, D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40,
                                                       A_alpha_00, A_alpha_05, A_alpha_10, A_alpha_15, A_alpha_25, A_alpha_30, A_alpha_35, A_alpha_40)
    pad_lower = (frequency_array < minimum_frequency).sum()
    pad_upper = (frequency_array > maximum_frequency).sum()
    nonGR_tol = np.pad(nonGR_tol, (pad_lower, pad_upper), constant_values=(0,))

    h_plus_nonGR =  h_plus *  np.exp(1j*nonGR_tol)
    h_cross_nonGR = h_cross * np.exp(1j*nonGR_tol)

    # print(frequency_array)
    # print(f_array)
    # print(frequency_bounds)
    # print(pad_lower)
    # print(pad_upper)
    # print(len(frequency_array))
    # print(len(f_array))
    # print(len(nonGR_tol))
    # print(np.exp(1j*nonGR_tol))

    return dict(plus=h_plus_nonGR, cross=h_cross_nonGR)

