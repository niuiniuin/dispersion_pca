import numpy as np
import scipy.integrate as integrate

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

conver_fac = (const.c.to('km/s')/Planck18.H0).to_value(u.Gpc)
OmegaM = Planck18.Om0
OmegaL = Planck18.Ode0

@np.vectorize
def D_alpha(redshift):
    '''
    returned value is in the unit of Mpc
    '''
    D_00 = (1+redshift)**(1-alpha_00) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_00-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_05 = (1+redshift)**(1-alpha_05) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_05-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_10 = (1+redshift)**(1-alpha_10) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_10-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_15 = (1+redshift)**(1-alpha_15) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_15-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_25 = (1+redshift)**(1-alpha_25) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_25-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_30 = (1+redshift)**(1-alpha_30) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_30-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_35 = (1+redshift)**(1-alpha_35) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_35-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    D_40 = (1+redshift)**(1-alpha_40) * conver_fac * integrate.quad(lambda z: (1+z)**(alpha_40-2) / np.sqrt(OmegaM*(1+z)**3+OmegaL), 0, redshift)[0]
    return D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40  # in the unit of Mpc

z_interp = np.linspace(0.01, 0.17, 10000)
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



# print(cosmos.z_at_value(Planck18.luminosity_distance, 50*u.Mpc))
# print(cosmos.z_at_value(Planck18.luminosity_distance, 800*u.Mpc))



# params_nonGR = dict(A_alpha_00=0.0, A_alpha_05=0.0, A_alpha_10=0.0, A_alpha_15=0.0, A_alpha_25=0.0, A_alpha_30=0.0, A_alpha_35=0.0, A_alpha_40=0.0)
# waveform_arguments = dict(waveform_approximant='IMRPhenomXAS', 
#                           reference_frequency=20., 
#                           minimum_frequency=20.,
#                           maximum_frequency=1024.,
#                           catch_waveform_errors=True)

# mass_1 = 30
# mass_2 = 20
# luminosity_distance = 200
# theta_jn = 0.0
# phase = 0.0
# a_1 = 0.0
# tilt_1 = 0.0
# phi_12 = 0.0
# a_2 = 0.0
# tilt_2 = 0.0
# phi_jl = 0.0

# frequency_array = np.arange(0, 1024, 1/4)
# h = waveform_alpha_all(frequency_array, 
#                         mass_1, mass_2, 
#                         luminosity_distance, theta_jn, phase, 
#                         a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, 
#                         params_nonGR['A_alpha_00'], params_nonGR['A_alpha_05'], params_nonGR['A_alpha_10'], params_nonGR['A_alpha_15'], params_nonGR['A_alpha_25'], params_nonGR['A_alpha_30'], params_nonGR['A_alpha_35'], params_nonGR['A_alpha_40'], 
#                         **waveform_arguments)

# plt.figure()
# para = ["A_alpha_00","A_alpha_05","A_alpha_10","A_alpha_15","A_alpha_25","A_alpha_30","A_alpha_35","A_alpha_40"]
# for i in range(len(para)):
#    params_temp = params_nonGR.copy()
#    params_temp[para[i]] = 1.0
#    h = waveform_alpha_all(frequency_array, 
#                                                      mass_1, mass_2, 
#                                                      luminosity_distance, theta_jn, phase, 
#                                                      a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, 
#                                                      params_temp['A_alpha_00'], params_temp['A_alpha_05'], params_temp['A_alpha_10'], params_temp['A_alpha_15'], params_temp['A_alpha_25'], params_temp['A_alpha_30'], params_temp['A_alpha_35'], params_temp['A_alpha_40'], 
#                                                      **waveform_arguments)
#    plt.semilogx(frequency_array, h['plus'].real, label=f'{para[i]}=1.0')

# h = waveform_alpha_all(frequency_array, 
#                                                   mass_1, mass_2, 
#                                                   luminosity_distance, theta_jn, phase, 
#                                                   a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, 
#                                                   params_nonGR['A_alpha_00'], params_nonGR['A_alpha_05'], params_nonGR['A_alpha_10'], params_nonGR['A_alpha_15'], params_nonGR['A_alpha_25'], params_nonGR['A_alpha_30'], params_nonGR['A_alpha_35'], params_nonGR['A_alpha_40'], 
#                                                   **waveform_arguments)
# plt.semilogx(frequency_array, h['plus'].real, label=f'GR')

# plt.legend()
# plt.savefig('real_all_Aalpha.png')



# print(interp_z_D_alpha(200))
# mass_1 = 5e5
# mass_2 = 3e5
# luminosity_distance = 3e3

mass_1 = 70
mass_2 = 50
luminosity_distance = 1000

total_mass = mass_1 + mass_2
eta = mass_1 * mass_2 / (total_mass**2)
chirp_mass = total_mass * eta**(3/5)
z, D_00, D_05, D_10, D_15, D_25, D_30, D_35, D_40 = interp_z_D_alpha(luminosity_distance)

dphi_Aalpha_00 = 1.0
dphi_Aalpha_05 = 1.0
dphi_Aalpha_10 = 1.0
dphi_Aalpha_15 = 1.0
dphi_Aalpha_25 = 1.0
dphi_Aalpha_30 = 1.0
dphi_Aalpha_35 = 1.0
dphi_Aalpha_40 = 1.0

M_sec = total_mass * lal.MTSUN_SI
Mc_sec = chirp_mass * lal.MTSUN_SI

alpha_00 = 0.0
alpha_05 = 0.5
alpha_10 = 1.0
alpha_15 = 1.5
alpha_25 = 2.5
alpha_30 = 3.0
alpha_35 = 3.5
alpha_40 = 4.0

mass_factor = 1e3

phase_nonGR_00 = dphi_Aalpha_00 * np.pi*(1+z)**(alpha_00-1)/(alpha_00-1) * (total_mass/mass_factor)**(1-alpha_00) * D_00
phase_nonGR_05 = dphi_Aalpha_05 * np.pi*(1+z)**(alpha_05-1)/(alpha_05-1) * (total_mass/mass_factor)**(1-alpha_05) * D_05
phase_nonGR_10 = dphi_Aalpha_10 * np.pi*D_10
phase_nonGR_15 = dphi_Aalpha_15 * np.pi*(1+z)**(alpha_15-1)/(alpha_15-1) * (total_mass/mass_factor)**(1-alpha_15) * D_15
phase_nonGR_25 = dphi_Aalpha_25 * np.pi*(1+z)**(alpha_25-1)/(alpha_25-1) * (total_mass/mass_factor)**(1-alpha_25) * D_25
phase_nonGR_30 = dphi_Aalpha_30 * np.pi*(1+z)**(alpha_30-1)/(alpha_30-1) * (total_mass/mass_factor)**(1-alpha_30) * D_30
phase_nonGR_35 = dphi_Aalpha_35 * np.pi*(1+z)**(alpha_35-1)/(alpha_35-1) * (total_mass/mass_factor)**(1-alpha_35) * D_35
phase_nonGR_40 = dphi_Aalpha_40 * np.pi*(1+z)**(alpha_40-1)/(alpha_40-1) * (total_mass/mass_factor)**(1-alpha_40) * D_40
print('phase_nonGR_00: ', phase_nonGR_00)
print('phase_nonGR_05: ', phase_nonGR_05)
print('phase_nonGR_10: ', phase_nonGR_10)
print('phase_nonGR_15: ', phase_nonGR_15)
print('phase_nonGR_25: ', phase_nonGR_25)
print('phase_nonGR_30: ', phase_nonGR_30)
print('phase_nonGR_35: ', phase_nonGR_35)
print('phase_nonGR_40: ', phase_nonGR_40)


f_array = 100
phase_nonGR_00 = dphi_Aalpha_00 * (M_sec*f_array)**(alpha_00-1) * np.pi*(1+z)**(alpha_00-1)/(alpha_00-1) * (total_mass/mass_factor)**(1-alpha_00) * D_00
phase_nonGR_05 = dphi_Aalpha_05 * (M_sec*f_array)**(alpha_05-1) * np.pi*(1+z)**(alpha_05-1)/(alpha_05-1) * (total_mass/mass_factor)**(1-alpha_05) * D_05
phase_nonGR_10 = dphi_Aalpha_10 * np.pi * D_10 * np.log(np.pi*Mc_sec*f_array)
phase_nonGR_15 = dphi_Aalpha_15 * (M_sec*f_array)**(alpha_15-1) * np.pi*(1+z)**(alpha_15-1)/(alpha_15-1) * (total_mass/mass_factor)**(1-alpha_15) * D_15
phase_nonGR_25 = dphi_Aalpha_25 * (M_sec*f_array)**(alpha_25-1) * np.pi*(1+z)**(alpha_25-1)/(alpha_25-1) * (total_mass/mass_factor)**(1-alpha_25) * D_25
phase_nonGR_30 = dphi_Aalpha_30 * (M_sec*f_array)**(alpha_30-1) * np.pi*(1+z)**(alpha_30-1)/(alpha_30-1) * (total_mass/mass_factor)**(1-alpha_30) * D_30
phase_nonGR_35 = dphi_Aalpha_35 * (M_sec*f_array)**(alpha_35-1) * np.pi*(1+z)**(alpha_35-1)/(alpha_35-1) * (total_mass/mass_factor)**(1-alpha_35) * D_35
phase_nonGR_40 = dphi_Aalpha_40 * (M_sec*f_array)**(alpha_40-1) * np.pi*(1+z)**(alpha_40-1)/(alpha_40-1) * (total_mass/mass_factor)**(1-alpha_40) * D_40
print('(with Mf)phase_nonGR_00: ', phase_nonGR_00)
print('(with Mf)phase_nonGR_05: ', phase_nonGR_05)
print('(with Mf)phase_nonGR_10: ', phase_nonGR_10)
print('(with Mf)phase_nonGR_15: ', phase_nonGR_15)
print('(with Mf)phase_nonGR_25: ', phase_nonGR_25)
print('(with Mf)phase_nonGR_30: ', phase_nonGR_30)
print('(with Mf)phase_nonGR_35: ', phase_nonGR_35)
print('(with Mf)phase_nonGR_40: ', phase_nonGR_40)