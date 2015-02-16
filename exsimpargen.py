
# coding: utf-8

# In[1]:

#!/usr/bin/env python
"""exsimpargen.py: utility for generating input parameter files for exsim simulation software"""

__author__ = "Massimiliano Pittore"
__copyright__ = "Copyright 2015, GFZ, Centre for Early Warning"
__license__ = "GPL"
__maintainer__ = "Massimiliano Pittore"
__email__ = "pittore@gfz-potsdam.de"


# In[5]:

from string import Template
import numpy as np
import random
from random import choice
import scipy.stats as stats


# In[12]:

#init parameter file text with placeholders for variables
#NOTE: don´t modify it
par_text=Template('!Input file for program EXSIM12\n!Title\n  $title\n!Write acc, psa, husid files for each site?\n $write_files\n!MW, Stress, flag (0=fmax; 1=kappa), fmax or kappa\n  $mw $stress  $kappa_flag  $fmax_kappa\n!lat and lon of upper edge of fault\n  $upedge_lat $upedge_lon\n!strike,dip, depth of fault\n  $strike $dip $depth_of_fault\n!fault type (S=strikeslip; R=reverse; N=normal; U=undifferentiated)\n! (Only used if Wells and Coppersmith is used to obtain FL and FW).\n  $fault_style\n!fault length and width, dl, dw, stress_ref\n!Note: Force program to use Wells and Coppersmith (WC) for FL and/or FW if\n! either entry = 0.0.\n! dl and dw are the subsource length and width\n! stress_ref is a reference to allow scaling of WC size as per Atkinson&Boore(2006BSSA)\n! If Wells and Coppersmith are used to obtain FL and/or FW, the WC values are\n! modified to account for the scaling implied by differences in the stress\n! specified above and a stress that is assumed to be valid for the generic WC\n! relations; this stress is stress_ref. The value of 70 bars is an educated\n! guess for stress_ref, but it is not based on a quantitative analysis.\n! The WC values of FL and/or FW are multiplied by the factor\n! (stress_ref/stress)^(1/3).\n! Note that four entries on the following line are needed as placeholders,\n! even if not used)\n  $fault_length $fault_width $fault_sub_length $fault_sub_width $stress_ref !fault length and width, dl, dw, stress_ref\n!vrup/beta\n  $v_rup\n!hypo location in along fault and down dip distance from the fault\n!reference point (an upper corner)(-1.0, -1.0 for a random location);\n!number of iterations over hypocenter (need an entry, but only used if\n!either of the first two values are -1.0, indicating a random location)\n  $hypo_loc_along  $hypo_loc_down  $hypo_iter\n!Enter type of risetime (1=original, 2=1/f0)\n $risetime_type\n!tpadl, tpadt, delta t (length of 0pads at front and back of time series, timestep)\n $tpad_l $tpad_t $delta_t\n!beta , rho\n  $beta $rho\n!Geometric spreading: this example is for bilinear with transition at 40km\n! r_ref, nseg (hinged line segments), (rlow(i), slope)\n! (Usually set r_ref = 1.0 km)\n    1.0\n    2\n      1.0 -1.0\n     40.0 -0.5\n!Quality factor: Qmin, Q0, and eta, Q=max(Qmin, Q0*F**eta)\n   $q_min  $q_0  $eta\n!path duration: example has duration increasing as 0.05R\n!(ndur_hinges,(rdur(i), dur(i), i = 1, ndur_hinges), durslope)\n    2\n    0.0 0.0\n   10.0 0.0\n  0.05\n!Type of window: 1 for Saragoni-Hart taper windows, 0 for tapered boxcar\n!window, epsilon, and eta values of Saragoni-Hart window\n  1    0.2    0.2\n!low-cut filter corner (Hz), nslope (0 ==> no filter)\n $lc_filter_corner $lc_filter_nslope\n! %damping of response spectra\n $damping\n!# of f and Min and Max F for response spectra\n\
  $res_spectra_num_f $res_spectra_min_f   $res_spectra_max_f\n\
!no. of frequencies for summary output (10 max):\n\
 $num_f_summary_output\n\
!frequency (-1.0, 99.0 for pgv, pga):\n\
 -1.0 99.0 0.5 5.0\n\
!Output file names stem:\n\
  $out_basename\n\
!Name of crustal amplification file:\n\
  $crust_amp_file\n\
!Name of site amplification file:\n\
  $site_amp_file\n\
!Name of empirical filter file:\n\
  $emp_amp_file\n\
!DynamicFlag (0=no; use 1 for dynamic corner freq), PulsingPercent (typical 50.)\n\
  $dyn_corner_freq   $pulsing_percent\n\
!iflagscalefactor (1=vel^2; 2=acc^2; 3=asymptotic acc^2 (dmb); typical=2)\n\
  $flagscalefactor\n\
!iflagfas_avg (1=arithmetic; 2=geometric, 3=rms: USE 3!)\n\
  $flagfas_avg\n\
!iflagpsa_avg (1=arithmetic; 2=geometric: USE 2!, 3=rms)\n\
  $flagpsa_avg\n\
!deterministic flag,gama,nu,t0, impulse peak (see Motazedian and Atkinson, 2005)\n\
  $det_flag   $gamma  $nu  $t0  $impulse_peak\n\
!iseed, # of trials\n\
  $seed  $num_trials\n\
!islipweight = -1  -> unity slip for all subfaults,\n\
!islipweight =  0  -> specify slips read from text file,\n\
!islipweight =  1  -> random weights\n\
   $slipweight\n\
! Text file containing matrix of slip weights (need a placeholder\n\
! even if do not assign the slip weights\n\
  $slip_weights_file\n\
!Number of Sites, site coord flag (1=lat,long; 2=R,Az; 3=N,E)\n\
  $num_sites  $site_coord_type\n\
!If "Y" below and strike = 0.0:\n\
!  if site coord flag = 2, move origin of the radial line to the midpoint of\n\
!                         the top edge of the fault\n\
!  if site coord flag = 3 and siteLocation(1) = 0, redefine\n\
!                         siteLocation(1) = 0 to be the midpoint of the\n\
!                         top edge of the fault (so that the sites will be\n\
!                         along a line normal to the midpoint)\n\
!  if site coord flag = 3 and siteLocation(2) = 0, redefine\n\
!                         siteLocation(1) = 0 to be the far end of the fault,\n\
!                         so that the sites are along a line along the\n\
!                         strike of the fault\n\
 Y\n\
!Coordinates of each site\n\
  $site_coords\n')


# In[13]:

#set defaults as a dictionary
# NOTE: define here the default values

sim_pars=dict(title='test',
write_files='Y',
mw=7,                 # simulation variable
stress=100,           # simulation variable
kappa_flag=1,
fmax_kappa=0.01,
upedge_lat=42.7,      # simulation variable
upedge_lon=73.95,     # simulation variable
strike=0.0,          # simulation variable
dip=50.0,             # simulation variable
depth_of_fault=0.0,
fault_style='R',      # simulation variable
fault_length=0,       #force program to use Well & Coppersmith
fault_width=0,        #force program to use Well & Coppersmith
fault_sub_length=1.5, 
fault_sub_width=1.5, 
stress_ref=70,
v_rup=0.85,
hypo_loc_along=99.0,
hypo_loc_down=18.0,  
hypo_iter=1,
risetime_type=2,
tpad_l=50.0, 
tpad_t=20.0, 
delta_t=0.01,
beta=3.7, 
rho=2.8,
q_min=60,  
q_0=180,  
eta=0.45,
lc_filter_corner=0.05, 
lc_filter_nslope=8,
damping=5.0,
res_spectra_num_f=5, 
res_spectra_min_f=0.1,   
res_spectra_max_f=99.0,
num_f_summary_output=2,
out_basename='default_basename',          # simulation variable
crust_amp_file='crustal_amps.txt', # simulation variable
site_amp_file='site_amps.txt',     # simulation variable
emp_amp_file='empirical_amps.txt',       # simulation variable
dyn_corner_freq=1,   
pulsing_percent=50.0,
flagscalefactor=2,
flagfas_avg=3,
flagpsa_avg=2,
det_flag=0,   
gamma=1.0,  
nu=90.0,  
t0=4.0,  
impulse_peak=10.0,
seed=309,  
num_trials=5,                              # simulation variable
slipweight=-1,
slip_weights_file='slip_weights.txt',      # simulation variable
num_sites=1,                               # simulation variable
site_coord_type=1,                         # simulation variable
site_coords='42.854302 74.533163')         # simulation variable


# In[2]:

# define (average) values for the simulation parameters
upedge_lat=42.7
upedge_lat_sd=0.1 # lat sampled from normal
upedge_lon=73.95
upedge_lon_sd=0.1 # lon sampled from normal
mw=7
mw_sd=0.15 # mw sampled from normal
stress=100
stress_min = 50 # stress sampled from uniform
stress_max = 300 
strike=0.0    
dip=50.0
depth_of_fault=0.0


# In[ ]:

# number of samples for each simulation parameter
nmags=10
nlats=10
nlongs=10
nstress=10
#tot_instances=nstress*nlongs*nlats*nmags
#tot_instances
#around 10'000 simulations with 10 draws each


# In[ ]:

#stress_array=np.linspace(50,300,nstress)
#stress_array
#stress_array=np.random.uniform.ppf(np.linspace(0.001,0.99,nstress),50,300)


# In[ ]:

#sample values from normal distribution, with equal percentile spacing

mag_array=stats.norm.ppf(np.linspace(0.001,0.99,nmags),mw,mw_sd)

#generate stress values evenly across a range
stress_array=np.linspace(50,300,nstress)
#sample values from uniform distrib
#stress_array=np.random.uniform(50,300,nstress)

lat_array=stats.norm.ppf(np.linspace(0.001,0.99,nlats),upedge_lat,upedge_lat_sd)
lon_array=stats.norm.ppf(np.linspace(0.001,0.99,nlongs),upedge_lon,upedge_lon_sd)

#basename for the acc (and other) output files
base_accfiles='test_bishkek'
#basename for the parameter files being generated
base_initfiles='test_bishkek'

#where to put the generated par files
outpath='/home/max/Desktop/Documents/GFZ_sync/workspace/EXSIM_inputgen/parfiles/'

#randomly subsample from the parameters´ arrays
for m in np.random.choice(mag_array,size=5,replace=False):
    for s in np.random.choice(stress_array,size=5,replace=False):
        for lat in np.random.choice(lat_array,size=5,replace=False):
            for lon in np.random.choice(lon_array,size=5,replace=False):
                sim_pars['upedge_lat']='{0:.3f}'.format(lat)
                sim_pars['upedge_lon']='{0:.3f}'.format(lon)
                sim_pars['mw']='{0:.1f}'.format(m)
                sim_pars['stress']='{0:.1f}'.format(s)
                sim_pars['out_basename']=base_accfiles+'_{0:.1f}'.format(m)+'_{0:.1f}_'.format(s)
                #should the seed be changed as well?

                #print (sim_pars['upedge_lat'],sim_pars['upedge_lon'],sim_pars['mw'],sim_pars['stress'])
               
                parfile = open(outpath+base_initfiles+'_{0:.1f}'.format(m)+'_{0:.1f}'.format(s)+'.par', 'w')
                parfile.write(par_text.substitute(sim_pars))
                parfile.close()
                
                #TODO: write the name of the acc file along with the simulation parameters, for later processing


# In[ ]:




# In[37]:




# In[43]:




# In[90]:




# In[ ]:



