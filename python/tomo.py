from __future__ import division
import math 
import sys
import numpy as np
import pylab as plt
import tomo_functions as tomo_fns
import os 
import glob


run_tomo_code = False

nturns = 58 #number of frames or turns of data
nbins = 500 #number of bins in each frame
dtbin = 1E-9 #duration of each bin

dturns = 40 #number of machine turns between frames

#number of bins to ignore at either end of profile
binslz = 0
binsuz = 20 

synch_index = 248 #number of bins from lower profile bound to the synchronous bin

recon_start = 1 #Number of the first profile at which to reconstruct
recon_stop = 58 #Number of the last profile at which to reconstruct
recon_step = 20 #Step between reconstructions

#Principal harmonic
V_rf = 7.0E3
harmonic = 1

#second harmonic system
V_rf2 = 0
hratio = 1
phi12 = 0 

Bref = 0.4585
Bdot = 2.2
Rnom = 25
rho = 8.239
gamma_tr = 4.1
alpha_c = (1/gamma_tr)**2

ndim = nbins - (binslz + binsuz)

data_file = "psbnotch.dat"

phi_s = None	
#specify Bdot or synchronous phase phi_s and set the other to None
if Bdot != None:
	phi_s = np.arcsin((2*rho*math.pi*Rnom*Bdot)/V_rf)
else:
	Bdot = (V_rf*np.sin(phi_s))/(rho*2*math.pi*Rnom)
#calculate synchrotron tune
#----------------------------
PROTON_MASS = 0.93827231E9
SPEED_OF_LIGHT = 299792458
mass = PROTON_MASS
#function relating B to total Energy given the bending radius and rest mass
BtoEnom = lambda B, rho, Erest, q: ((q*B*rho*SPEED_OF_LIGHT)**2 + (Erest**2))**0.5
#synchrotron tune
Qs_fn = lambda phi0, V0, h, beta, E0, eta: np.sqrt(np.abs(eta*np.cos(phi0)) * h* V0 / (2 * math.pi * (beta**2) * E0))

Etot = BtoEnom(Bref, rho, mass, 1) #total energy
print ("K.E (MeV) ",1e-6*(Etot-mass))

p0 = np.sqrt(Etot**2 - mass**2)
betarel = p0/Etot
betagam = p0/mass
gamma = betagam/betarel

eta = (alpha_c - gamma**(-2)) #phase_slip
print ("gamma, eta ",gamma, eta)

print ("phi_s ",phi_s)
Qs = Qs_fn(phi_s, V_rf, harmonic, betarel, Etot, eta)
Ts = 1/Qs
yfac = np.sqrt(1 - 0.5*(math.pi - 2*phi_s)*np.tan(phi_s))
bh = ((2*Qs)/(harmonic*abs(eta))*yfac) #bucket height

print ("Synchrotron tune, period ",Qs, 1/Qs)
print ("bh in dp/p ",bh)	

	
#----------------------------	



#set tomography output directory
tomo_outdir = "tomo_output"
files = os.listdir(os.curdir)
if tomo_outdir not in files:
	print("create output directory "+tomo_outdir)
	os.mkdir(tomo_outdir)
	
#set tomography output directory for image files
tomo_outdir_images = "tomo_images"
files = os.listdir(os.curdir)
if tomo_outdir_images not in files:
	print("create output directory "+tomo_outdir_images)
	os.mkdir(tomo_outdir_images)
	


if run_tomo_code:
	#remove existing output files
	output_files =  glob.glob(tomo_outdir+'/*') 
	for f in output_files:
		os.remove(f)
		
		
tomo_code_param = [data_file, nturns, nbins, synch_index, dtbin, dturns, binslz, binsuz, recon_start, recon_stop, recon_step, Bref, Bdot, V_rf, V_rf2, harmonic, hratio, phi12, Rnom, rho, gamma_tr]

if run_tomo_code:
	tomo_fns.run_tomo_code(tomo_outdir, tomo_code_param)

outfiles = os.listdir(tomo_outdir)
if any("image" in f1 for f1 in outfiles):
	run_success = True
else:
	run_success = False
print("run_sucess ",run_success)


if run_success:

	#list files alphabetically
	output_files = tomo_fns.list_files(tomo_outdir)
	
	#print ("store_files ",store_files)
	
	image_files = []
	plotinfo_files = []

	for f1 in output_files:			
		if "image" in f1:
			image_files.append(f1)
		if "plotinfo" in f1:
			plotinfo_files.append(f1)
	
	print("image files ",image_files)
	#print(store_plotinfo_files)
	
	rec_img_l = []
	for file1 in image_files:
		fspl = (file1.split("."))
		rec_img_id_str = fspl[0][5:]
		rec_img_id = int(rec_img_id_str)
		rec_img_l.append(rec_img_id)
		
	print ("rec_img_l ",rec_img_l)
		
	plotinfo = tomo_fns.get_plotinfo(input_file=tomo_outdir+"/plotinfo.data")
	
	x0 = plotinfo['xat0'][0]
	y0 = plotinfo['yat0'][0]
	dEbin =  plotinfo['dEbin'][0]
	dtbin =  plotinfo['dtbin'][0]
	phi_s_out = plotinfo['phi_s'][0]
	
	print ("x0, y0, dEbin, dtbin, phi_s_out ",x0, y0, dEbin, dtbin, phi_s_out)

	
	#input data 
	fname_prof = 'profiles.data'
	profiles_norm_raw = tomo_fns.file_read(tomo_outdir+'/'+fname_prof)
	profiles_norm_data_a = np.array(np.split(profiles_norm_raw, nturns))

	#reconstructed data 
	fname_prof = 'profiles_recon.data'
	profiles_norm_recon_raw = tomo_fns.file_read(tomo_outdir+'/'+fname_prof)
	profiles_norm_recon_a = np.array(np.split(profiles_norm_recon_raw, nturns))
	
	dim = profiles_norm_data_a.shape
	
	tomo_fns.mountain_range(profiles_norm_data_a, nturns, dim[1], fmount='mountainrange_data')
	tomo_fns.mountain_range(profiles_norm_recon_a, nturns, dim[1], fmount='mountainrange_recon')
	
	inputp = tomo_fns.get_input()
	E0, gamma0, eta0, beta0, omega_rev0, Vpeak, dtbin = tomo_fns.set_param(inputp)
	print ("E0 ",E0, Etot)
	
	T_rf = 2*np.pi*harmonic/omega_rev0
	T_rf_ns = T_rf*1e9
	
	nloop = len(rec_img_l)

	for iloop in range(nloop):

		irc = rec_img_l[iloop]

		print("index, image and turn no. ",iloop, irc, (irc-1)*dturns + 1)

		fname_discrep = tomo_outdir+'/d000'[:-len(str(irc))]+ str(irc) + '.data'
		discrep_data = tomo_fns.file_read(fname_discrep, sel_col=1)
		discrep_final = discrep_data[-1]
		print("discrepancy ",discrep_final)  
			
		fname_img = tomo_outdir+'/image000'[:-len(str(irc))]+ str(irc) + '.data'
		print("image file name ",fname_img)
	
		#reconstruction data
		image_data_a = tomo_fns.file_read(fname_img) #1D data
		ndim_image =  int(np.sqrt(len(image_data_a)))
		
		#split into array
		image_data_spl = np.split(image_data_a, ndim_image)
		#transpose to get phase on horizontal axis, energy on vertical
		image_data_tr = np.transpose(image_data_spl)

		#vertical axis
		delta_E_a = dEbin*(np.array(list(range(ndim_image))) - y0)

		#horizontal axis
		t_synch = phi_s_out*T_rf/(2*math.pi)
		t_tomo_ns = 1e9*dtbin*(np.array(list(range(ndim))) - x0) + 1e9*t_synch
		phase_tomo = t_tomo_ns*(2*math.pi)/T_rf_ns
		phase_tomo_deg = phase_tomo*180/math.pi		
		
		plt.contour(phase_tomo_deg, delta_E_a*1e-6, image_data_tr)
		plt.ylabel(r'$\Delta$E [MeV]')
		plt.xlabel(r'$\phi$ [ns]')
		plt.show()
	
	
