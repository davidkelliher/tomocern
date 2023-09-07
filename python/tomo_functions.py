from __future__ import division
import math 
import sys
import numpy as np
from subprocess import call
import os
import pylab as plt

PROTON_MASS = 938.27203e6
SPEED_OF_LIGHT = 299792458
BtoEnom = lambda B, rho, Erest, q: ((q*B*rho*SPEED_OF_LIGHT)**2 + (Erest**2))**0.5

def get_input(input_file="input_v2.dat"):

	f1 = open(input_file,'r')

	float_l = []
	for line in f1:
		if line[0] != '!':
			try:
				fl = float(line)
				float_l.append(fl)
			except:
				pass

	f1.close()
	
	input_struct = np.zeros(1, dtype=[('numframes','i2'),('numframesignore','i2'),('numbins','i2'),('framebinwidth','f8'),('binslowerignore','i2'),('binsupperignore','i2'),('binslowerempty','i2'),('binsupperempty','i2'),('synchbin','i2'),('profilecount','i2'),('profilelength','i2'),('dtbin','f8'),('dturns','i2'),('rebin','i2'),('y0','f2'),('dEmax','f8'),('filmstart','i2'),('filmstop','i2'),('filmstep','i2'),('iterations','i2'),('npartroot','i2'),('extend_phasespace','i2'),('beamref','i2'),('machineref','i2'),('VRF1ref','f8'),('VRF1dot','f8'),('VRF2ref','f8'),('VRF2dot','f8'),('h','i2'),('hratio','f8'),('phi12','f8'),('Bref','f8'),('Bdot','f8'),('Rnom','f8'),('rhonom','f8'),('gammatnom','f8'),('Erest','f8'),('q','i2'),('selffields','i2'),('coupling','f8'),('Zovern','f8'),('pickup','f8')])
	
	input_struct['numframes'] = float_l[0]
	input_struct['numframesignore'] = float_l[1]
	input_struct['numbins'] = float_l[2]
	input_struct['framebinwidth'] = float_l[3]
	input_struct['binslowerignore'] = float_l[5]
	input_struct['binsupperignore'] = float_l[6]
	input_struct['binslowerempty'] = float_l[7]
	input_struct['binsupperempty'] = float_l[8]
	input_struct['synchbin'] = float_l[10]

	input_struct['profilecount'] = float_l[0] - float_l[1]
	input_struct['profilelength'] = math.ceil((float_l[2] - (float_l[5] + float_l[6]))/float_l[9])
	input_struct['dtbin'] = float_l[3]
	input_struct['dturns'] = float_l[4]
	input_struct['rebin'] = float_l[9]
	input_struct['y0'] = input_struct['profilelength']/2
	input_struct['dEmax'] = float_l[11]
	input_struct['filmstart'] = float_l[12]
	input_struct['filmstop'] = float_l[13]
	input_struct['filmstep'] = float_l[14]
	input_struct['iterations'] = float_l[15]
	input_struct['npartroot'] = float_l[16]
	input_struct['extend_phasespace'] = float_l[17]
	input_struct['beamref'] = float_l[18]
	input_struct['machineref'] = float_l[19]
	input_struct['VRF1ref'] = float_l[20]
	input_struct['VRF1dot'] = float_l[21]
	input_struct['VRF2ref'] = float_l[22]
	input_struct['VRF2dot'] = float_l[23]
	input_struct['h'] = float_l[24]
	input_struct['hratio'] = float_l[25]
	input_struct['phi12'] = float_l[26]
	input_struct['Bref'] = float_l[27]
	input_struct['Bdot'] = float_l[28]
	input_struct['Rnom'] = float_l[29]
	input_struct['rhonom'] = float_l[30]
	input_struct['gammatnom'] = float_l[31]
	input_struct['Erest'] = float_l[32]
	input_struct['q'] = float_l[33]

	if int(float_l[34]) == 1:
		input_struct['selffields'] = True
	else:
		input_struct['selffields'] = False

	input_struct['coupling'] = float_l[35]
	input_struct['Zovern'] = float_l[36]
	input_struct['pickup'] = float_l[37]

	return input_struct

	
	
def write_inputfile(input_file, descrip, datafile, inputp):
	"""Write input_v2.dat file"""

	print ("write input file ",input_file)	
	f1 = open(input_file,'w')
	print(descrip, file=f1)
	for i in range(10):
		print('', file=f1)
	print("! Input data file =", file=f1)
	print(datafile, file=f1)
	print("! Directory in which to write all output =", file=f1)
	print("./", file=f1)
	print("! Number of frames in input data =", file=f1)
	print(inputp["numframes"][0], file=f1)
	print("! Number of frames to ignore =", file=f1)
	print(inputp["numframesignore"][0], file=f1)
	print("! Number of bins in each frame =", file=f1)
	print(inputp["numbins"][0], file=f1)
	print("! Width (in s) of each frame bin =", file=f1)
	print(inputp["framebinwidth"][0], file=f1)
	print("! Number of machine turns between frames =", file=f1)
	print(inputp["dturns"][0], file=f1)
	print("! Number of frame bins before the lower profile bound to ignore =", file=f1)
	print(inputp["binslowerignore"][0], file=f1)
	print("! Number of frame bins after the upper profile bound to ignore =", file=f1)
	print(inputp["binsupperignore"][0], file=f1)
	print("! Number of frame bins after the lower profile bound to treat as empty", file=f1)	
	print("! at the reconstructed time = ", file=f1)
	print(inputp["binslowerempty"][0], file=f1)
	print("! Number of frame bins after the lower profile bound to treat as empty", file=f1)	
	print("! at the reconstructed time = ", file=f1)
	print(inputp["binsupperempty"][0], file=f1)
	print("! Number of frame bins to rebin into one profile bin = ", file=f1)
	print(inputp["rebin"][0], file=f1)
	print("! Time (in frame bins) from the lower profile bound to the synchronous phase", file=f1)
	print("! (if <0, a fit is performed) in the bunch reference frame = ", file=f1)
	print(inputp["synchbin"][0], file=f1)
	print("! Max energy (in eV) of reconstructed phase space (if >0) = ", file=f1)
	print(inputp["dEmax"][0], file=f1)	
	print("! Number of the first profile at which to reconstruct =", file=f1)
	print(inputp["filmstart"][0], file=f1)
	print("! Number of the last profile at which to reconstruct =", file=f1)
	print(inputp["filmstop"][0], file=f1)
	print("! Step between reconstructions =", file=f1)
	print(inputp["filmstep"][0], file=f1)
	print("! Number of iterations for each reconstruction =", file=f1)
	print(inputp["iterations"][0], file=f1)
	print("! Square root of the number of test particles to track per cell =", file=f1)
	print(inputp["npartroot"][0], file=f1)
	print("! Flag to extend the region in phase space of map elements (if =1) =", file=f1)
	print(inputp["extend_phasespace"][0], file=f1)
	print("! Reference frame for bunch parameters (synchronous phase, baseline, integral) =", file=f1)
	print(inputp["beamref"][0], file=f1)
	print("! Reference frame for machine parameters (RF voltages, B-field) =", file=f1)
	print(inputp["machineref"][0], file=f1)
	print("!", file=f1)
	print("! Machine and Particle Parameters:", file=f1)
	print("! Peak RF voltage (in V) of principal RF system =", file=f1)
	print(inputp['VRF1ref'][0], file=f1)
	print("! and its time derivative (in V/s) =", file=f1)
	print(inputp['VRF1dot'][0], file=f1)
	print("! Peak RF voltage (in V) of higher-harmonic RF system =", file=f1)
	print(inputp['VRF2ref'][0], file=f1)
	print("! and its time derivative (in V/s) =", file=f1)
	print(inputp['VRF2dot'][0], file=f1)
	print("! Harmonic number of principal RF system =", file=f1)
	print(inputp['h'][0], file=f1)
	print("! Ratio of harmonics between RF systems =", file=f1)
	print(inputp['hratio'][0], file=f1)
	print("! Phase difference (in radians of the principal harmonic) between RF systems =", file=f1)
	print(inputp['phi12'][0], file=f1)
	print("! Dipole magnetic field (in T) =", file=f1)
	print(inputp["Bref"][0], file=f1)
	print("! and its time derivative (in T/s) =", file=f1)
	print(inputp["Bdot"][0], file=f1)
	print("! Machine radius (in m) =", file=f1)
	print(inputp["Rnom"][0], file=f1)
	print("! Bending radius (in m) =", file=f1)
	print(inputp["rhonom"][0], file=f1)
	print("! Gamma transition =", file=f1)
	print(inputp["gammatnom"][0], file=f1)
	print("! Rest mass (in eV/c**2) of accelerated particle =", file=f1)
	print(inputp["Erest"][0], file=f1)
	print("! Charge state of accelerated particle =", file=f1)
	print(inputp["q"][0], file=f1)
	print("!", file=f1)
	print("! Space Charge Parameters:", file=f1)
	print("! Flag to include self-fields in the tracking (if =1) =", file=f1)
	print("0", file=f1)
	print("! Geometrical coupling coefficient = ", file=f1)
	print(inputp['coupling'][0], file=f1)
	print("! Reactive impedance (in Ohms per mode number) over a machine turn =", file=f1)
	print(inputp['Zovern'][0], file=f1)
	print("! Effective pick-up sensitivity (in digitizer units per instantaneous Amp) =", file=f1)
	print(inputp['pickup'][0], file=f1)
	
	

def run_tomo_code(tomo_outdir, input_param):
	"""Prepare input file for tomography code. Run tomography code"""

	file_data, nturns, nbins, synch_index, dtbin, dturns, binslz, binsuz, recon_start, recon_stop, recon_step, Bref, Bdot, V_rf, V_rf2, harmonic, hratio, phi12, Rnom, rho, gamma_tr = input_param

	input_settings = get_input("input_default.dat")

	fname = 'input_v2.dat'

	descrip = 'Example'
	datafile = file_data#'rawdata.txt'
	input_settings["numframes"] = nturns
	input_settings["numbins"] = nbins
	input_settings["synchbin"] = synch_index
	input_settings["framebinwidth"] = dtbin
	input_settings["dturns"] = dturns
	input_settings["binslowerignore"] = binslz
	input_settings["binsupperignore"] = binsuz
	input_settings["filmstart"] = recon_start
	input_settings["filmstop"] = recon_stop
	input_settings["filmstep"] = recon_step
	input_settings['Bref'] = Bref
	input_settings['Bdot'] = Bdot
	input_settings['VRF1ref'] = V_rf
	input_settings['VRF2ref'] = V_rf2
	input_settings['h'] = harmonic
	input_settings['hratio'] = hratio
	input_settings['phi12'] = phi12
	input_settings['Rnom'] = Rnom
	input_settings['rhonom'] = rho
	input_settings['gammatnom'] = gamma_tr
	input_settings['pickup'] = 1.0
	input_settings['selffields'] = 0
	
	
	write_inputfile(fname, descrip, datafile, input_settings)

	temp_files = ['plotinfo.data','profiles.data','d001.data','d100.data','jmax.data','image001.data','image100.data']


	#run tomography code
	call(["tomo"])

	files = os.listdir(os.curdir)
	for fn in files:
		if fn[-4:] == 'data':
			os.rename(fn,tomo_outdir+"/"+fn)

	return
	
	
def get_plotinfo(input_file="plotinfo.data"):

	f1 = open(input_file,'r')

	float_l = []
	for line in f1:
		if line[0] == ' ':
			#print line
			lspl = line.split("=")
			try:
				
				fl = float(lspl[1])
				float_l.append(fl)
			except:
				pass

	f1.close()


	input_struct = np.zeros(1, dtype=[('profilecount','i2'),('profilelength','i2'),('dtbin','f8'),('dEbin','f8'),('eperimage','f8'),('xat0','f8'),('yat0','f8'),('phi_s','f8')])


	input_struct['profilecount'] = float_l[0]
	input_struct['profilelength'] = float_l[1]
	input_struct['dtbin'] = float_l[2]
	input_struct['dEbin'] = float_l[3]
	input_struct['eperimage'] = float_l[4]
	input_struct['xat0'] = float_l[5]
	input_struct['yat0'] = float_l[6]
	input_struct['phi_s'] = float_l[10]

	return input_struct


def file_read(filen, sel_col=0):
	""" filen: file path.
		sel_col: indext of column in data file to read. Default 0"""
	
	f1 = open(filen,'r')
	data_l = []
	for data in f1:
		dataspl = data.split()
		data_l.append(float(dataspl[sel_col]))
	data_a = np.array(data_l)

	return data_a
	
def mountain_range(data, nt, nslices, bucketlim_low=None, bucketlim_high=None, sync_bin = None, incl_xaxis = False, xaxis = None, xlabel='time bin',fmount='mountainrange', extra_data = None):
	"""Plot 2D data which has shape (nslices, nt) in mountain range style
		extra_data - 1D data to superimpose on figure """

	#nd = data.shape
	nd = len(data)
	ndx = nslices

	#establish range
	max_data_t = np.array([max(data[it]) for it in range(nt)])
	min_data_t = np.array([min(data[it]) for it in range(nt)])
	range_data_t = max_data_t - min_data_t
	max_range = np.max(range_data_t)
	
	ampfac = 2000
	
	#plt.plot(range_data_t)
	#plt.show()

	lw_a = np.linspace(0.8, 0.0, nt)
	lw_a = 0.4 - np.array(list(range(nt)))/6000
	

	fig = plt.figure(figsize=(8, 8), facecolor='w')
	ax = plt.subplot(111, frameon=False)
	
	if not incl_xaxis:
		xaxis = list(range(1,ndx+1))

	for it in range(nt):
		ax.plot(xaxis, ampfac*data[it] + it,'k',lw = lw_a[it])

	ax.set_ylabel('turn number')

	if bucketlim_low != None:
		plt.axvline(x=bucketlim_low+0.5,color='gray')

	if bucketlim_high != None:
		plt.axvline(x=bucketlim_high+0.5,color='gray')

	if sync_bin != None:
		plt.axvline(x=sync_bin, color='r', linestyle='--')

	if extra_data != None:
		col_l = ['r','m','b']
		i = 0
		for ex_data in extra_data:
			plt.plot(ex_data, np.array(list(range(nt))), col_l[i])
			i = i + 1

	#show turn number of right axis
	ax.set_xlabel(xlabel)		
	ax.set_ylabel('turn number')
	plt.savefig(fmount)
	plt.show()
	
	return
	
def set_param(inputp):
	"""Calculate some parameters"""

	Erest = inputp['Erest'][0]
	gammat = inputp['gammatnom'][0]
	Rnom = inputp['Rnom'][0]
	E0 =  BtoEnom(inputp['Bref'][0], inputp['rhonom'][0], Erest, inputp['q'][0]) #total energy
	gamma0 = E0/Erest
	eta0 = 1/(gamma0**2) - 1/(gammat**2) 
	beta0 = (1 - (1/gamma0**2))**0.5
	omegarev0 = beta0*SPEED_OF_LIGHT/Rnom
	Vpeak = inputp['VRF1ref'][0]
	dtbin = inputp['dtbin'][0]
	
	return E0, gamma0, eta0, beta0, omegarev0, Vpeak, dtbin
	
def list_files(path):
	"""list all files in data directory"""
	def sorted_ls(path):
		mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
		return list(sorted(os.listdir(path), key=mtime))

	def first_chars(x):
		return(x[:40])

	sort_by_name = True
	if sort_by_name:
		files = sorted(os.listdir(path), key = first_chars)
	else:
		files = sorted_ls(path)

	return files
