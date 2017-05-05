import numpy as np 
import matplotlib.pyplot as plt
import math
from matplotlib.backends.backend_pdf import PdfPages
import operator as op


N_avgad = 6.022e+23
k_boltz = 1.38e-23

pp = PdfPages('A1_results.pdf')

def ncr(n, r):    #http://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom


def pdf_velf(vx, vy, vz, mol_wt, rho, temp):
	nden = rho/mol_wt
	v = (vx**2 + vy**2 + vz**2)**0.5
	return nden*(mol_wt/(2*np.pi*k_boltz*temp))**1.5*math.exp(-1*mol_wt*v**2/(2*k_boltz*temp))

def pdf_speedf(v, mol_wt, rho, temp):
	#print mol_wt
	nden = rho/mol_wt
	return nden*(mol_wt/(2*np.pi*k_boltz*temp))**1.5*4*np.pi*v**2*math.exp(-1*mol_wt*v**2/(2*k_boltz*temp))

def pdf_energyf(eps, mol_wt, rho, temp):
	nden = rho/mol_wt
	return nden*np.pi/(np.pi*k_boltz*temp)**1.5*eps**0.5*math.exp(-eps/(k_boltz*temp))

def vmeanf(v, pdf_speed, nden):
	vmean = 0
	for i in range(len(v)):
		vmean += v[i]*(v[i] - v[i-1]) * pdf_speed[i]
	return vmean/nden


def pressuref(mol_wt, v, pdf_speed):
	pressure = 0
	for i in range(len(pdf_speed)):
		pressure += mol_wt * pdf_speed[i] * (v[i] - v[i-1]) * v[i]**2 / 3
	return pressure

def int_energyf(eps, pdf_energy):
	energy = 0
	for i in range(len(pdf_energy)):
		energy += eps[i] * pdf_energy[i]*(eps[i] - eps[i-1])
	return energy

def entropyf(v, pdf_speed):
	entropy = 0
	for i in range(len(v)):
		if(pdf_speed[i] != 0.0):
			entropy += pdf_speed[i]*math.log(pdf_speed[i])*(v[i] - v[i-1])
	return entropy

def exponential(v):
	return np.exp(-1*v)

def exp_log(x, p, beta):
	return -1/math.log(p)*beta*(1-p)*np.exp(-beta*x)/(1 - (1 - p)*np.exp(-beta*x))

def levy(v, mu, c):
	return (c/2*np.pi)**0.5*math.exp(-c/(2*(v - mu)))/ (v - mu)**1.5


def plotgraph(title, name, xdata, ydata, vlines, vlinmax, label, vlincolor, xlabel, ylabel):
	plt.figure()
	plt.title(title)
	plt.plot(xdata, ydata)
	for i in range(len(vlines)):
		plt.vlines(vlines[i], 0, vlinmax[i], label = label[i], color = vlincolor[i], linestyle = 'dashed')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.legend()
	plt.savefig(name)
	plt.close()

def plotq8(title, name, xdata, ydata, xlabel, ylabel, flag):
	fig, ax1 = plt.subplots()
	#ax1.title(title)
	ax1.plot(xdata[0], ydata[0], 'b-')
	ax1.set_ylabel(ylabel[0])
	ax1.set_xlabel(xlabel)
	if flag == 0:
		ax1.set_xscale('symlog')
	for tl in ax1.get_yticklabels():
		tl.set_color('b')

	ax2 = ax1.twinx()
	ax2.plot(xdata[1], ydata[1], 'r-')
	ax1.set_ylabel(ylabel[1])
	ax1.set_xlabel(xlabel)
	if flag == 0:
		ax1.set_xscale('symlog')
	for tl in ax2.get_yticklabels():
		tl.set_color('r')

	plt.savefig(name)

def question2_3_4_5(mol_wt, rho, temp, ngrids):    #mol_wt in amu
	mol_wt = mol_wt/N_avgad *10**-3
	sigma = (k_boltz*temp/mol_wt)**0.5
	vx = np.linspace(-5*sigma, 5*sigma, ngrids)
	vy = np.linspace(-5*sigma, 5*sigma, ngrids)
	vz = np.linspace(-5*sigma, 5*sigma, ngrids)
	v = np.linspace(0, 10*sigma, ngrids*100)
	eps = 0.5*mol_wt*v**2
	
	pdf_vel = np.zeros((len(vx), len(vy), len(vz)))
	pdf_speed = np.zeros(len(v))
	pdf_energy = np.zeros(len(eps))
	for i in range(len(vx)):
		for j in range(len(vy)):
			for k in range(len(vz)):
				pdf_vel[i,j,k] = pdf_velf(vx[i], vy[j], vz[k], mol_wt, rho, temp)

	for i in range(len(v)):
		pdf_speed[i] = pdf_speedf(v[i], mol_wt, rho, temp)
		pdf_energy[i] = pdf_energyf(eps[i], mol_wt, rho, temp)

	v_mean = (8*k_boltz*temp/(np.pi * mol_wt))**0.5
	y1 = pdf_velf(v_mean, 0, 0, mol_wt, rho, temp)
	y2 = pdf_speedf(v_mean, mol_wt, rho, temp)
	v_rms = (3*k_boltz*temp / mol_wt)**0.5
	y3 = pdf_velf(v_rms, 0, 0, mol_wt, rho, temp)
	y4 = pdf_speedf(v_rms, mol_wt, rho, temp)
	e_mean = 1.5*k_boltz*temp
	y5 = pdf_energyf(e_mean, mol_wt, rho, temp)
	plotgraph("Velocity Distribution function","f_vel.png", vx, pdf_vel[:, ngrids/2, ngrids/2], [v_mean, v_rms], [y1, y3], ["Vmean", "Vrms"], ['g', 'r'], "vx", "f(vx, vy, vz)")
	plotgraph("Speed Distribution function","f_speed.png", v, pdf_speed, [v_mean, v_rms], [y2, y4], ["Vmean", "Vrms"], ['g', 'r'], "v", "f(v)")
	plotgraph("Energy Distribution function","f_energy.png", eps, pdf_energy, [e_mean], [y5], ["Emean"], ['g'], "energy", "f(e)")
	pressure = pressuref(mol_wt, v, pdf_speed)
	energy = int_energyf(eps, pdf_energy)
	nden = rho / mol_wt
	entropy = entropyf(v, pdf_speed)
	vmean = vmeanf(v, pdf_speed, nden)
	print  pressure, energy, vmean, entropy
	return pressure, energy, vmean, entropy

def question6(mol_wt, rho, temp):       #mol_wt in amu
	pressure_id = rho*k_boltz*N_avgad*temp/ mol_wt * 10**3
	#print pressure_id
	#print energy_id
	v_mean_id = (8*k_boltz*temp/(np.pi * mol_wt)*N_avgad*10**3)**0.5
	nden = rho/mol_wt * N_avgad * 1000
	energy_id = 1.5*nden*k_boltz*temp
	#entropy_id = nden*mol_wt/N_avgad*10**-3 / (2*k_boltz*temp)*(np.log(nden*(mol_wt*10**-3 / (N_avgad * 2*np.pi * k_boltz*temp))**1.5)+1.5)
	#print entropy_id
	ngrids = [50, 100, 200]
	error_p = np.zeros(len(ngrids))
	error_e = np.zeros(len(ngrids))
	error_v = np.zeros(len(ngrids))
	#error_s = np.zeros(len(ngrids))
	for i in range(len(ngrids)):
		pressure, energy, vmean, entropy = question2_3_4_5(mol_wt, rho, temp, ngrids[i])
		error_p[i] = abs(pressure - pressure_id)
		error_e[i] = abs(energy - energy_id)
		error_v[i] = abs( vmean- v_mean_id)
		#error_s[i] = abs(entropy - entropy_id)
	plotgraph("Error in Pressure", "Error_p.png", ngrids, error_p, [], [], [], [], "ngrids", "error")
	plotgraph("Error in Energy", "Error_e.png", ngrids, error_e, [], [], [], [], "ngrids", "error")
	plotgraph("Error in mean velocity", "Error_v.png", ngrids, error_v, [], [], [], [], "ngrids", "error")
	#plotgraph("Error in Entropy", "Error_s.png", ngrids, error_s, [], [], [], [], "ngrids", "error")

def question7(mol_wt, rho, temp):
	ngrids = 1000
	mol_wt = mol_wt/N_avgad *10**-3
	nden = rho / mol_wt
	sigma = (k_boltz*temp/mol_wt)**0.5
	v = np.linspace(0.001, 10*sigma, ngrids*100)
	pdf_max_boltz = np.zeros(len(v))
	pdf_exp =  np.zeros(len(v))
	pdf_exp_log = np.zeros(len(v))
	pdf_levy =  np.zeros(len(v))
	for i in range(len(v)):
		pdf_max_boltz[i] = pdf_speedf(v[i], mol_wt, rho, temp)
		pdf_exp[i] = exponential(v[i])
		pdf_exp_log[i] = exp_log(v[i], 0.5, 2)
		pdf_levy[i] = levy(v[i], 0, 2)

	entropy_max = entropyf(v, pdf_max_boltz)
	entropy_exp = entropyf(v, pdf_exp*nden)
	entropy_exp_log = entropyf(v, pdf_exp_log*nden)
	entropy_levy = entropyf(v, pdf_levy*nden)

	return entropy_max, entropy_exp, entropy_exp_log, entropy_levy

def question8(rho, temp):
	m_p = 1.67e-27
	m_e = 9.11e-31
	sigma = (k_boltz*temp/m_p)**0.5
	sigma2 = (k_boltz*temp/m_e)**0.5
	ngrids = 100
	vx1 = np.linspace(-5*sigma, 5*sigma, ngrids)
	vx2 = np.linspace(-5*sigma2, 5*sigma2, ngrids)
	vy = np.linspace(-5*sigma, 5*sigma, ngrids)
	vz = np.linspace(-5*sigma, 5*sigma, ngrids)
	v = np.linspace(0, 10*sigma, ngrids*100)
	v2 = np.linspace(0, 10*sigma2, ngrids*100)
	eps1 = 0.5*m_p*v**2
	eps2 = 0.5*m_e*v2**2
	
	pdf_vel1 = np.zeros((len(vx1), len(vy), len(vz)))
	pdf_vel2 = np.zeros((len(vx2), len(vy), len(vz)))
	pdf_speed1 = np.zeros(len(v))
	pdf_speed2 = np.zeros(len(v))
	pdf_energy1 = np.zeros(len(eps1))
	pdf_energy2 = np.zeros(len(eps2))
	for i in range(len(vx1)):
		for j in range(len(vy)):
			for k in range(len(vz)):
				pdf_vel1[i,j,k] = pdf_velf(vx1[i], vy[j], vz[k], m_p, rho, temp)
				pdf_vel2[i,j,k] = pdf_velf(vx2[i], vy[j], vz[k], m_e, rho, temp)

	for i in range(len(v)):
		pdf_speed1[i] = pdf_speedf(v[i], m_p, rho, temp)
		pdf_speed2[i] = pdf_speedf(v2[i], m_e, rho, temp)
		pdf_energy1[i] = pdf_energyf(eps1[i], m_p, rho, temp)
		pdf_energy2[i] = pdf_energyf(eps2[i], m_e, rho, temp)

	plotq8("Velocity distribution", "Plasma_vel.png", [vx1, vx2], [pdf_vel1[:, 50, 50], pdf_vel2[:, 50, 50]], "log(v)", ["f(v_p)", "f(v_e)"], 0)
	plotq8("Speed Distribution", "Plasma_speed.png", [v, v2], [pdf_speed1, pdf_speed2], "log(v)", ["f(v_p)", "f(v_e)"], 0)
	plotq8("Energy distribution", "Plasma_energy.png", [eps1[0:ngrids*50], eps2[0:ngrids*50]], [pdf_energy1[0:ngrids*50], pdf_energy2[0:ngrids*50]], "e", ["f(ve)", "f(e)"], 1)





question6(28.1, 1.225, 300)
e1, e2, e3, e4 = question7(28.1, 1.225, 300)
question8(1.225, 300)
