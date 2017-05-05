import numpy as np
import matplotlib.pyplot as plt

def rk2integrate(eta, Mse, del_eta, func):
	n = int(eta/del_eta)
	v = np.zeros(n+1)
	eps = 0.01
	for i in range(1, n+1):
		k1 = func(Mse, v[i-1], eps)*del_eta
		k2 = func(Mse, v[i-1] + k1, eps)*del_eta
		v[i] = v[i-1] + (k1 + k2)*0.5
	return v

def gen_e(Mse, V, eps):
	return (2*Mse**2*((1+2*V/Mse**2)**0.5 -1) + 2*np.exp(-V)- 2 + eps**2)**0.5


def gennormv(eta, del_eta, Mse = 1):
	return rk2integrate(eta, Mse, del_eta, gen_e)

def bohmnormv(eta, del_eta, Mse = 1):
	n = int(eta/del_eta)
	v = np.zeros(n+1)
	eps = 0.01
	for i in range(1, n+1):
		v[i] = del_eta*i*eps
	return v

def chlanormv(eta, Mse = 1):
	n = len(eta)
	v = np.zeros(n)
	for i in range(0, n):
		v[i] = 3.0**(4.0/3.0)/(2.0**(5.0/3.0))*eta[i]**(4.0/3.0)
	return v



def plotfig(name, title, xlabel, ylabel, xdata, ydata, legend):
	plt.figure()
	plt.title(title)
	for i in range(len(xdata)):
		plt.plot(xdata[i], ydata[i], label = legend[i])
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.legend(loc = 2)
	plt.savefig(name)
	plt.close()

def question3():
	etas = np.arange(0, 20.001, 0.01)
	v_gen = gennormv(20, 0.01)
	e_gen = gen_e(1, v_gen, 0.01)
	vbohm = bohmnormv(20, 0.01)
	e_bohm = np.linspace(0.01, 0.01, len(vbohm))
	vchla = chlanormv(etas)
	e_chla = 2**(3.0/4.0)*vchla**(1.0/3.0)
	plotfig("Gen_Vnorm.png", "Normalised Sheath potential", "eta", "V", [etas, etas, etas], [v_gen, vbohm, vchla], ["Generalised Case", "Bohm's approx.", "Child-Langmuir approx."])
	plotfig("Gen_Enorm.png", "Normalised Sheath Electric Field", "Eta", "E", [etas, etas, etas], [e_gen, e_bohm, e_chla], ["Generalised Case", "Bohm's approx.", "Child-Langmuir approx."])

def question4():
	del_eta = 0.01
	for eps in [0.1, 0.01, 0.001]:
		eta = 0
		v = del_eta*eps
		while (v < 2.8):
			k1 = gen_e(1, v, eps)*del_eta
			k2 = gen_e(1, v+k1, eps)*del_eta
			v += 0.5*(k1 + k2)
			eta +=del_eta
		print "For epsilon = "+str(eps)+" sheath thickness(normalised) =", eta
	return eta

def question5():
	print gen_e(1, 2.8, 0.1), gen_e(1, 2.8, 0.01), gen_e(1, 2.8, 0.001)

def question6():
	eps = 0.01
	del_eta = 0.01
	for i in [1, 10, 100]:
		vmax = 1.6021e-19*i/(1.3806e-23*1e4)
		eta = 0
		v = del_eta*eps
		while (v < vmax):
			k1 = gen_e(1, v, eps)*del_eta
			k2 = gen_e(1, v+k1, eps)*del_eta
			v += 0.5*(k1 + k2)
			eta += del_eta
			if (abs(v - eta*eps) > abs(v - 3.0**(4.0/3.0)/(2.0**(5.0/3.0))*eta**(4.0/3.0))):
				print "Child- Langmuir region", eta, v
		print "For electrode potential = "+str(eps)+" sheath thickness(normalised) =", eta
	return eta

question3()
question4()
question5()
question6()