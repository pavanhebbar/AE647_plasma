import numpy as np
import matplotlib.pyplot as plt

k_boltz = 1.38064e-23
e = 1.60271e-19
eps = 8.85419e-12
mp = 1.67262e-27
me = 9.10938e-31
mu0 = 4*np.pi*1e-7

class plasma:
	def __init__(self, n, kT, B, mass):
		self.n = n
		self.kT = kT
		self.b0 = B
		self.mass = mass

	def len_1(self):
		return debye_len(self.n, self.kT)/rho_c(self.mass, e, self.kT, self.b0)

	def len_2(self):
		return debye_len(self.n, self.kT)/mean_fp(self.n)

	def len_3(self):
		return self.n*debye_len(self.n, self.kT)

	def time_1(self):
		return omega_p(self.n, self.mass)/omega_c(e, self.mass, self.b0)

	def time_2(self):
		return beta(self.n, self.kT, self.b0)

def debye_len(n, kT): #kT in eV
	return (eps*kT*e/(n*e**2)*0.5)**0.5

def rho_c(mass, charge, kT, b0):
	return mass/(charge*b0)*(2*kT*e/mass)**0.5

def mean_fp(n):
	return (3/(4*np.pi*n))**(1/3)

def omega_p(n, mass):
	return (n*e**2/(mass*eps))**0.5

def omega_c(charge, mass, b0):
	return charge*b0/mass

def beta(n, kT, b0):
	return n*kT*e/(b0**2/(2*mu0))

def boris_push(charge, mass, e_field_amp, e_field_freq, b0, v0, tstep, ntime):
    v = np.zeros((ntime+1, 3))
    x = np.zeros((ntime + 1, 3))
    mu = np.zeros(ntime+1)
    mu[0] = 0.5*mass*np.sum(v0**2)/(2*b0)
    vd_constE = np.array([0, e_field_amp[0]/b0, 0])
    x[0] = [0, 0, 0]
    v[0] = v0
    for i in range(1, ntime+1):
        vp = np.zeros(3)
        alp = float (charge*b0/(2*mass)*tstep)
        e_field = e_field_amp*np.cos(e_field_freq/(2*np.pi)*i*tstep)
        vm = v[i-1] + alp*e_field/b0
        vp[0] = vm[0]*(1-alp**2)/(1+alp**2) + 2*alp/(1+alp**2)*vm[1]
        vp[1] = vm[1]*(1-alp**2)/(1+alp**2) - 2*alp/(1+alp**2)*vm[0]
        vp[2] = 0
        v[i] = vp + alp*e_field/b0
        x[i] = x[i-1] + v[i]*tstep
        vper = v[i] - vd_constE
        mu[i] = 0.5*mass*np.sum(vper**2)/(2*b0)
    return x, v, mu

def rk2_q3(charge, mass, e_field, n, temp, tstep, ntime):
    v = np.zeros((ntime+1, 3))
    x = np.zeros((ntime + 1, 3))
    vc = 4.8e-14*n*15*temp**-1.5
    x[0] = [0, 0, 0]
    v[0] = np.zeros(3)
    for i in range(1, ntime + 1):
        kv1 = (charge*e_field - vc*mass*v[i-1])/mass*tstep
        kx1 = v[i-1]*tstep
        kv2 = (charge*e_field - vc*mass*(v[i-1]+kv1/2))/mass*tstep
        kx2 = (v[i-1] + kv1)*tstep
        v[i] = v[i-1] + kv2
        x[i] = x[i-1] + kx2
    return x, v

def plotfig(name, title, xlabel, ylabel, xdata, ydata, legend, scale = "lin", loct = 2):
    plt.figure()
    plt.title(title)
    if scale == "log":
        plt.xscale('log')
        plt.yscale('log')
    for i in range(len(xdata)):
        plt.plot(xdata[i], ydata[i], label = legend[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc = loct)
    plt.savefig(name)
    plt.close()


def ques1():
	solar_e = plasma(1e7, 10, 2.8e-9, me)
	sunspot_e = plasma(1e20, 1, 3, me)
	fu_reac_e = plasma(1e20, 10000, 5, me)
	solar_p = plasma(1e7, 10, 2.8e-9, mp)
	sunspot_p = plasma(1e20, 1, 3, mp)
	fu_reac_p = plasma(1e20, 10000, 5, mp)
	lens_no = [1, 2, 3]
	solarp_len = [solar_p.len_1(), solar_p.len_2(), solar_p.len_3()]
	sunspotp_len = [sunspot_p.len_1(), sunspot_p.len_2(), sunspot_p.len_3()]
	fup_len = [fu_reac_p.len_1(), fu_reac_p.len_2(), fu_reac_p.len_3()]
	solare_len = [solar_e.len_1(), solar_e.len_2(), solar_e.len_3()]
	sunspote_len = [sunspot_e.len_1(), sunspot_e.len_2(), sunspot_e.len_3()]
	fue_len = [fu_reac_e.len_1(), fu_reac_e.len_2(), fu_reac_e.len_3()]
	times_no = [1, 2]
	solarp_time = [solar_p.time_1(), solar_p.time_2()]
	sunspotp_time = [sunspot_p.time_1(), sunspot_p.time_2()]
	fup_time = [fu_reac_p.time_1(), fu_reac_p.time_2()]
	solare_time = [solar_e.time_1(), solar_e.time_2()]
	sunspote_time = [sunspot_e.time_1(), sunspot_e.time_2()]
	fue_time = [fu_reac_e.time_1(), fu_reac_e.time_2()]
	plt.figure()
	plt.title("Comparision between length scales of ions")
	plt.xlabel("Length scale number")
	plt.ylabel("value")
	plt.xlim(0, 4)
	plt.yscale('log')
	plt.plot(lens_no, solarp_len, 'o', label = 'solar wind')
	plt.plot(lens_no, sunspotp_len, 'o', label = 'sunspot')
	plt.plot(lens_no, fup_len, 'o', label = 'fusion reactor')
	plt.legend(loc = 2)
	plt.savefig('ques1_len.png')
	plt.close()
	plt.figure()
	plt.title("Comparision between time scales of ions")
	plt.xlabel("Time scale number")
	plt.ylabel("value")
	plt.xlim(0, 4)
	plt.yscale('log')
	plt.plot(times_no, solarp_time, 'o', label = 'solar wind')
	plt.plot(times_no, sunspotp_time, 'o', label = 'sunspot')
	plt.plot(times_no, fup_time, 'o', label = 'fusion reactor')
	plt.legend(loc = 2)
	plt.savefig('ques1_time.png')
	plt.close()
	plt.figure()
	plt.title("Comparision between length scales of electrons")
	plt.xlabel("Length scale number")
	plt.ylabel("value")
	plt.xlim(0, 4)
	plt.yscale('log')
	plt.plot(lens_no, solare_len, 'o', label = 'solar wind')
	plt.plot(lens_no, sunspote_len, 'o', label = 'sunspot')
	plt.plot(lens_no, fue_len, 'o', label = 'fusion reactor')
	plt.legend(loc = 2)
	plt.savefig('ques1_len2.png')
	plt.close()
	plt.figure()
	plt.title("Comparision between time scales of electrons")
	plt.xlabel("Time scale number")
	plt.ylabel("value")
	plt.xlim(0, 4)
	plt.yscale('log')
	plt.plot(times_no, solare_time, 'o', label = 'solar wind')
	plt.plot(times_no, sunspote_time, 'o', label = 'sunspot')
	plt.plot(times_no, fue_time, 'o', label = 'fusion reactor')
	plt.legend(loc = 2)
	plt.savefig('ques1_timee.png')
	plt.close()


def expec2b(charge, mass, e_field_amp, e_field_freq, b0, v0, tstep, ntime):
    v = np.zeros((ntime+1, 3))
    x = np.zeros((ntime + 1, 3))
    omega_c = charge*b0/mass
    omega = e_field_freq/(2*np.pi)
    x[0] = [0, 0, 0]
    v[0] = v0
    v0 = v0[0]
    e_field_amp = e_field_amp[0]
    for i in range(ntime + 1):
    	x[i][0] = v0/omega_c*np.sin(omega_c*i*tstep) - charge/mass*e_field_amp/(omega**2 - omega_c**2)*np.cos(omega*i*tstep)+ charge/mass*e_field_amp/(omega**2 - omega_c**2)
    	x[i][1] = v0/omega_c*np.cos(omega_c*i*tstep) + charge/mass*e_field_amp/(omega**2 - omega_c**2)*np.sin(omega*i*tstep)
    	v[i][0] = v0*np.cos(omega_c*i*tstep) - charge/mass*e_field_amp*omega/(omega**2 - omega_c**2)*np.sin(omega*i*tstep)
    	v[i][1] = v0*np.sin(omega_c*i*tstep) - charge/mass*e_field_amp*omega/(omega**2 - omega_c**2)*np.cos(omega*i*tstep)
    return x, v

def ques2a_b():
    charge = 5
    mass = 10
    tf = 200
    tstep = 0.01
    b_field = np.array([0, 0, 1.0])
    e_amp = np.array([10**5, 0, 0])
    e_freq = 0.2
    v0 = np.array([10**4, 0, 0])
    time = int (tf/tstep)
    ntime = np.linspace(0,tf, time + 1)
    radius = mass*np.sum(v0**2)**0.5/(charge*b_field[2])
    mu_0 = boris_push(charge, mass, np.zeros(3), 0, b_field[2], v0, tstep, time)[2]
    plotfig("ques2a_1.png", "Magnetic moment under no E", "time", "mu", [ntime], [mu_0], ["magnetic moment"])
    mu = boris_push(charge, mass, e_amp, 0, b_field[2], v0, tstep, time)[2]
    plotfig("ques2a_2.png", "Magnetic moment under const E", "time", "mu", [ntime], [mu], ["magnetic moment"])
    boris2 = boris_push(charge, mass, e_amp, e_freq, b_field[2], v0, tstep, time)
    boris_pos2 = boris2[0]

    plotfig("ques2b_pos.png", "Trajectory of particle", "x", "y", [boris_pos2[:, 0]], [boris_pos2[:, 1]], ["trajectory"])
    plotfig("ques2b_vel.png", "Velocity of particle", "time", "v", [ntime, ntime], [boris2[1][:, 0], boris2[1][:,1]], ["vx", "vy"])

def ques2c():
    charge = e
    mass = mp
    tf = 2
    tstep = 0.01
    time = int (tf/tstep)
    e_field = np.array([10, 0, 0])
    n = 1e20
    temp = 1*e/k_boltz
    time = int (tf/tstep)
    ntime = np.linspace(0,tf, time + 1)
    xf, vf  = rk2_q3(charge, mass, e_field, n, temp, tstep, time)
    plotfig("ques2c_pos", "Position of particle", "time", "x", [ntime], [xf], ['position'])
    plotfig("ques2c_vel", "Velocity of particle", "time", "x", [ntime], [vf], ['velocity'], 'lin', 4)

ques1()
ques2a_b()
ques2c()