import numpy as np
import matplotlib.pyplot as plt

k_boltz = 1.38064e-23
e = 1.60271e-19
eps = 8.85419e-12
mp = 1.67262e-27

class beam:
    def __init__(self, velocity, density, temp, mol_wt): #temp in ev
        self.den = density
        self.vel = velocity
        self.temp = temp*e/k_boltz
        self.wt = mol_wt*mp

    def getbeam(self):
        return self.den, self.vel, self.temp, self.wt

    def vel_dist(self, vel):
        return self.den*(self.wt/(2*np.pi*k_boltz*self.temp))**0.5*np.exp(-self.wt*(vel - self.vel)**2/(2*k_boltz*self.temp))

    def v_thermal(self):
        return (k_boltz*self.temp/self.wt)**0.5

    def col_freq(self):
        return self.den*e**4/(16*np.pi*eps**2*self.wt**2*self.v_thermal()**3)

def delf(f, dx):
    delf = np.zeros(len(f))
    delf[0] = (f[1] - f[0])/dx
    for i in range(len(f)-1):
        delf[i] = (f[i+1] - f[i-1])/(2*dx)
    delf[len(f)-1] = (f[-1] - f[-2])/dx
    return delf

def entropyf(v, pdf_speed):
    entropy = 0
    for i in range(len(v)):
        if(pdf_speed[i] != 0.0):
            entropy += -1*pdf_speed[i]*np.log(pdf_speed[i])*(v[i] - v[i-1])
            #print entropy, vel[i], pdf_speed[i]
    return entropy

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

def ques1(beams):
    freqmax = 0
    for i in range(len(beams)):
        freqmax = max(freqmax, beams[i].col_freq())
    return freqmax

def ques2(beams):
    minvel = beams[0].getbeam()[1] - 10*beams[0].v_thermal()
    maxvel = beams[0].getbeam()[1] + 10*beams[0].v_thermal()
    for i in range(len(beams)):
        minvel = min(minvel, beams[i].getbeam()[1] - 10*beams[i].v_thermal())
        maxvel = max(maxvel, beams[i].getbeam()[1] + 10*beams[i].v_thermal())
    
    vel = np.linspace(4*minvel, 4*maxvel, 10000)
    f_dist = np.zeros(len(vel))
    for i in range(len(beams)):
        f_dist += beams[i].vel_dist(vel)
    plotfig("Init_dist.png", "Initial Distribution Function in velocity space", "v(m/s)", "n(v)", [vel], [f_dist], ["combined pdf"], "lin", 1)
    return f_dist, vel

def ques3(f_dist, vel, wt):
    dv = vel[1] - vel[0]
    den = np.sum(f_dist*dv)
    mean = np.sum(vel*f_dist*dv)/den
    sigma_sq = np.sum((vel - mean)**2*f_dist*dv)/den
    temp = sigma_sq*wt/e*mp
    return den, mean, temp

def ques4(beam, vel):
    plotfig("Eqb_pdf.png", "Equiibrium Distribution function in velocity space", "v(m/s)", "n(v)", [vel], [beam.vel_dist(vel)], ["equilibrium pdf"])
    return beam.vel_dist(vel)

def ques5(f_in, fmean, dv, col_freq, mass, ext_force):
    delt = 0.10/col_freq
    f_evol = np.zeros((1, len(f_in)))
    f_evol[0, :] = f_in.copy()
    while (np.sum((f_evol[-1] - fmean)**2)**0.5/np.sum(fmean) >= 0.01):
        f = f_evol[-1].copy()
        bgk_rhs = col_freq*(fmean - f)
        extf_term = ext_force/(mass*mp)*delf(fmean, dv)
        f += (bgk_rhs + ext_force)*delt
        f_evol = np.vstack([f_evol, f.copy()])
    print f_evol.shape[1]*delt
    return f_evol

def ques6(f, vel, den, delt):
    entropy = np.zeros(len(f[:, 0]))
    for i in range(len(f)):
        entropy[i] = entropyf(vel, f[i]/den)
    time = np.arange(0, f.shape[0])*delt
    plotfig("Entropy.png", "Variation of entropy with time", "time", "entropy", [time], [entropy], ["Entropy"])

if __name__== "__main__":
    beamhot = beam(0, 10**10, 1, 1)
    v_therm = beamhot.v_thermal()
    beamcold1 = beam(50*v_therm, 10**10, 0.01, 1)
    beamcold2 = beam(-50*v_therm, 10**10, 0.01, 1)
    col_freq = ques1([beamhot, beamcold1, beamcold2])
    f_dist, vel = ques2([beamhot, beamcold1, beamcold2])
    den, mean, temp = ques3(f_dist, vel, 1)
    print den, mean, temp
    beam_eqb = beam(mean, den, temp, 1)
    fmean = ques4(beam_eqb, vel)
    f_evol = ques5(f_dist, fmean, vel[1] - vel[0], col_freq, 1, 0)
    ques6(f_evol, vel, den, 0.1/col_freq)
