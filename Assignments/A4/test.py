import numpy as np
import matplotlib.pyplot as plt

k_boltz = 1.38064e-23
e = 1.60271e-19
eps = 8.85419e-12
mp = 1.67262e-27

def entropyf(v, pdf_speed):
    entropy = 0
    for i in range(len(v)):
        if(pdf_speed[i] != 0.0):
            entropy += -1*pdf_speed[i]*np.log(pdf_speed[i])*(v[i] - v[i-1])
            #print entropy, vel[i], pdf_speed[i]
    return entropy

def testf(mean, sigma):
    vel = np.linspace(mean-6*sigma, mean + 6*sigma, 1000)
    f = 1/(2*np.pi*sigma**2)**0.5*np.exp(-(vel-mean)**2/(2*sigma**2))
    entropy = entropyf(vel, f)
    print entropy
    print np.log(sigma*(2*np.pi*np.e)**0.5)

testf(0, (1667*e/mp)**0.5)