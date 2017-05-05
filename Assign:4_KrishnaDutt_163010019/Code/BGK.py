### Python script to implement BGK model for a 1D plasma ###
import pylab as plt
from scipy.integrate import simps
from copy import deepcopy

def PDF(m,KT,N,V,Vo) :
    """ Method to implement Boltzmann distribution """
    
    F = (N*(m/(2*plt.pi*KT))**0.5)*plt.exp(-(m/(2*KT))*(V-Vo)**2)
    #F = ((m/(2*plt.pi*KT))**0.5)*plt.exp(-(m/(2*KT))*(V-Vo)**2)
    return F
    

def PlotPDF(m) :
    """ Method to plot PDF of ions """
    
    N = 10.0**10
    KT = 1.0  # in ev , so corresponding value of normalised m is considered
    Velocity = plt.linspace(-100,100,1001)
    #Warm Plasma
    FW = PDF(m,KT,N,Velocity,0.0)
    
    #Beam 1
    Vth = 50.0*(KT/m)**0.5
    FB1 = PDF(m,KT/100.0,N,Velocity,Vth)
    
    #Beam 2
    FB2 = PDF(m,KT/100.0,N,Velocity,-Vth)
    
    plt.figure()
    plt.grid()
    plt.plot(Velocity,FW,'b--')
    plt.plot(Velocity,FB1,'r--')
    plt.plot(Velocity,FB2,'r--')
    plt.legend(['Warm Plasma','Cold Beam 1','Cold Beam 2'], loc=0)
    plt.title("Distribution Function of ions in velocity space")
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("PDF")
    plt.show()
    

def Integrate(X,Y) :
          """ Method to find and return integral of Entropy using Simpson's method """
          
          length = len(X) # Should be odd for simpson's to work 
          if length%2 == 0 :
             length = length -1
          I = 0.0
          for i in range(0,len(Y)) :
              if Y[i] == 0 :
                 pass
              else :
                 Y[i] = Y[i]*plt.log(Y[i])
                    
          for i in range(1,length,2) :
              """ Find coefficient and integral """
              
              if X[i+1] == 0 :
                 #print "Zero detected"
                 c = Y[i+1]
                 a = (Y[i]*X[i-1] -Y[i-1]*X[i])/(X[i-1]*X[i]*(X[i]-X[i-1]))
                 b = (Y[i-1] -Y[i+1] - a*X[i-1]**2)/X[i-1]
              else :
                a = ((X[i]-X[i-1])*Y[i+1] - X[i+1]*(Y[i]-Y[i-1]))/( (X[i]- X[i-1])*X[i+1]*(X[i+1] - X[i] - X[i-1]))
                b = ( Y[i+1] - a*X[i+1]**2 )/X[i+1]
                c = Y[i-1] - a*X[i-1]**2 -b*X[i-1]
              #print i,X[i+1],(X[i]- X[i-1]),(X[i+1] - X[i] - X[i-1])
              I += (a/3.0)*( X[i+1]**3 - X[i-1]**3 ) + (b/2.0)*(X[i+1]**2 - X[i-1]**2) + c*(X[i+1]-X[i-1])
                
          return I



def RK4(Fn,h,Feqb,nuC) :
    """ Method to implement RK4 """
    
    K1 = nuC*(Feqb-Fn)
    K2 = nuC*(Feqb -(Fn + K1*h*0.5 ))
    K3 = nuC*(Feqb -(Fn + K2*h*0.5 ))
    K4 = nuC*(Feqb -(Fn + K3*h))
    
    F = Fn + (K1+2*K2+2*K3+K4)*(h/6.0)
    return F
        

def BGK(m) :
    """ Method to implement BGK model """ 
    
    N = 10.0**10
    KT = 1.0  # in ev , so corresponding value of normalised m is considered
    Velocity = plt.linspace(-80,80,1000)
    #Warm Plasma
    FW = PDF(m,KT,N,Velocity,0.0)
    #Beam 1
    Vth = 50.0*(KT/m)**0.5
    FB1 = PDF(m,KT/100.0,N,Velocity,Vth)
    #Beam 2
    FB2 = PDF(m,KT/100.0,N,Velocity,-Vth)
    
    # Entropy
    S = []
    
    # Numerically computing Density and temperature 
    NDensity = simps(FW,Velocity) + simps(FB1,Velocity) + simps(FB2,Velocity)
    Energy = simps(0.5*m*Velocity**2*FW,Velocity) + simps(0.5*m*(Velocity-Vth)**2*FB1,Velocity) + simps(0.5*m*(Velocity+Vth)**2*FB2,Velocity)
    #Energy = simps(0.5*m*Velocity**2*FW,Velocity) + simps(0.5*m*(Velocity)**2*FB1,Velocity) + simps(0.5*m*(Velocity)**2*FB2,Velocity)
    Temperature = 2*Energy/(NDensity)
    VelAvg = (simps(Velocity*FW,Velocity) + simps(Velocity*FB1,Velocity) + simps(Velocity*FB2,Velocity))/NDensity
    
    print NDensity , Temperature, VelAvg
    
    
    # Equilibrium Boltzmann Distribution
    Feqb = PDF(m,Temperature,NDensity,Velocity,VelAvg)
    
    # BGK Model 
    nuC = (NDensity/(m**0.5*Temperature**1.5))*2.89*10**-27
    print nuC
    
    Error = 1.0
    F = FW+FB1+FB2
    
    h = 0.1 # Normalised wrt to nuc
    time = 0.0
    
    # Finding Entropy
    S.append(abs(Integrate(Velocity,deepcopy(FW)) + Integrate(Velocity,deepcopy(FB1)) + Integrate(Velocity,deepcopy(FB2))))
    
    while (Error > 0.01) :
          F = RK4(F,h,Feqb,1.0)
          Error =   (sum(((F - Feqb))**2))**0.5/(sum(Feqb**2))**0.5
          time = time + h
          # Find entropy 
          S.append(abs(Integrate(Velocity,deepcopy(F)))) 
          
    print "Time required to reach equilibrium : %f"%(time)
    
    """
    plt.figure()
    plt.grid()
    plt.plot(Velocity,Feqb,'b.--')
    plt.plot(Velocity,F,'r--')
    plt.plot(Velocity,FW,'g--')
    plt.plot(Velocity,FB1,'c--')
    plt.plot(Velocity,FB2,'c--')
    plt.show()
    """
    # Entropy Plot
    Time = plt.linspace(0.0,time,len(S))
    plt.figure()
    plt.grid()
    plt.plot(Time,S,'g--')
    plt.title("Evolution of Entropy with time")
    plt.xlabel('Time ( t*(collision frequency))')
    plt.ylabel("Entropy")
    plt.show()
    


""" Test Case """
if __name__ == "__main__" :
   mass = 1.0
   PlotPDF(mass)
   BGK(mass)       
