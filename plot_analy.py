import numpy as np 
import matplotlib.pylab as plt

ma = np.logspace(-12,-5,100) 

t = 1e-9
m0 = ma[0]
a = 1.0
N0 = 1.0/m0


g = np.exp(-a*t)
a_sol = N0/m0 * (2.0/(a*N0*t))**2 * np.exp((2.0/(a*N0*t)) * (1.0 - ma/m0))
#a_sol = N0/(2.0*np.sqrt(np.pi)*m0**2)*(g/(1.0-g)**0.75) * np.exp(-(ma/m0)*(1.0 - np.sqrt(1.0 - g))**2)

t = 2e-9
g = np.exp(-a*t)
a_sol_2 = N0/m0 * (2.0/(a*N0*t))**2 * np.exp((2.0/(a*N0*t)) * (1.0 - ma/m0))
#a_sol_2 = N0/(2.0*np.sqrt(np.pi)*m0**2)*(g/(1.0-g)**0.75) * np.exp(-(ma/m0)*(1.0 - np.sqrt(1.0 - g))**2)


fname = 'test.txt'

data = np.loadtxt(fname)

q_init = data[0,:]
q = data[1,:]
m = data[2,:]

fig = plt.figure()

plt.plot(ma,a_sol*ma**2)
plt.plot(ma,a_sol_2*ma**2)
plt.plot(m,q_init*m**2)
plt.plot(m,q*m**2)

plt.ylim(1e-6,1e3)

plt.yscale('log')
plt.xscale('log')

plt.show()