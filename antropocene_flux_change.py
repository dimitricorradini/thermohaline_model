import numpy as np
import matplotlib.pyplot as plt
rate =1/1000
k = 1
alpha = 0.5
T = 30
beta = -1
e = 30

h = 0.01
tmax = 200000

X = np.zeros((2*tmax,))
t = np.zeros((2*tmax,))
X[0]=16.697

def f1(a,t):
  return rate*(-2*a*np.abs(a)+2*k*alpha*T*np.abs(a)-2*k*beta*e)
betas = []
betas.append(0)
#runge kutta of order 4
for i in range(tmax):
  t[i]=h*i
  
  t1 = h*f1(X[i], t[i])

  t2 = h*f1(X[i]+t1/2, t[i]+h/2)

  t3 = h*f1(X[i]+t2/2, t[i]+h/2)

  t4 = h*f1(X[i]+t3, t[i]+h)

  X[i+1]=X[i]+1/6*(t1+2*t2+2*t3+t4)
  beta +=2e-5
  betas.append(beta)

for i in range(tmax, 2*tmax-1):
  t[i]=h*i
  
  t1 = h*f1(X[i], t[i])

  t2 = h*f1(X[i]+t1/2, t[i]+h/2)

  t3 = h*f1(X[i]+t2/2, t[i]+h/2)

  t4 = h*f1(X[i]+t3, t[i]+h)

  X[i+1]=X[i]+1/6*(t1+2*t2+2*t3+t4)
  beta -= 2e-5
  betas.append(beta)


X_plot =[]
for i, x in enumerate(X):
  X_plot.append(x/(k*alpha*T))
  

plt.plot(betas, X_plot)
plt.ylabel('Adimensional flux')
plt.xlabel('Salinity forcing')
plt.show()





