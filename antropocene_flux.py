import numpy as np
import matplotlib.pyplot as plt
rate =1/1000
k = 1
alpha = 0.5
T = 30
beta = 0.5
e = 30

h = 0.01
tmax = 500000

X = np.zeros((2*tmax,))
t = np.zeros((2*tmax,))
X[0]=20

def f1(a,t):
	return rate*(-2*a*np.abs(a)+2*k*alpha*T*np.abs(a)-2*k*beta*e)

#runge kutta of order 4
for i in range(tmax):
  t[i]=h*i
  
  t1 = h*f1(X[i], t[i])

  t2 = h*f1(X[i]+t1/2, t[i]+h/2)

  t3 = h*f1(X[i]+t2/2, t[i]+h/2)

  t4 = h*f1(X[i]+t3, t[i]+h)

  X[i+1]=X[i]+1/6*(t1+2*t2+2*t3+t4)
  beta +=1e-5

for i in range(tmax, 2*tmax-1):
  t[i]=h*i
  
  t1 = h*f1(X[i], t[i])

  t2 = h*f1(X[i]+t1/2, t[i]+h/2)

  t3 = h*f1(X[i]+t2/2, t[i]+h/2)

  t4 = h*f1(X[i]+t3, t[i]+h)

  X[i+1]=X[i]+1/6*(t1+2*t2+2*t3+t4)
  beta -=2e-5


X_plot = []

for i, x in enumerate(X):
  if i % 10 == 0:
    X_plot.append(x/(k*alpha*T))
  

    

plt.plot( np.arange(len(X_plot))*h*10, X_plot)
plt.ylabel('Adimensional flux')
plt.xlabel('Time')
plt.show()


