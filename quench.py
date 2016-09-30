import numpy as np
import scipy as sp
import scipy.integrate as spi
import matplotlib.pyplot as plt
import functions as fct
import cmath
from cmath import pi as pi
from cmath import tan as tan
from cmath import atan as atan


def f(x,e):
	g=0
	for _,en in enumerate(e):
		g+=1/(x-en)	
	return g
#f_v= np.vectorize(f)



def epsilon_n(n,d,delta):return delta*(n+d)
epsilon_n_v= np.vectorize(epsilon_n)

def alfa_m(m,d,A,N): return 1/pi*(atan(A*N/(m+d+1e-8)).real)
alfa_m_v= np.vectorize(alfa_m)

def E_m(m,d,delta,A,N): return epsilon_n(m,d,delta)+delta*alfa_m(m,d,A,N)
E_m_v= np.vectorize(E_m)

#def g_m(m,d,delta,A,N,V): return np.sqrt(1/(1+pi**2*V**2*A**2*N/(2*(A**2*delta**2*N**2+epsilon_n(m,d,delta)**2))))
def g_m(m,d,delta,A,N,V): return np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*np.sin(pi*(alfa_m(m,d,A,N)))**2)))
#g_m_v= np.vectorize(g_m)

def lambda_mn(m,n,d,delta,A,N,V): return g_m(m,d,delta,A,N,V)*V/(E_m(m,d,delta,A,N)-epsilon_n(n,d,delta))
#lambda_mn_v= np.vectorize(lambda_mn)

def main():
	Delta=1.
	N=5
	delta=Delta/N
	d=0.5
	V=1
	A=pi*V**2/(2*Delta**2)

	epsilon_v=np.zeros((2*N))
	alfa_v=np.zeros((2*N))
	E_v=np.zeros((2*N))
	g_v=np.zeros((2*N))
	lambda_v=np.zeros((2*N,2*N))
	m_v=range(-N,N)
	n_v=range(-N,N)

#	print(lambda_v)

#	epsilon_v=epsilon_n_v(n_v,d,delta)
#	alfa_v=alfa_m_v(m_v,d,A,N)
#	E_v=E_m_v(m_v,d,delta,A,N)
#	g_v=g_m_v(m_v,d,delta,A,N,V)
	#lambda_v=lambda_mn_v(m_v,n_v,d,delta,A,N,V)


#	for m in range(-N,N):
#		g_v[m+N]=g_m(m,d,delta,A,N,V)
#		epsilon_v[m+N]=epsilon_n(m,d,delta)
#		alfa_v[m+N]=alfa_m(m,d,A,N)
#		E_v[m+N]=E_m(m,d,delta,A,N)
#		for n in range(-N,N):
#			lambda_v[m+N,n+N]=lambda_mn(m,n,d,delta,A,N,V)
#			lambda_v[m+N,n+N]=

	for m in range(-N,N):
		epsilon_v[m+N]=delta*(m+d)
		alfa_v[m+N]=1/pi*(atan(A*N/(m+d+1e-8)).real)
		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*np.sin(pi*(alfa_v[m+N]))**2)))
#		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*N**2*A**2/(N**2*A**2+(m+d+alfa_v[m+N])**2))))	#exact
#		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*N**2*A**2/(N**2*A**2+(m+d)**2))))		#approx
		E_v[m+N]=epsilon_v[m+N]+delta*alfa_v[m+N]
		for n in range(-N,N):
			lambda_v[m+N,n+N]=g_v[m+N]*V/(E_v[m+N]-epsilon_v[n+N])



#	print(lambda_v)

#	plt.plot(m_v,epsilon_v,'r',label='epsilon')
#	plt.plot(m_v,E_v,'b',label='E')
#	plt.plot(m_v,alfa_v,'b',label='alfa')
#	plt.plot(m_v,g_v,'b',label='g')
#	plt.legend()
#	plt.show()


	print(np.sum(g_v**2))
	print(np.sum(lambda_v**2,axis=1))

#starting with H_0 and quenching to H_1
	E_old=0
	for n in range(-N,0):	
		E_old+=epsilon_v[n+N]
	E_old=E_old/N
	print(E_old)
	E_new=0
	for n in range(-N,0):	
		for m in range(-N,N):
			E_new+=E_v[n+N]*lambda_v[m+N,n+N]**2
	E_new=E_new/N

	print(E_old,E_new)


#	E_new=sum(epsilon_v[epsilon_v<0])/N


#starting with H_1 and quenching to H_0
	E_old1=0
	for m in range(-N,0):	
		E_old1+=E_v[m+N]
	E_old1=E_old1/N

	E_new1=0
	for n in range(-N,N):	
		for m in range(-N,0):
			E_new1+=epsilon_v[n+N]*lambda_v[m+N,n+N]**2
	E_new1=E_new1/N


	print(E_old1,E_new1)



	x=np.linspace(-1.1,1.1,1000)
	y=V**2/(2*N*Delta)*f(x,epsilon_v)


	plt.plot(x,y,'b',label='g')
	plt.plot(x,x,'r',label='g')
	plt.vlines(epsilon_v,-1000,1000,linewidth=4)
	plt.vlines(E_v,-1000,1000,linewidth=3,linestyles='dashed',colors='g')
	plt.ylim(-1.2,1.2)
	plt.xlim(-1.1,1.1)
	plt.legend()
	plt.show()







if __name__ == '__main__':
	main()
