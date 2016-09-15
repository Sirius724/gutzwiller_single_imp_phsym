import numpy as np
import scipy as sp
import scipy.integrate as spi
import matplotlib.pyplot as plt
import functions as fct
import cmath
from cmath import pi as pi
from cmath import tan as tan

#
def fft(a,dt):
	n_l=len(a)
	b=np.zeros(n_l,dtype=complex)
	for n in range(n_l):
		b[n]=a[n]*(-1)**n	
	c= np.fft.fft(b)
	for n in range(n_l):
		c[n]*=((-1)**n)*dt	#(1j**n_l)*
	return c
#
def ifft(a,dw):
	n_l=len(a)
	b=np.zeros(n_l,dtype=complex)
	for n in range(n_l):
		b[n]=a[n]*(-1)**n	
	c= np.fft.ifft(b)
	for n in range(n_l):
		c[n]*=((-1)**n)*dw		#((-1j)**n_l)
	return c

#
def J_t(t,J_1,J_2):
	if (t>0):
		J_t=J_2 
	elif (t<0):
		J_t=J_1
	else:
		J_t=0
	return J_t
J_t_v = np.vectorize(J_t)

#
def G0_t(t,gamma0,eta):
	if (t>0):
		G0_t=-gamma0/(t-eta*1j)
	elif (t<0):
		G0_t=-gamma0/(t+eta*1j)
	else:
		G0_t=0
	return G0_t
G0_t_v = np.vectorize(G0_t)

#
def GV_t_v(g0_v,delta,dt,V):
	g0w_v=np.zeros(n_t,dtype=complex)	#G_0(omega)
	gvw_v=np.zeros(n_t,dtype=complex)	#G_V(omega)
	if (abs(delta)<1e-5):
		gv_v=g0_v
	elif (abs(delta-pi/2)<1e-5):
		gv_v=0
	else:
		g0w_v=fft(g0_v,dt)
		gvw_v=g0w_v/(1-V*g0w_v)
		gv_v =ifft(gvw_v,1/dt)
	return gv_v,g0w_v,gvw_v
#
def f_tip(t,alfa):
	return alfa/(2.)/(1j*t+alfa)

#
#def g0V_t(t1,t2):
#	if (abs((J_t(t2+1e-4)-J_2)<1e-5)):
#		G0V_t=GVt(t1-t2)
#	if (abs((J_t(t2+1e-4)-J_1)<1e-5)):
#		G0V_t=G0t(t1-t2)
#	return G0V_t
#
def exp_delta(t1,t2,delta,eta):
	if (abs(delta)<1e-5):
		exp_delta=1
	elif (abs(delta-pi/2)<1e-5):
		exp_delta=np.sqrt((abs(t1)-1j*eta)/(abs(t2)-1j*eta)) 
	else:
	 	exp_delta=((abs(t1)-1j*eta)/(abs(t2)-1j*eta))**(delta/pi) 
	return exp_delta

# chi_jj(t1,t2) where t1=0
def Chi_J_t(t,gamma0,eta):
	if (t>=0):
		Chi_t=0
	if (t<0):
		Chi_t=2*(f_tip(-t,alfa)(-t,gamma0,eta)*G0_t(t,gamma0,eta)).imag
	return Chi_t
Chi_J_t_v = np.vectorize(Chi_J_t)

#chi_rhorho(t1,t2) where t1=0
def Chi_V_t(t,gamma0,eta):
	if (t>=0):
		Chi_t=0
	if (t<0):
		Chi_t=2*(G0_t(-t,gamma0,eta)*G0_t(t,gamma0,eta)).imag
	return Chi_t
Chi_V_t_v = np.vectorize(Chi_V_t)

# index of -t
def minus_t(n,n_t):
	if(n==0):
		minus_t=0
	else:
		minus_t=n_t-n
	return minus_t

#########
def main():

	if (icase==1):
		n_jumps=1
		t_jumps[0]=0	

	n_t=2**18		#points for t dicretization
	t_m=np.sqrt(n_t*pi/2)	#sqrt(n_t*pi/2) optimal choice for best accuracy when eta=1
	t_v = np.linspace(-t_m,t_m,n_t)
	dt=t_v[1]-t_v[0]	#discretization step
	j_v = J_t_v(t_v,J_1,J_2)	#current as a fct of time

	w_v=np.zeros(n_eta)
	eta_v=np.linspace(eta_min,eta_max,n_eta)

	for index,eta in enumerate(eta_v):
		gamma0=eta/2.
		V=(tan(delta)/pi/gamma0).real
		print eta,gamma0,delta,V
		g0_v = G0_t_v(t_v,gamma0,eta)
		#chi_t_v = Chi_t_v(t_v,gamma0,eta)
		#w_v[index]=spi.quad(Chi_t, -np.inf, 0,args=(gamma0,eta))[0]

		#Computes G_V(t)
		(gv_v,g0w_v,gvw_v)=GV_t_v(g0_v,delta,dt,V)

#	print t_jumps
#	plt.plot(t_v,j_v)
#	plt.plot(t_v,g0_v.real)
#	plt.plot(t_v,g0_v.imag)
#		plt.plot(t_v,chi_t_v.real)
#	plt.plot(t_v,chi_t_v.imag)
#	print w_v

	plt.plot(t_v,g0w_v.real,label='g0w real')
	plt.plot(t_v,g0w_v.imag,label='g0w imag')
	plt.plot(t_v,gvw_v.real,label='gvw real')
	plt.plot(t_v,gvw_v.imag,label='gvw imag')

	plt.legend()
	plt.show()


	plt.plot(t_v,g0_v.real,label='g0 real')
	plt.plot(t_v,g0_v.imag,label='g0 imag')
	plt.plot(t_v,gv_v.real,label='gV real')
	plt.plot(t_v,gv_v.imag,label='gV imag')



#	plt.plot(eta_v,w_v)

	plt.legend()
	plt.show()




if __name__ == '__main__':
	#intervals
	delta_t=0.01
	n_t=2000001
	#
	eta_min,eta_max,delta_eta,n_eta=(1,1,0.1,1)
	tau_min,tau_max,delta_tau,n_tau=(0,5000,100,100)
	delta_min,delta_max,delta_delta,n_delta=(0,pi/2,pi/2-1e-10,2)
	#std values
	alfa=1.
#	eta=1.
	delta=0.45*pi
	tau=10
#	gamma0=1
#	gamma0=eta/2.
#	V=(tan(delta)/pi/gamma0).real
	#
	label="test"
	#
	icase=1
	J_1=1
	J_2=3
	n_jumps=1
	t_jumps=np.zeros(n_jumps)

	main()



