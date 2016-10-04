import numpy as np
import scipy as sp
import scipy.integrate as spi
import scipy.optimize as spo
import matplotlib.pyplot as plt
import functions as fct
import cmath
from cmath import pi as pi
from cmath import tan as tan
from cmath import atan as atan


def f(x,e):
	"""	sum_n 1/(x-e[n])	"""
	g=0
	for _,en in enumerate(e):
		g+=1/(x-en)	
	return g

def f2(x,e,pend):
	"""	sum_n 1/(x-e[n])-pend*x	"""
	g=f(x,e)-pend*x
	return g


#def epsilon_n(n,d,delta):return delta*(n+d)
#epsilon_n_v= np.vectorize(epsilon_n)

#def alfa_m(m,d,A,N): return 1/pi*(atan(A*N/(m+d+1e-8)).real)
#alfa_m_v= np.vectorize(alfa_m)

#def E_m(m,d,delta,A,N): return epsilon_n(m,d,delta)+delta*alfa_m(m,d,A,N)
#E_m_v= np.vectorize(E_m)

#def g_m(m,d,delta,A,N,V): return np.sqrt(1/(1+pi**2*V**2*A**2*N/(2*(A**2*delta**2*N**2+epsilon_n(m,d,delta)**2))))
#def g_m(m,d,delta,A,N,V): return np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*np.sin(pi*(alfa_m(m,d,A,N)))**2)))
#g_m_v= np.vectorize(g_m)

#def lambda_mn(m,n,d,delta,A,N,V): return g_m(m,d,delta,A,N,V)*V/(E_m(m,d,delta,A,N)-epsilon_n(n,d,delta))
#lambda_mn_v= np.vectorize(lambda_mn)

def main():
	Delta=1.
	N=50
	delta=Delta/N
	d=0.5
	V=1
	A=pi*V**2/(2*Delta**2) #/(2*N)
	tau=0.4578


	epsilon_v=np.zeros((2*N))
	alfa_v=np.zeros((2*N))		#energy o phase shift for cond. electrons
	E_v=np.zeros((2*N+1))		# approximate energies
	g_v=np.zeros((2*N+1))
	lambda_v=np.zeros((2*N+1,2*N))
	m_v=range(-N,N+1)
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

	for n in range(-N,N):
		epsilon_v[n+N]=delta*(n+d)
		alfa_v[n+N]=1/pi*(atan(A*N/(n+d+1e-8)).real)
#		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*N**2*A**2/(N**2*A**2+(m+d+alfa_v[m+N])**2))))	#exact
#		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*N**2*A**2/(N**2*A**2+(m+d)**2))))		#approx
	for m in range(-N,N):
		E_v[m+N]=epsilon_v[m+N]+delta*alfa_v[m+N]
	E_v[2*N]=0		#bound state at zero energy

#computes energies numerically with fsolve
	E_true=np.zeros(2*N+1)
	for i,ei in enumerate(E_v):
		E_true[i]=spo.fsolve(f2,E_v[i],args=(epsilon_v,2*N*Delta/V**2))

#	alfa_true=(E_true-epsilon_v)/delta


#	print (E_true)

	for m in range(-N,N+1):
		sum_n=0.
		for n in range(-N,N):
			sum_n+=1./(E_true[m+N]-epsilon_v[n+N])**2
		g_v[m+N]=np.sqrt(1./(1.+(V**2)/(2.*N)*sum_n))

	for m in range(-N,N+1):
#		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*np.sin(pi*(alfa_v[m+N]))**2)))
#		g_v[m+N]=np.sqrt(1/(1+pi**2*V**2/(2*N*delta**2*np.sin(pi*(alfa_true[m+N]))**2)))
		for n in range(-N,N):
#			lambda_v[m+N,n+N]=g_v[m+N]*V/(E_v[m+N]-epsilon_v[n+N])
			lambda_v[m+N,n+N]=g_v[m+N]*V/(E_true[m+N]-epsilon_v[n+N])	#/np.sqrt(2*N)



#	print(np.transpose(np.asmatrix(g_v)).shape)
	rot_mat=np.concatenate((lambda_v/np.sqrt(2.*N),np.transpose(np.asmatrix(g_v))),axis=1)


#	print(lambda_v.shape)
#	print(rot_mat.shape)

	print("det=",np.linalg.det(rot_mat))
#	print(lambda_v)

#	plt.plot(m_v,epsilon_v,'r',label='epsilon')
#	plt.plot(m_v,E_v,'b',label='E')
#	plt.plot(m_v,alfa_v,'b',label='alfa')
#	plt.plot(m_v,g_v,'b',label='g')
#	plt.legend()
#	plt.show()

#	print (g_v**2)
#Check normalization
#	print(np.sum(g_v**2))
#	print(np.sum(lambda_v**2,axis=1)/(2.*N)+g_v**2)		#sum over n
#	print(np.sum(lambda_v**2,axis=0)/(2.*N))		#sum over m


#######################
#starting with H_0 and quenching to H_1
#initial energy <psi_0IH_0Ipsi_0>
	E_old=0
	for n in range(-N,0):		#occupied states
		E_old+=epsilon_v[n+N]
	E_old=E_old/N			#E per electron

#<psi_0 I H_1 I psi_0>
	E_new=0
	for n in range(-N,0):	
		for m in range(-N,N+1):
			E_new+=E_true[m+N]*lambda_v[m+N,n+N]**2/(2.*N)
	E_new=E_new/N

	print(E_old,E_new)
	print((E_old-E_new)*N)


#	E_new=sum(epsilon_v[epsilon_v<0])/N




#starting with H_1 and quenching to H_0
#<psi_1 I H_1 I psi_1>
	E_old1=0
	for m in range(-N,0):	
		E_old1+=E_true[m+N]
	E_old1=E_old1/N

#<psi_1 I H_0 I psi_1>
	E_new1=0
	for n in range(-N,N):	
		for m in range(-N,0):
			E_new1+=epsilon_v[n+N]*lambda_v[m+N,n+N]**2/(2.*N)
	E_new1=E_new1/N


	print(E_old1,E_new1)
	print((E_old1-E_new1)*N)


	x=np.linspace(-1.4,1.4,1000)
	y=V**2/(2*N*Delta)*f(x,epsilon_v)


#Plotting

	if (n<30):

		plt.plot(x,y,'b',label='g')
		plt.plot(x,x,'r',label='g')
		plt.vlines(epsilon_v,-1000,1000,linewidth=4)
		plt.vlines(E_v,-1000,1000,linewidth=3,linestyles='dashed',colors='g')
		plt.vlines(E_true,-1000,1000,linewidth=2,colors='y')
		plt.ylim(-1.4,1.4)
		plt.xlim(-1.4,1.4)
#	plt.legend()
		plt.show()



	epsilon_vm=np.zeros((2*N+1))
	for n in range(-N,N):
		epsilon_vm[n+N]=epsilon_v[n+N]
	epsilon_vm[2*N]=0			# I add the d state at zero energy


###########Energies at the third interval H0 H1 H0

	U_mm=np.zeros((2*N+1,2*N+1),dtype=complex)
	for m in range(-N,N+1):
		U_mm[m+N,m+N]=np.exp(-1j*E_v[n+N]*tau)


	U_npn=np.dot(np.transpose(rot_mat),np.dot(U_mm,rot_mat))

#	print(U_mm)



	E_new2=0.
	for m in range(-N,N+1):
		for n in range(-N,0):
			E_new2+=np.absolute(U_npn[m+N,n+N])**2*epsilon_vm[m+N]
	E_new2=E_new2/N


	print(E_old,E_new,E_new2)	#H_0,H_1,H_0


###########Energies at the third interval H1 H0 H1 

	U_mm2=np.zeros((2*N+1,2*N+1),dtype=complex)
	for m in range(-N,N+1):
		U_mm2[m+N,m+N]=np.exp(-1j*epsilon_v[n+N]*tau)


	U_npn2=np.dot(rot_mat,np.dot(U_mm,np.transpose(rot_mat)))

#	print(U_mm)



	E_new3=0.
	for m in range(-N,0):
		for m2 in range(-N,N+1):
			E_new3+=np.absolute(U_npn2[m2+N,m+N])**2*E_v[m2+N]
	E_new3=E_new3/N


	print(E_old1,E_new1,E_new3)	#H_1,H_0,H_1


###########Energies at the fourth interval H0 H1 H0 H1
#All in basis n of H0

# H0
	H_00=np.zeros((2*N+1,2*N+1))
	for m in range(-N,N+1):
		H_00[m+N,m+N]=epsilon_vm[m+N]

	H_01=np.dot(rot_mat,np.dot(H_00,np.transpose(rot_mat)))

# H1
	H_11=np.zeros((2*N+1,2*N+1))
	for m in range(-N,N+1):
		H_11[m+N,m+N]=E_v[m+N]
	H_10=np.dot(np.transpose(rot_mat),np.dot(H_11,rot_mat))

#time evolution in H0
	U_00=np.zeros((2*N+1,2*N+1),dtype=complex)
	for m in range(-N,N+1):
		U_00[m+N,m+N]=np.exp(-1j*epsilon_vm[m+N]*tau)
	U_01=np.dot(rot_mat,np.dot(U_00,np.transpose(rot_mat)))

	

#time evolution in H1
	U_11=np.zeros((2*N+1,2*N+1),dtype=complex)
	for m in range(-N,N+1):
		U_11[m+N,m+N]=np.exp(-1j*E_v[m+N]*tau)
	U_10=np.dot(np.transpose(rot_mat),np.dot(U_11,rot_mat))

#	U_20=np.dot(U_10,U_00)

	n_step=50
	E_step=np.zeros((n_step))
	U_rot=np.eye(2*N+1)
	step=0

	while (step < n_step):
		A=np.dot(U_rot.conj().T,np.dot(H_00,U_rot))
		for n in range(-N,0):
			E_step[step]+=A[n+N,n+N]
		step+=1
		U_rot=np.dot(U_00,U_rot)
		A=np.dot(U_rot.conj().T,np.dot(H_10,U_rot))
		for n in range(-N,0):
			E_step[step]+=A[n+N,n+N]
		step+=1
		U_rot=np.dot(U_10,U_rot)

	E_step=E_step/N
	print(E_step)


	plt.plot(E_step,'b',label='g')
#	plt.legend()
	plt.show()


#	E_new4=0
#	for n in range(-N,N+1):
#		for n2 in range(-N,0):
#			E_new4+=np.absolute(U_20[n+N,n2+N])**2*epsilon_vm[n+N]
#	E_new4=E_new4/N
#	print(E_old,E_new,E_new2,E_new4)	#H_0,H_1,H_0



#	print(U_mm)




if __name__ == '__main__':
	main()
