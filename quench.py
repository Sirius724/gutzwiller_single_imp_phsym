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

def occup(E,mu,T):
#	return 1./(np.exp((E-mu)/T)+1)
	if(E-mu >=0):	return np.exp(-(E-mu)/T)/(np.exp(-(E-mu)/T)+1)
	if(E-mu <0):	return 1/(np.exp((E-mu)/T)+1)
occup_v=np.vectorize(occup)

def main():
	Delta=1.
	N=200
	delta=Delta/N
	d=0.5
	V=1
	A=pi*V**2/(2*Delta**2) #/(2*N)
	tau=10
	mu=0
	T=0.001
	n_step=1000

	epsilon_v=np.zeros((2*N))
	alfa_v=np.zeros((2*N))		#energy o phase shift for cond. electrons
	E_v=np.zeros((2*N+1))		# approximate energies
	g_v=np.zeros((2*N+1))
	lambda_v=np.zeros((2*N+1,2*N))
	m_v=range(-N,N+1)
	n_v=range(-N,N)



#Energies
	for n in range(-N,N):
		epsilon_v[n+N]=delta*(n+d)
		alfa_v[n+N]=1/pi*(atan(A*N/(n+d+1e-8)).real)
	for m in range(-N,N):
		E_v[m+N]=epsilon_v[m+N]+delta*alfa_v[m+N]
	E_v[2*N]=0		#bound state at zero energy


#Hamiltonian
	Ham=np.zeros((2*N+1,2*N+1))	#last state is d
	for n in range(-N,N):
		Ham[n+N,n+N]=epsilon_v[n+N]
		Ham[2*N,n+N]=V/np.sqrt(2*N)
		Ham[n+N,2*N]=V/np.sqrt(2*N)


#diagonalizes Hamiltonian

	energies,rot_mat0=np.linalg.eig(Ham)

	#	print(energies)



#computes energies numerically with fsolve satrting with approx solution E_v
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



# Creates matrix of eigenvectors
	rot_mat=np.concatenate((lambda_v/np.sqrt(2.*N),np.transpose(np.asmatrix(g_v))),axis=1)


#Checks that it diagonalizes Hamiltonian
	Ham_diag=np.dot(rot_mat,np.dot(Ham,np.transpose(rot_mat)))
	Ham_diag2=np.dot(np.transpose(rot_mat0),np.dot(Ham,rot_mat0))
#	plt.matshow(Ham)
#	plt.matshow(Ham_diag)
#	plt.matshow(Ham_diag2)
#	plt.show()


#	print(rot_mat-np.transpose(rot_mat0))


#	print(lambda_v.shape)
#	print(rot_mat.shape)

	print("det=",np.linalg.det(rot_mat))
	print("det0=",np.linalg.det(rot_mat0))
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


#Plotting graphic method

	if (N<30):

		plt.plot(x,y,'b',label='g')
		plt.plot(x,x,'r',label='g')
		plt.vlines(epsilon_v,-1000,1000,linewidth=3)
		plt.vlines(E_v,-1000,1000,linewidth=3,linestyles='dashed',colors='g')
		plt.vlines(E_true,-1000,1000,linewidth=2,colors='y')
		plt.vlines(energies,-1000,1000,linewidth=2,linestyles='dotted',colors='c')
		plt.ylim(-1.4,1.4)
		plt.xlim(-1.4,1.4)
	#	plt.legend()
		plt.show()




#Checks that the two methods are consistent
#	print(np.sort(energies))
#	print(np.sort(E_true))
	assert(np.all(abs(np.sort(energies)-np.sort(E_true))<1e-8))

	epsilon_vm=np.zeros((2*N+1))
	for n in range(-N,N):
		epsilon_vm[n+N]=epsilon_v[n+N]
	epsilon_vm[2*N]=0			# I add the d state at zero energy


#############
# H0
	H_00=np.zeros((2*N+1,2*N+1))			#in original basis nd diagonal
	for m in range(-N,N+1):
		H_00[m+N,m+N]=epsilon_vm[m+N]
#	H_00=np.copy(Ham)

#	Ham_diag=np.dot(rot_mat,np.dot(Ham,np.transpose(rot_mat)))

	H_01=np.dot(rot_mat,np.dot(H_00,np.transpose(rot_mat)))	#in  basis m

# H1
	H_11=np.zeros((2*N+1,2*N+1))		#in basis m diagonal
	for m in range(-N,N+1):
		H_11[m+N,m+N]=E_true[m+N]

	assert(np.all(abs(H_11-Ham_diag)<1e-8))

	H_10=np.dot(np.transpose(rot_mat),np.dot(H_11,rot_mat))	#in original basis diagonal, so this is theoriginal Ham
	assert(np.all(abs(Ham-H_10)<1e-8))


#	plt.matshow(H_11-Ham_diag)
#	plt.colorbar()
#	plt.matshow(H_11)
#	plt.colorbar()
#	plt.matshow(Ham_diag)
#	plt.colorbar()
#	plt.show()



#	plt.matshow(H_10)
#	plt.colorbar()
#	plt.show()


#	print(H_10)

##########################Time evolution
#time evolution in H0
	U_00=np.zeros((2*N+1,2*N+1),dtype=complex)
	for m in range(-N,N+1):
		U_00[m+N,m+N]=np.exp(-1j*H_00[m+N,m+N]*tau)
	U_01=np.dot(rot_mat,np.dot(U_00,np.transpose(rot_mat)))
###	Ham_diag=np.dot(rot_mat,np.dot(Ham,np.transpose(rot_mat)))



#time evolution in H1
	U_11=np.zeros((2*N+1,2*N+1),dtype=complex)
	for m in range(-N,N+1):
		U_11[m+N,m+N]=np.exp(-1j*H_11[m+N,m+N]*tau)
	U_10=np.dot(np.transpose(rot_mat),np.dot(U_11,rot_mat))


#Time evolution starting from eigenstates of H0
#H0 H1 H0 H1 H0 ...

	E_step=np.zeros((n_step))
	U_rot=np.eye(2*N+1)
	step=0

	while (step < n_step):
		A=np.dot(U_rot.conj().T,np.dot(H_00,U_rot))
		for n in range(-N,N+1):
			E_step[step]+=A[n+N,n+N].real*occup(epsilon_vm[n+N],mu,T)
		step+=1
		U_rot=np.dot(U_00,U_rot)
		A=np.dot(U_rot.conj().T,np.dot(H_10,U_rot))
		for n in range(-N,N+1):
			E_step[step]+=A[n+N,n+N].real*occup(epsilon_vm[n+N],mu,T)
		step+=1
		U_rot=np.dot(U_10,U_rot)

	E_step=E_step/N
	print(E_step)


#Time evolution starting from eigenstates of H1
#nothing H1 H0 H1 H0 ...
#	n_step=n_step
	E_step1=np.zeros((n_step))
	U_rot=np.eye(2*N+1)
	step=0
	
	while (step < n_step):
		A=np.dot(U_rot.conj().T,np.dot(H_11,U_rot))
		for n in range(-N,0):
			E_step1[step]+=A[n+N,n+N].real*occup(E_true[n+N],mu,T)
		step+=1
		U_rot=np.dot(U_11,U_rot)
		A=np.dot(U_rot.conj().T,np.dot(H_01,U_rot))
		for n in range(-N,0):
			E_step1[step]+=A[n+N,n+N].real*occup(E_true[n+N],mu,T)
		step+=1
		U_rot=np.dot(U_01,U_rot)

	E_step1=E_step1/N
	print(E_step1)


#	x=np.linspace(-1,1,1001)
#	y=occup_v(x,0,0.01)
#	plt.plot(x,y)

	plt.plot(E_step[1:],'b',label='0')	#start from 1 so that they are in the same cycle, starting with H1
	plt.plot(E_step1,'r',label='1')
#	plt.plot(E_step[1:]-E_step[1:],'b',label='0')	#start from 1 so that they are in the same cycle, starting with H1
	plt.legend()
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
