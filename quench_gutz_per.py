#f2py -L/path/to/lapack -llapack -m module -c module.f

import scipy.linalg.blas as blas
import scipy.linalg.lapack as lapack
import numpy as np
import scipy as sp
import matplotlib
from matplotlib.ticker import ScalarFormatter
import scipy.integrate as spi
import scipy.optimize as spo
from scipy.optimize import fsolve
from scipy.optimize import minimize,minimize_scalar
import matplotlib.pyplot as plt
#import functions as fct
import cmath
from cmath import pi as pi
from cmath import tan as tan
from cmath import atan as atan
from numpy import sqrt as sqrt
from numpy import log as log
import sys
from power import fft
#import routines
from numpy import sin,cos
from scipy.integrate import odeint
from scipy.integrate import ode
from scipy.integrate import complex_ode
from odeintw import odeintw
#from odelab import Solver

np.set_printoptions(linewidth=500,precision=6)

def prim_(x,a,label):
	if (label=="sqrt"): return 1+x*sqrt(1-x**2)/pi-np.arccos(x)/pi-a
#	if (label=="arcos"): return 1-sqrt(1-x**2)/pi+x*np.arccos(x)/pi-a
	if (label=="cos"): return 0.5+0.5*sin(pi*x/2)-a

def set_energies(label,N):
	energies=np.zeros(2*N,dtype=float)
	for n in np.arange(2*N):
		integral=float((n+0.5)/(2*N))
#		print(n,integral)
		if (label=="sqrt"): energies[n]=float(fsolve(prim_, 0,(integral,"sqrt"))[0])
#		if (label=="arcos"): energies[n]=float(fsolve(prim_, 0,(integral,"arcos"))[0])
		if (label=="cos"): energies[n]=float(fsolve(prim_, 0,(integral,"cos"))[0])
		if (label=="theta"): energies[n]=(n+0.5-N)/(N)
	return energies				#in units Delta

def f_to_str(x,a=15,b=10):
	''' float to string'''
	str_f="{:+0"+str(a)+"."+str(b)+"f}"
	return str(str_f.format(x))
#	return str("{:+06.3f}".format(x))

def occup(E,mu,T):
#	return 1./(np.exp((E-mu)/T)+1)
	if(E-mu >=0):	return np.exp(-(E-mu)/T)/(np.exp(-(E-mu)/T)+1)
	if(E-mu <0):	return 1/(np.exp((E-mu)/T)+1)
occup_v=np.vectorize(occup)

def R(P,theta=0):
	if (P<1 and P>0): R=2*sqrt(P*(1-P))*cos(theta)
	else:R=1e-8
	return R
R_v=np.vectorize(R)

#def dR_dP(P,theta=0):
#	if (P>1e-10 and P<1-1e-6): return (1-2*P)/sqrt(P*(1-P))*cos(theta)
#	else: return 1e3*cos(theta)

def dR_dP(P,theta=0):
	return (1-2*P)/sqrt(P*(1-P))*cos(theta)


def dR_dtheta(P,theta):
	return -2*sqrt(P*(1-P))*sin(theta)

def V_av_formula(Gamma,P):
	z=max(4*P*(1-P),1e-10)
	return 4*Gamma*sqrt(z)/pi*log(Gamma*z)

def E_tot_formula(P,Gamma,U,theta=0):		#OK
	if (P<0 or P>1):
		return 1e8
#	z=max(4*P*(1-P)*cos(theta)**2,1e-10)
	else:
		z=4*P*(1-P)*cos(theta)**2
		return 2*Gamma*z/pi*log(Gamma*z/np.e)+U/2*P

def E_tot_numeric(P,U,H0,V_0,T,N):
		if (P<0 or P>1):
			return 1e8
		else:
			R_=R(P)
			Ham=H0+V_0*R_
			energies,rot_mat0,info=lapack.dsyev(Ham)
			occup_m=occup_v(energies,0,T)
			E_tot=0
			for m in range(-N,N+1):		#eigenvalues
				E_tot+=energies[m+N]*occup_m[m+N]
			E_tot*=2		#spin

			return E_tot+U/2*P

def E_tot_V_numeric(P,U,H0,V_0,T,N):
		R_=R(P)
		Ham=H0+V_0*R_
		energies,rot_mat0,info=lapack.dsyev(Ham)
		occup_m=occup_v(energies,0,T)
		E_tot=0
		V_av=0
		T_av=0
		for m in range(-N,N+1):		#eigenvalues
			E_tot+=energies[m+N]*occup_m[m+N]
			for n in range(-N,N):		#eigenvalues
				V_av+=rot_mat0[n+N,m+N]*rot_mat0[2*N,m+N]*occup_m[m+N]
		E_tot*=2		#spin
		V_av*=2*V_0[2*N,0]	#*R_
		V_av*=2			#spin

		return (E_tot+U/2*P,V_av)


def error_formula(P,Gamma,U):
	return U/2+dR_dP(P)*V_av_formula(Gamma,P)

def error_formula_z(z,Gamma,U):
#	return sqrt(1-z)*(log(z*Gamma)+1)+pi*U/(16*Gamma)
	return (log(z*Gamma)+1)+pi*U/(16*Gamma)

def z_P(P,theta=0):
	return 4*P*(1-P)*cos(theta)**2
z_P_v=np.vectorize(z_P)

def P_z(z):
	return 0.5-0.5*sqrt(1-z)

def dz_dP(P,theta=0):
	return (4-8*P)*cos(theta)**2

def dz_dtheta(P,theta):
	return 4*P*(1-P)*(-2*sin(theta)*cos(theta))

def y0(psi,P,theta,N):
	y=np.zeros((2*N+1)**2+2,dtype=complex)	#just create it
	for m in range(0,2*N+1):		#columnwise
		for n in range(0,2*N+1):
			y[n+(2*N+1)*m]=psi[n,m]+0j
	y[(2*N+1)**2]=P+0j
	y[(2*N+1)**2+1]=theta+0j
	return y



#def dy_dt(y,t0,H0,V_0,U,occup,epsilon,N,V,epsilon_vm):	#for odeintw
def dy_dt(t0,y,H0,V_0,U,occup,epsilon,N,V,epsilon_vm):	#for ode
	#unpacks
	psi_old,P_old,theta_old=psi_vec_to_mat(y,N)	#psi is complex

	#computes V_av, needed for dp and dtheta
	T_av,V_av,occup_old,occupd_old=E_T_V_n(psi_old,occup,N,V,epsilon_vm)
	R_old=R(P_old,theta_old)
	E_old=R_old*V_av+T_av+U/2*P_old

#	print (E_old)
	H_t=H0+R_old*V_0
	iHpsi=-1j*blas.zgemm(alpha=1.0, a=H_t, b=psi_old, trans_b=False)
	dP=-V_av*dR_dtheta(P_old,theta_old)
	dtheta=U/2+V_av*dR_dP(P_old,theta_old)
#	print(dP,dtheta)
	#packs into vector
	dy_dt=psi_mat_to_vec(iHpsi,dP,dtheta,N)
	return dy_dt

def dy_dt_real(t0,y,H0,V_0,U,occup,epsilon,N,V,epsilon_vm):
	dy_complex=dy_dt(t0,y_real_to_complex(y,N),H0,V_0,U,occup,epsilon,N,V,epsilon_vm)
	dy_real=y_complex_to_real(dy_complex,N)
#	dy_complex2=y_real_to_complex(dy_real,N)
#	error=np.amax(abs((dy_complex-dy_complex2)))
#	if(error>1e-6): print("dy_dt real error", error)
	return dy_real

def y_real_to_complex(y_real,N):
	M=(2*N+1)**2
	y_complex=np.zeros(M+2,dtype=complex)
	for m in range(0,M):
		y_complex[m]=y_real[m]+1j*y_real[m+M]
	y_complex[M]=y_real[2*M]
	y_complex[M+1]=y_real[2*M+1]
	return y_complex

def y_complex_to_real(y_complex,N):
	M=(2*N+1)**2
	y_real=np.zeros(2*M+2,dtype=float)
	for m in range(0,M):
		y_real[m]=np.real(y_complex[m])
		y_real[m+M]=np.imag(y_complex[m])
	y_real[2*M]=np.real(y_complex[M])
	y_real[2*M+1]=np.real(y_complex[M+1])
	return y_real

def E_T_V_n(psi,occup,N,V,epsilon_vm):
	"""returns  <T>, <V> and <n>"""
	V_av,T_av,n_occ,nd_occ=(0,0,0,0)
	for m in range(-N,N+1):		#eigenvalues
		for n in range(-N,N):		##n conduction states
			V_av+=2*np.real(np.conj(psi[n+N,m+N])*psi[2*N,m+N])*occup[m+N]
			T_av+=epsilon_vm[n+N]*(np.abs(psi[n+N,m+N]))**2*occup[m+N]
		for n in range(-N,N+1):			#n conduction states + d
			n_occ+=(np.abs(psi[n+N,m+N]))**2*occup[m+N]
		nd_occ+=(np.abs(psi[2*N,m+N]))**2*occup[m+N]
	V_av*=V/sqrt(2.*N)
	V_av*=2		#spin
	n_occ*=2	#spin
	T_av*=2
	nd_occ*=2	#spin
	return(T_av,V_av,n_occ,nd_occ)

def psi_vec_to_mat(psi_vec,N):
	M=2*N+1
	psi_mat=np.zeros((M,M),dtype=complex)	#just create it
	for m in range(0,M):
		for n in range(0,M):
			psi_mat[n,m]=psi_vec[n+M*m]
	P=psi_vec[M**2].real
	theta=psi_vec[M**2+1].real
	return(psi_mat,P,theta)

def psi_mat_to_vec(psi_mat,P,theta,N):
	M=2*N+1
	psi_vec=np.zeros((M*M+2),dtype=complex)	#just create it
	for m in range(0,M):
		for n in range(0,M):
			psi_vec[n+M*m]=psi_mat[n,m]
	psi_vec[M**2]=P
	psi_vec[M**2+1]=theta
	return(psi_vec)



def main():
	delta_t=0.01
	t_0=0
	l_deltat=True
	Delta=1.
	N=100
	d=0.5
	Gamma=0.1
	tau=10
	mu=0
	T=0.0001
	U=1
	n_step=1000
	periodic=True
#	lH0=0		#computes <H_0> or <H_1>
#	epsilon_d=0.7
	lfft=True
	l_edinh0=False		#here not used
	dos="theta"
	epsilon_d_0=0	#if l_edinh0=False, the d state must nevertheless exist with an epsilon_d_0 energy>>0
	Gamma0=0
	#reads parameters
	if(len(sys.argv)>1):
		Gamma=float(sys.argv[1])
		U=float(sys.argv[2])
		N=int(sys.argv[3])
		n_step=int(sys.argv[4])
		P_init=float(sys.argv[5])
		T=float(sys.argv[6])
		tau=float(sys.argv[7])
	if(len(sys.argv)>7):
		if(int(sys.argv[8])==1): periodic=True
		if(int(sys.argv[8])==0): periodic=False
	if(len(sys.argv)>9):
		if(int(sys.argv[9])==1): l_deltat=True
		if(int(sys.argv[9])==0): l_deltat=False
		delta_t=float(sys.argv[10])
	if(len(sys.argv)>10):
		epsilon_d_0=float(sys.argv[11])
	if(len(sys.argv)>11):
		dos=sys.argv[12]
	if(len(sys.argv)>12):
		Gamma0=float(sys.argv[13])

	print("Gamma0=",Gamma0)
	if(n_step==0): step_av=0
	if(n_step==0): tau=0

	print(sys.argv[8],periodic)

	step_period=0
	if (delta_t>0): step_period=2*int(tau/delta_t)	#number of steps per period, is even


	delta=Delta/N
	step_av=int(n_step/2) #100 #int(n_step/2)
	V=np.sqrt(2*Gamma/np.pi)
	V0=np.sqrt(2*Gamma0/np.pi)
	#Gamma=np.pi*V**2/2
	print ("V=",V)

#	print(l_edinh0)
#	if(periodic): lfft=False
	print(lfft)

#	if (not l_deltat or not periodic): delta_t=tau
	if (not l_deltat): delta_t=tau

	print("Setting up")
	print ("N=",N)
	print ("ldelta_t=",l_deltat)
	print ("delta_t=",delta_t)
	print ("Tau=",tau)
	print ("Periodic=",periodic)
	if (periodic and l_deltat): print ("steps per period=", step_period)
#Energies
#	epsilon_v=np.zeros((2*N))
#	for n in range(-N,N):
#		epsilon_v[n+N]=delta*(n+d)

	print ("Setting energies")
	epsilon_v=set_energies(dos,int(N))
#	print(epsilon_v)

	epsilon_vm=np.zeros((2*N+1))
	epsilon_all=np.zeros((2*N+1))
	for n in range(-N,N):
		epsilon_vm[n+N]=epsilon_v[n+N]
		epsilon_all[n+N]=epsilon_v[n+N]
#	if(l_edinh0): epsilon_vm[2*N]=epsilon_d #0			# I add the d state at zero energy to be used in H0
#	if(not l_edinh0): epsilon_vm[2*N]=epsilon_d_0 #0			# I add the d state at zero energy to be used in H0
#	epsilon_all[2*N]=epsilon_d					# I add the d state at e_d energy for phase shifts


#H0
	H0=np.zeros((2*N+1,2*N+1))	#last state is d, with zero energy
	for n in range(-N,N+1):
		H0[n+N,n+N]=epsilon_vm[n+N]
#V
	V_0=np.zeros((2*N+1,2*N+1))	#last state is d, with zero energy
	for n in range(-N,N):
		V_0[2*N,n+N]=V/np.sqrt(2.*N)
		V_0[n+N,2*N]=V/np.sqrt(2.*N)
#V0
	V0_0=np.zeros((2*N+1,2*N+1))	#last state is d, with zero energy
	for n in range(-N,N):
		V0_0[2*N,n+N]=V0/np.sqrt(2.*N)
		V0_0[n+N,2*N]=V0/np.sqrt(2.*N)

#	if(not l_edinh0): V_0[2*N,2*N]=(epsilon_d-epsilon_d_0)	# if epsilon_d is not in H_0, I put it in V, minus the epsilon_d_0 energy, so that H has epsilon_d

#	Ham=H0+V_0

#commutator [T,V]
#	CTV2=blas.dgemm(alpha=1.0, a=H0, b=V_0, trans_b=False)-blas.dgemm(alpha=1.0, a=V_0, b=H0, trans_b=False)
	CTV=np.zeros((2*N+1,2*N+1))
	for n in range(-N,N):
		CTV[2*N,n+N]=-V/np.sqrt(2.*N)*epsilon_vm[n+N]
		CTV[n+N,2*N]=V/np.sqrt(2.*N)*epsilon_vm[n+N]

#	print(CTV)
#	print(CTV2)
#	print (np.amax(CTV-CTV2))
#	print (np.amin(CTV-CTV2))
#	print(Ham[2*N,2*N])
#	print (H0)
#	print (V_0)
#	print (Ham)

#diagonalizes Hamiltonian

#	print dir(routines)
#	print routines.diagonalize.__doc__



	print("SCF")
#dE/dP=0
#(1-2P)/2P(1-P) <RV>+U/2=0

	P_old=P_init

#for gamma
#	result=minimize(E_tot_formula,P_init,args=(Gamma,U),method='Nelder-Mead')
	result=minimize(E_tot_numeric,P_old,args=(U,H0,V_0,T,N),method='Nelder-Mead')
	P_new=result.x[0]
	E_tot=result.fun
	R_new=R(P_new)
	z=z_P(P_new)
	print(Gamma,U,P_new,E_tot,R_new,z)
	P_old=P_new
	if not(result.success):
		print(result.message)
		sys.exit()

#for gamma0
#	result0=minimize(E_tot_formula,P_init,args=(Gamma0,U),method='Nelder-Mead')
	result0=minimize(E_tot_numeric,0.1,args=(U,H0,V0_0,T,N),method='Nelder-Mead')
	P_new0=result0.x[0]
	E_tot0=result0.fun
	R_new0=R(P_new0)
	z0=z_P(P_new0)
	print(Gamma0,U,P_new0,E_tot0,R_new0,z0)
	P_old0=P_new0
#	if (P_new0<1e-8):P_old0=0
	if not(result0.success):
		print(result0.message)
		sys.exit()



##################################
	Ham=H0+V_0*R_new
	energies,rot_mat0,info=lapack.dsyev(Ham)
	occup_m=occup_v(energies,mu,T)
	E_tot=0
	for m in range(-N,N+1):
		E_tot+=energies[m+N]*occup_m[m+N]
	E_tot*=2 #spin
#eigenvalues and vectors for the steady state with V and U
	rot_mat_inf=rot_mat0
	energies_inf=energies
	occup_inf=occup_m
	P_inf=P_new

#eigenvalues and vectors for the state with V0 and U
	Ham0=H0+V0_0*R_new0
	energies0,rot_mat00,info0=lapack.dsyev(Ham0)
	occup_m0=occup_v(energies0,mu,T)
	E_tot0=0
	for m in range(-N,N+1):
		E_tot0+=energies0[m+N]*occup_m0[m+N]
	E_tot0*=2 #spin


#non-interacting problem U=0, finite V (R=1)
	Ham=H0+V_0
	energies_U0,rot_mat_U0,info=lapack.dsyev(Ham)
	occup_U0=occup_v(energies_U0,mu,T)
	E_tot_U0=0
	for m in range(-N,N+1):		#eigenvalues
		E_tot_U0+=energies_U0[m+N]*occup_U0[m+N]
	E_tot_U0*=2

#	z_old=z_P(P_old)
#	for U in np.arange(0,4,0.01):
#		z_new=fsolve(error_formula_z,z_old,args=(Gamma,U))[0]
#		P_new=P_z(z_new)
#		E_tot=E_tot_formula(Gamma,P_new,U)
#		R_new=R(P_new)
#		z=R_new**2
#		print(Gamma,U,P_new,E_tot,R_new,z)
#		P_old=P_new

#####Print to file
	N_str=str(N)
	T_str=str(T)
	V_str=f_to_str(V,6,3)	#str("{:+06.3f}".format(V))
	Gamma_str=f_to_str(Gamma,6,3)	#str("{:+06.3f}".format(V))
	Gamma0_str=f_to_str(Gamma0,6,3)	#str("{:+06.3f}".format(V))
	tau_str=str(tau)
	per_str=str(periodic)
	step_str=str(n_step)
	stepav_str=str(step_av)
	U_str=f_to_str(U,6,3)	#str("{:+06.3f}".format(epsilon_d))
	mu_str=f_to_str(mu,6,3)	#str("{:+06.3f}".format(mu))
	if(l_edinh0): H_str='_edH0'
	if(not l_edinh0): H_str='_edH1'
	name="N"+N_str+"_G"+Gamma_str+"_0G"+Gamma0_str+"_T"+T_str+"_tau"+tau_str+"_step"+step_str+"_per"+per_str+"_U"+U_str+H_str+"_dt"+str(l_deltat)+f_to_str(delta_t,6,3)+dos
	file_out=name+".txt"
	file_e=name+"_e.txt"

######## print static quantities
	out_file = open(file_e,"w",1)
	out_file.write(str(P_new)+"\t"+str(E_tot)+"\t"+str(z)+"\n")
	out_file.close()

	P_file = open("P_new","w",1)
	P_file.write(str(P_new))
	P_file.close()


#occupation for psi_0
	occup_n0_1d=occup_v(epsilon_vm,mu,T)
	occup_n0=np.reshape(occup_n0_1d,(2*N+1,1))	#makes this a matrix

	occup_n0_1d_inf=occup_v(energies_inf,mu,T)
	occup_n0_inf=np.reshape(occup_n0_1d_inf,(2*N+1,1))	#makes this a matrix


###initial energy <psi_0 | H_0 | psi_0>
	E_000=0
	for n in range(-N,N+1):		
		E_000+=epsilon_vm[n+N]*occup_n0[n+N,0]
	E_000*=2

#with the eigenvalues for the steady state
	E_000_inf=0
	for n in range(-N,N+1):		
		E_000_inf+=energies_inf[n+N]*occup_n0_inf[n+N,0]
	E_000_inf*=2

#with the eigenvalues for V0
	E_0000_inf=0
	for n in range(-N,N+1):		
		E_0000_inf+=energies0[n+N]*occup_m0[n+N]
	E_0000_inf*=2


	print ("E_000(V=0)=", E_000)
	print ("E_000_inf=", E_000_inf)
	print ("E_tot_U0(U=0)=", E_tot_U0)
	print ("E_0000(V=V0)=", E_0000_inf)



	file_out=name+".txt"
	file_out2=name+"_2.txt"
	file_out3=name+"_3.txt"
	file_power=name+"_power.txt"
	file_power2=name+"_power2.txt"
	file_fft=name+"_fft.txt"
	file_fft_half=name+"_fft_half.txt"
	file_fft_power=name+"_power_fft.txt"
	file_e=name+"_e.txt"



################################################
###Choice of the initial state
#eigenstates at V=0 (U=0 or finite is the same)
#	psi_old=np.eye(2*N+1,dtype=complex)	#use eigenstates of the not hybridized system  -- each column is |m(t)>
#	occup=occup_m
#	P_old=1e-16
#	energies_old=epsilon_vm
#eigenstates at finite V and U
#	psi_old=rot_mat_inf.astype(complex)	#use eigenstates with V and U
#	occup=occup_n0_1d_inf
#	P_old=P_inf
#	energies_old=energies_inf
#eigenstates at U=0 (interaction quench)
#	psi_old=rot_mat_U0.astype(complex)	#use eigenstates with V, but U=0
#	occup=occup_U0
#	P_old=0.5
#	energies_old=energies_U0
#eigenstates at V0 and U
	psi_old=rot_mat00.astype(complex)	#use eigenstates with V0 and U
	occup=occup_m0
	P_old=P_old0
	energies_old=energies0

#Other values
	theta_old=0
	H_t=np.zeros((2*N+1,2*N+1))
	psi_new=np.zeros((2*N+1,2*N+1),dtype=complex)	#just create it
	R_old=R(P_old,theta_old)

##Energy at t=0+
#from eigenvalues
	E_init=0
	occup_init=0
	for m in range(-N,N+1):		
		E_init+=energies_old[m+N]*occup[m+N]
		occup_init+=occup[m+N]
	E_init*=2
	occup_init*=2
	E_init+=U/2*P_old

	T_av,V_av,occup_old,occupd_old=E_T_V_n(psi_old,occup,N,V,epsilon_vm)
	E_old=R_old*V_av+T_av+U/2*P_old

	print ("P_init,theta_init,R_init=",P_old,theta_old,R(P_old,theta_old))
	print ("E_init=",E_init,E_old,E_old-E_init)
	print ("T_init=",T_av, "V_init=", V_av)
	print ("Occup_init=",occup_init,occup_old)


#############
	if (n_step>0):

		t_step=np.zeros((n_step))
		for step in range (0,n_step):
			t_step[step]=step*delta_t+t_0

		print ("y_init")
		y_init=psi_mat_to_vec(psi_old,P_old,theta_old,N)

		#check
		psi_old2,P_old2,theta_old2=psi_vec_to_mat(y_init,N)
		print(P_old-P_old2,abs(psi_old-psi_old2).max(),theta_old-theta_old2)

		t_init=0
		y_init_real=y_complex_to_real(y_init,N)

		out_file = open(file_out3,"w",1)
		out_file_power = open(file_power,"w")

		print ("ode")


#(1)odeintw
#########this works, good E conservation, but has memory issues (stores all results)
#		"""
###		y_t,result=odeint(dy_dt,y_init,t_step,args=(N,occup,H0,V_0,U))
###		y_t, infodict = odeintw(zfunc, z0, t, args=(K, L, M), Dfun=zjac,full_output=True)
#		y_t=odeintw(dy_dt,y_init,t_step,args=(H0,V_0,U,occup,energies_old,N,V,epsilon_vm))
#		"""


#(2)ode(zvode)
#########problems with energy conservation,  good memory (one step at a time)
#(3)complex_ode 
#########ValueError: operands could not be broadcast together with shapes (51,) (50,) (problem with passing arguments to dy_dt)
#(4)ode with real y

		
#(2)
#		r = ode(dy_dt).set_integrator('zvode', method='Adams')
#		r.set_initial_value(y_init, t_init)
#(3)		
#		r = complex_ode(dy_dt).set_integrator('vode', method='Adams')
#		r.set_initial_value(y_init, t_init)
#(4)
		#def dy_dt_real(t0,y,N,occup,H0,V_0,U,epsilon_vm,V):
		r = ode(dy_dt_real).set_integrator('dopri5')		#lsoda, vode: bad energy cons., dopri5, dop853
		r.set_initial_value(y_init_real,t_init)

		r.set_f_params(H0,V_0,U,occup,energies_old,N,V,epsilon_vm)
		t1=n_step*delta_t
		dt=delta_t

		V_external=np.zeros(n_step+1,dtype=int)
		constant_deltat=int(tau/delta_t)		#number of deltat's  during which V is constant
		print(tau,delta_t,constant_deltat)

		i=0
		while i<n_step+1:
			for j in range(0,constant_deltat):
				if (i<n_step): V_external[i]=1
				i+=1
			i+=constant_deltat


		i=0
		while r.successful() and r.t < t1:
#			print(r.t)

			y_t=y_real_to_complex(r.y,N)
			y_mat,P_t,theta_t=psi_vec_to_mat(y_t,N)

			R_t=R(P_t,theta_t)
			z_t=z_P(P_t,theta_t)
			if (V_external[i]==0): T_t,V_t,n_t,nd_t=E_T_V_n(y_mat,occup,N,V0,epsilon_vm)
			if (V_external[i]==1): T_t,V_t,n_t,nd_t=E_T_V_n(y_mat,occup,N,V,epsilon_vm)
			E_t=T_t+R_t*V_t+U/2*P_t
			print(r.t,P_t,theta_t,E_t,n_t,nd_t,T_t,V_t,V_external[i])
			out_file.write(f_to_str(r.t)+"\t"+f_to_str(P_t)+"\t"+f_to_str(theta_t)+"\t"+f_to_str(R_t)+"\t"+\
					f_to_str(E_t)+"\t"+f_to_str(V_t)+"\t"+f_to_str(T_t)+"\t"+f_to_str(n_t)+"\t"+f_to_str(nd_t)+"\t"+f_to_str(z_t)+"\t"+f_to_str(V_external[i])+"\n")

			if (V_external[i]==0): r.set_f_params(H0,V0_0,U,occup,energies_old,N,V0,epsilon_vm)
			if (V_external[i]==1): r.set_f_params(H0,V_0, U,occup,energies_old,N,V, epsilon_vm)
			r.integrate(r.t+dt)
			i=i+1




		print ("end ode")



if __name__ == '__main__':
	main()
