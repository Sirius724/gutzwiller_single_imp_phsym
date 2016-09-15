!sudo systemctl disable cups-browsed.service

program power
  use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

IMPLICIT NONE

INTEGER, PARAMETER              :: idp=KIND(1.0D0), ids=KIND(1.0)
INTEGER, PARAMETER              :: idpc=kind((1.0d0,0.0d0)),idsc=kind((1.0,0.0))
REAL(KIND=idp), PARAMETER       :: Pi = 3.141592653589793D0
complex(KIND=idpc), parameter	:: I=(0.0d0,1.0d0)
REAL(KIND=idp)			:: alfa,eta,gamma0,V,delta, j_1, J_2, w_a, w_b,t
REAL(KIND=idp)			:: t_max,t_min,delta_t,  eta_min, eta_max, delta_eta, delta_min,delta_max, delta_delta
REAL(KIND=idp)			:: w_min, w_max, delta_w,w, w_a_0, w_b_0, w_a_pi, w_b_pi, tau, w_tot, w_tot_0, w_tot_pi, tau_min,tau_max, delta_tau
integer				:: l,m,n,n_t, n_eta,n_delta,j, n_jumps,icase, n_tau, i_tau
 character(LEN=100)::label
integer :: status = 0, ignored_status
integer :: l_dfti(1)
complex(KIND=idpc), allocatable	:: g0t_v(:), gw_v(:), gvt_v(:),g0w_v(:),gvw_v(:)
REAL(KIND=idp),allocatable:: t_v(:), chia_v(:), chib_v(:), w_v(:), t_jumps(:), chitot_v(:), jt_v(:)
type(DFTI_DESCRIPTOR), POINTER :: hand
 CHARACTER(7) :: a_delta
logical:: j1j2



!intervals
delta_t=0.01
n_t=2000001
!
eta_min=1
eta_max=1
delta_eta=0.1
!
tau_min=0
tau_max=5000
delta_tau=100
!
delta_min=0
delta_max=pi/2
delta_delta=pi/2-1e-10
!std values
alfa=1
eta=1
delta=0
tau=10
!gamma0=eta/2
gamma0=1
V=tan(delta)/pi/gamma0
!
label="test"
j_1=3
j_2=5

!t_max=100
!t_min=-t_max
!n_t=int((t_max-t_min)/delta_t)	!should be even
t_max=(n_t)*delta_t/2
t_min=-t_max
delta_w=2*pi/delta_t/n_t
w_min=-(n_t)*delta_w/2
!w_min=-pi/delta_t
w_max=-w_min


n_eta=int((eta_max-eta_min)/delta_eta)
n_delta=2 !int((delta_max-delta_min)/delta_delta)+1
n_tau=int((tau_max-tau_min)/delta_tau)
n_eta=max(1,n_eta)
n_delta=max(1,n_delta)
n_tau=max(1,n_tau)



write(*, "(a10,10f15.5)"), "t_min=", t_min
write(*, "(a10,10f15.5)"), "t_max=", t_max
write(*, "(a10,10f15.5)"), "delta_t=", delta_t
write(*, "(a10,i10)"), "n_t=", n_t
write(*, "(a10,10f15.5)"), "w_min=", w_min
write(*, "(a10,10f15.5)"), "w_max=", w_max
write(*, "(a10,10f15.5)"), "delta_w=", delta_w
write(*, "(a10,10f15.5)"), "D_w*D_t=", delta_w*delta_t
write(*, "(a10,i10)"), "n_tau=", n_tau
write(*, "(a10,i10)"), "n_eta=", n_eta
write(*, "(a10,i10)"), "n_delta=", n_delta




 call system('mkdir ' //trim(label))
OPEN(unit=11,file=trim(label)//'/w',status='unknown')
OPEN(unit=13,file=trim(label)//'/w0pi',status='unknown')

allocate(gw_v(0:n_t-1))
allocate(g0w_v(0:n_t-1))
allocate(gvw_v(0:n_t-1))
allocate(g0t_v(0:n_t-1))
allocate(gvt_v(0:n_t-1))
allocate(t_v(0:n_t-1))
allocate(jt_v(0:n_t-1))
allocate(w_v(0:n_t-1))
allocate(chia_v(0:n_t-1))
allocate(chib_v(0:n_t-1))
allocate(chitot_v(0:n_t-1))

do j=0, n_eta-1						!eta loop
eta=eta_min+delta_eta*j
gamma0=eta/(2.0d0)
!gamma0=1
!delta=atan(pi*gamma0*V)

do i_tau=0,n_tau-1					!tau loop
tau=tau_min+delta_tau*i_tau


do l=0,n_delta-1					!delta loop
delta=delta_min+l*delta_delta
V=tan(delta)/pi/gamma0

icase=3						!icase loop

!print *, "a", icase,l,i_tau,j

if (icase==1) then
n_jumps=1
!if (.not. allocated(t_jumps)) 
allocate(t_jumps(0:n_jumps-1))
t_jumps(0)=tau
j1j2=.true.	!first j1, then J2
elseif (icase==2) then
n_jumps=1
!if (.not. allocated(t_jumps)) 
allocate(t_jumps(0:n_jumps-1))
t_jumps(0)=0
j1j2=.false.	!first j2, then J1
elseif (icase==3) then
n_jumps=2
!if (.not. allocated(t_jumps)) 
allocate(t_jumps(0:n_jumps-1))
t_jumps(0)=-tau
t_jumps(1)=+tau
j1j2=.true.	!first j1, then J2, then j1
elseif (icase==4) then
n_jumps=2
!if (.not. allocated(t_jumps))  
allocate(t_jumps(0:n_jumps-1))
t_jumps(0)=-tau
t_jumps(1)=+tau
j1j2=.false.	!first j2, then j1, then j2
endif

write (*, "(a10,i10)"), "case=", icase
write (*, "(a10,i10)"), "n_jumps=", n_jumps

!print *, "b", icase,l,i_tau,j


!if (i_tau==0 .and. icase==1 .and. j==0 .and. l==0) write(*,"(20a10)"),"alfa","eta","gamma0", "delta/pi","j_1","j_2","tau"
!write(*,"(10f10.5)"),alfa,eta,gamma0, delta/pi,j_1,j_2,tau


!!!finds gv(t)
!
do n=0,n_t-1
t=t_min+(n)*delta_t
w=w_min+(n)*delta_w
t_v(n)=t
w_v(n)=w
jt_v(n)=j_t(t)
enddo


!print *, "c", icase,l,i_tau,j

WRITE(UNIT=a_delta,FMT='(SPf7.4)') delta/pi
OPEN(unit=10,file=trim(label)//'/chi'//a_delta,status='unknown')
OPEN(unit=12,file=trim(label)//'/gw'//a_delta,status='unknown')

!delta=0
if (abs(delta)<1e-5) then
do n=0,n_t-1
t=t_v(n)
g0t_v(n)=g0t(t)
enddo
gvt_v=g0t_v

!delta=pi/2
elseif (abs(delta-pi/2)<1e-5) then
do n=0,n_t-1
t=t_v(n)
g0t_v(n)=g0t(t)
enddo
gvt_v=0

!general case
else

 call gvt_fft

write(12,"(i10,20a15)"),"#n","w_v(n)", "re g0w_v(n)","im g0w_v(n)", "re gvw_v(n)", "im gvw_v(n)"
do n=0,n_t-1
!t=t_v(n)
!w=w_v(n)
!write(10,"(i10,20e15.5)"),n,t_v(n), g0t(t), gvt_v(n), f_tip(t), chia_v(n), chib_v(n)
write(12,"(i10,20e15.5)"),n,w_v(n), g0w_v(n), gvw_v(n)
enddo

endif !delta 0 pi 

!print *, "d", icase,l,i_tau,j

!chi
 chia_v=0
 chib_v=0
 chitot_v=0
do n=0,(n_t-1)/2
t=t_v(n)
 chia_v(n)=chia(t) !2*real(f_tip(-t)*(gvt_v(minus_t(n))/exp_delta(t,0.0d0)-g0t(t)*exp_delta(t,0.0d0)))	!n_t-1-n gives -t
 chib_v(n)=chib(t) !2*real(f_tip(-t)*(g0t(-t)*exp_delta(t,0.0d0)-gvt_v(n)/exp_delta(t,0.0d0)))
enddo

do n=0,n_t-1
t=t_v(n)
 do m=0,n_jumps-1
if (j1j2) then
  if (mod(m,2)==0)  chitot_v(n)=chitot_v(n)+chi_tot(t_jumps(m),t)
  if (mod(m,2)==1)  chitot_v(n)=chitot_v(n)-chi_tot(t_jumps(m),t)
else
  if (mod(m,2)==1)  chitot_v(n)=chitot_v(n)+chi_tot(t_jumps(m),t)
  if (mod(m,2)==0)  chitot_v(n)=chitot_v(n)-chi_tot(t_jumps(m),t)
endif
 enddo
enddo


!writes gv
write(10,"(a10,20a15)"),"#n","t_v(n)", "jt_v(n)", "re g0t(t)","im g0t(t)", "re gvt_v(n)","im gvt_v(n)", "re f_tip(t)","im f_tip(t)","chia_v(n)", "chib_v(n)", "chitot_v(n)"
do n=0,n_t-1
t=t_v(n)
w=w_v(n)
!write(10,"(10e15.5)"),t, chia(t), chib(t), f_tip(t), g0t(t), gvt(t)
!write(10,"(2i10,10e15.5)"),n,minus_t(n), t_v(n),t_v(minus_t(n)), f_tip(t), g0t(t), gvt_v(n)
write(10,"(i10,20e15.5)"),n,t_v(n), jt_v(n), g0t(t), gvt_v(n), f_tip(t), chia_v(n), chib_v(n), chitot_v(n)
!write(12,"(i10,20e15.5)"),n,w_v(n), g0w_v(n), gvw_v(n)
enddo

!print *, "e", icase,l,i_tau,j


!W
w_a=0
w_b=0
w_tot=0
do n=0,n_t-1
w_a=w_a+chia_v(n)
w_b=w_b+chib_v(n)
w_tot=w_tot+chitot_v(n)*jt_v(n)
enddo
w_a=w_a*delta_t*J_1*(J_2-J_1)
w_b=w_b*delta_t*J_2*(J_1-J_2)
w_tot=w_tot*delta_t*(J_2-J_1)

!print *, "f", icase,l,i_tau,j


!writes w
if (j==0 .and. i_tau==0 .and. l==0) write(11,"(20a15)"),"#alfa","eta","v","gamma0", "delta/pi","tau", "j_1","j_2", "w_a","w_b","w_a+w_b", "w_tot"
write(11,"(20e15.5)"),alfa,eta,v,gamma0, delta/pi,tau,j_1,j_2, w_a,w_b,w_a+w_b, w_tot
!write(*,"(10f10.5)"),alfa,eta,gamma0, delta/pi,j_1,j_2, w_a,w_b, w_tot
!if (i_tau==0 .and. icase==1 .and. j==0 .and. l==0) write(*,"(20a10)"),"alfa","eta","gamma0", "delta/pi","j_1","j_2","tau", "w_a", "w_b", "w_tot"
!write(*,"(10f10.5)"),alfa,eta,gamma0, delta/pi,j_1,j_2,tau, w_a,w_b, w_tot




if (abs(delta)<1e-5) w_a_0=w_a
if (abs(delta)<1e-5) w_b_0=w_b
if (abs(delta)<1e-5) w_tot_0=w_tot
if (abs(delta-pi/2)<1e-5) w_a_pi=w_a
if (abs(delta-pi/2)<1e-5) w_b_pi=w_b
if (abs(delta-pi/2)<1e-5) w_tot_pi=w_tot




 close(unit=10)
 close(unit=12)


!if (allocated(t_jumps))  
deallocate(t_jumps)


if (j==0 .and. i_tau==0 .and. l==0) write(13,"(20a13)"),"#alfa","eta","gamma0", "j_1","j_2", "tau", "w_a_0","w_b_0","w_a_0+w_b_0","w_tot_0", "w_a_pi","w_b_pi","w_a_pi+w_b_pi", "w_tot_pi"
 if (abs(delta-pi/2)<1e-5) write(13,"(20e13.5)"),alfa,eta,gamma0, j_1,j_2, tau, w_a_0,w_b_0,w_a_0+w_b_0,w_tot_0, w_a_pi,w_b_pi,w_a_pi+w_b_pi, w_tot_pi

! print *,"a", icase,l,i_tau,j

!if (i_tau==0 .and. j==0 .and. l==0) write(*,"(12a10)"),"alfa","eta","gamma0", "j_1","j_2", "tau", "w_a_0","w_b_0","w_tot_0", "w_a_pi","w_b_pi", "w_tot_pi"
!if (abs(delta-pi/2)<1e-5) write(*,"(12f10.5)"),alfa,eta,gamma0, j_1,j_2, tau, w_a_0,w_b_0,w_tot_0, w_a_pi,w_b_pi, w_tot_pi


!enddo !icase
enddo	!delta
enddo	!tau
enddo	!eta

 close(unit=11)
 close(unit=13)




contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex function f_tip(t)
real(kind=idp),intent(in):: t

 f_tip=alfa/(2.0d0)/(I*t+alfa)

end function f_tip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=idp) function J_t(t)
real(kind=idp),intent(in):: t
integer::n
!integer,intent(in):: icase


J_t=J_1

if (icase==1 .and. t>t_jumps(0)) J_t=J_2 !case A
if (icase==2 .and. t<t_jumps(0)) J_t=J_2 !case B

if (icase>2) then		!general case
 if (j1j2) then
  if (mod(n_jumps,2)==0) then !i.e. j1j2j1
   do n=0,n_jumps/2-1
    if (t>t_jumps(n) .and. t<t_jumps(n+1)) j_t=J_2
   enddo
  else !i.e. j1j2 j1j2j1j2
   do n=0,n_jumps-1/2-1
    if (t>t_jumps(n) .and. t<t_jumps(n+1)) j_t=J_2
   enddo
   if (t>t_jumps(n_jumps-1)) j_t=j_2
  endif	!mod
 else	!j2j1
J_t=J_2
  if (mod(n_jumps,2)==0) then !i.e. j2j1j2
   do n=0,n_jumps/2-1
    if (t>t_jumps(n) .and. t<t_jumps(n+1)) j_t=J_1
   enddo
  else !i.e. j1j2 j1j2j1j2
   do n=0,n_jumps-1/2-1
    if (t>t_jumps(n) .and. t<t_jumps(n+1)) j_t=J_1
   enddo
   if (t>t_jumps(n_jumps-1)) j_t=j_1
  endif	!mod
 endif
endif



end function J_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

complex function exp_delta(t1,t2)
real(kind=idp),intent(in):: t1,t2

if (abs(delta)<1e-5) then
exp_delta=1
elseif (abs(delta-pi/2)<1e-5) then
 exp_delta=sqrt((abs(t1)-I*eta)/(abs(t2)-I*eta)) !-delta**2/pi**2/2)
else
 exp_delta=((abs(t1)-I*eta)/(abs(t2)-I*eta))**(delta/pi) !-delta**2/pi**2/2)
endif

end function exp_delta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex function G0t(t)
real(kind=idp),intent(in):: t

if (t>0) G0t=-gamma0/(t-i*eta)
if (t<0) G0t=-gamma0/(t+i*eta)
!if (t>0) G0t=-gamma0*(exp(-I*t/eta)-1)/t
!if (t<0) G0t=-gamma0*(exp(+I*t/eta)-1)/t

if (t==0) G0t=0

end function G0t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex function GVt(t)
real(kind=idp),intent(in):: t

!if (delta>1e-5) then
!GVt=G0t(t)*cos(delta)**2-eta/V/pi/(t**2+eta**2)*sin(delta)**2
if (abs(delta<1e-5)) then
GVt=G0t(t)
elseif (abs(delta-pi/2<1e-5)) then
GVt=0
else
print *, " I do not have Gv(t) "
stop
endif

end function GVt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex function G0V_t(t1,t2)
real(kind=idp),intent(in):: t1,t2
!integer,intent(in):: icase

if (abs((J_t(t2+1e-4)-J_2)<1e-5)) G0V_t=GVt(t1-t2)
if (abs((J_t(t2+1e-4)-J_1)<1e-5)) G0V_t=G0t(t1-t2)
!G0V_t=G0t(t)

end function G0V_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex function GT_t(t1,t2)
real(kind=idp),intent(in):: t1,t2
!integer,intent(in):: icase
integer::n

GT_t=G0V_t(t1,t2)
!GT_t=G0V_t(-t2)

if (j1j2) then
do n=0,n_jumps-1
if (mod(n,2)==1) GT_t=GT_t*exp_delta(t1-t_jumps(n),t2-t_jumps(n))
if (mod(n,2)==0) GT_t=GT_t*exp_delta(t2-t_jumps(n),t1-t_jumps(n))
enddo
else
do n=0,n_jumps-1
if (mod(n,2)==0) GT_t=GT_t*exp_delta(t1-t_jumps(n),t2-t_jumps(n))
if (mod(n,2)==1) GT_t=GT_t*exp_delta(t2-t_jumps(n),t1-t_jumps(n))
enddo
endif

end function GT_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function chia(t)				!one jump
real(kind=idp),intent(in):: t

 chia=2*real(f_tip(-t)*(gvt(-t)*exp_delta(0.0d0,t)-g0t(t)*exp_delta(t,0.0d0)))

end function chia
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function chib(t)				!one jump
real(kind=idp),intent(in):: t

 chib=2*real(f_tip(-t)*(g0t(-t)*exp_delta(t,0.0d0)-gvt(t)*exp_delta(0.0d0,t)))

end function chib
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function chi_tot(t1,t2)		!two jumps
real(kind=idp),intent(in):: t1,t2

 chi_tot=0
if (t1>t2)  chi_tot=2*real(f_tip(t1-t2)*(gT_t(t1,t2)-gT_t(t2,t1)))

end function chi_tot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function minus_t(n)
integer,intent(in):: n

if(n==0) then
minus_t=0
else
minus_t=n_t-n
endif

end function minus_t


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gvt_fft

!with corrections
do n=0,n_t-1
g0t_v(n)=g0t(t)*(-1)**n
enddo

gw_v=g0t_v

l_dfti(1)=n_t

  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 1,l_dfti)
 if (0 /= status) print *, "error"
!print *, status

  status = DftiCommitDescriptor(hand)
 if (0 /= status) print *, "error"

  status = DftiComputeBackward(hand, gw_v)	!+		t->w exp(iwt)
 if (0 /= status) print *, "error"

  status = DftiFreeDescriptor(hand)
 if (0 /= status) print *, "error"

!corrections
do n=0,n_t-1
gw_v(n)=gw_v(n)!*exp(i*n*delta_w*t_min)!*exp(i*t_min*w_min)
gw_v(n)=gw_v(n)*(i**n_t)*((-1)**n)
enddo

g0w_v=gw_v

!gv
gw_v=gw_v/(1-gw_v*V)

gvw_v=gw_v

!corrections
do n=0,n_t-1
gw_v(n)=gw_v(n)!*exp(-i*t_min*n*delta_w)
gw_v(n)=gw_v(n)*((-1)**n)
enddo

!fft
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 1,l_dfti)
 if (0 /= status) print *, "error"

  status = DftiCommitDescriptor(hand)
 if (0 /= status) print *, "error"

  status = DftiComputeForward(hand, gw_v)	!-		w->t exp(-iwt)
 if (0 /= status) print *, "error"

  status = DftiFreeDescriptor(hand)
 if (0 /= status) print *, "error"

!corrections
do n=0,n_t-1
gvt_v(n)=gw_v(n)!*exp(-i*w_min*n*delta_t)!*exp(-i*w_min*t_min)
gvt_v(n)=gw_v(n)*((-1)**n)*((-i)**n_t)
gvt_v(n)=gvt_v(n)/n_t
enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!ft by hand
!do n=0,n_t-1
!t=t_min+(n)*delta_t
!w=w_min+(n)*delta_w
!t_v(n)=t
!w_v(n)=w
!g0t_v(n)=g0t(t)
!enddo

!gw_v=0
!do n=0,n_t-1
!do m=0,n_t-1
!t=t_v(n)
!w=w_v(m)
!gw_v(m)=gw_v(m)+g0t_v(n)*exp(i*w*t)
!enddo
!enddo

!gw_v=gw_v/(1-gw_v*V)

!gvt_v=0
!do n=0,n_t-1
!do m=0,n_t-1
!t=t_v(n)
!w=w_v(m)
!gvt_v(n)=gvt_v(n)+gw_v(m)*exp(-i*w*t)
!enddo
!enddo

!gvt_v=gvt_v/n_t
!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine gvt_fft

end program power
