C    program modelE4_1.f    by Bartlomiej Borek , Sept. 2013

C    **this is basic simulation of model E4 in 1D

      implicit none
      integer, parameter :: nx=261
      integer, parameter :: ntime_end=50000000
      double precision, parameter :: pi=3.14159265358979323846
      double precision L(0:nx+1), H(0:nx+1), A(0:nx+1), P(0:nx+1)
      double precision Lt(0:nx+1), Ht(0:nx+1), At(0:nx+1), Pt(0:nx+1)
      double precision lapL, lapH, dL, dH, dA, dP
      double precision rspace, rtime, r1, pert
      double precision D,dx,dt
      double precision a1,a2,g1,g2,a3,g3
      double precision d1,d2
      double precision c1,c2,c3,c4,c5
      double precision Pm,k1,k2
      double precision f1,f2,f3
      integer i1,i2,i,ii,iq,ntime,fntime

C--------> set other parameters
C%% dt=1.0d-6 doesn't give very different fixed point.
C%% D is ratio of the 2 diffusion coefficients, so here space is rescaled.
      data D,dx,dt/
     &    1.0d2,5.0d-1,1.0d-5/
      data a1,a2,d1,d2,g1,g2/
     &    4.5d0,1.6667d0,0.004d0,0.00202d0,3.0d0,2.0d0/
      data c1,c2,c3,c4/
     &    0.675d0,0.675d0,1.0d-2,0.675d0/
      data a3,g3,c5/
     &    3.2d0,1.1d0,0.675d0/
      data Pm, k1,k2/
     &    1.0d0,1.0d2,2.4815d2/

C--------> open files to write to
      open(unit=90,file='L_E4_9_261.txt', status='unknown')
      open(unit=91,file='H_E4_9_261.txt', status='unknown')
      open(unit=92,file='A_E4_9_261.txt', status='unknown')
      open(unit=93,file='P_E4_9_261.txt', status='unknown')

C--------> set initial conditions
         do ii=0,nx+1
	      rspace=(ii-1)*dx
c           L(ii)=0.005445572791697402d0
c            H(ii)=0.005194499371405287d0
C % the stable fixed point of the homogenous system is numerically found for this set of parameters (see modelE1_simplest1_4.nb code). For the dx and dt chosen here the convergence is to L=0.005213230849236,H=0.005080131526230
C % for 3V model A=0.00521323084940534

c------> random noise w.r.t. space IC
           call random_number(r1)
c		  2*Pi*rspace/(dx*(nx-1))
		  pert=0.001*SIN(rspace*0.1/2*pi) 	
c		  write(99,*) pert, rspace, pi
            L(ii)=0.0099537412903138415697630222462534+pert
            H(ii)=0.0073825689719981407350966878264083
            A(ii)=0.019670567928752395134518190011574
            P(ii)=0.0039951539638696398637446589771259
c     +0.1+0.0002*(r1-0.5)
c           L(ii)=0.0099537412903138415697630222462534
         enddo

c         L(100)=0.05
c         r(135:165)=0.5


C************ Start main iteration loop *************
         do ntime=1,ntime_end
            fntime=ntime-1
            rtime=(fntime-1.0)*dt

C----->write spacetime files
            if(mod(fntime,100000)==0)then
c               write(99,*) fntime
               write(90,*) L(1:nx)
               write(91,*) H(1:nx)
               write(92,*) A(1:nx)
               write(93,*) P(1:nx)
            endif
         
C----------> integrate the equations (forward Euler and first order Laplacian)
  
            do i=1,nx
       
              lapL=(L(i+1)-L(i))+(L(i-1)-L(i))
c               lapL=0.0             
	         lapH=(H(i+1)-H(i))+(H(i-1)-H(i))
c               lapH=0.0
	       
	         f1=a1*(d1+P(i))/(P(i)+c1)
	         f2=g1*A(i)*(L(i)/(L(i)+c3))
		    f3=k1*(Pm-P(i))*L(i)-k2*P(i)

	         dL=f1-f2-f3
	         dH=a2*(d2+P(i))/(P(i)+c4)-g2*H(i)
	         dA=(a3*H(i))/(H(i)+c5)-g3*A(i)
	         dP=f3

              Lt(i)=L(i)+dt*dL+dt*lapL/(dx**2)
              Ht(i)=H(i)+dt*dH+dt*D*lapH/(dx**2)
              At(i)=A(i)+dt*dA
              Pt(i)=P(i)+dt*dP

            enddo 

C---------> Periodic boundary conditions for L,H,A,P

            Lt(nx+1)=Lt(1)
            Lt(0)=Lt(nx)
            Ht(nx+1)=Ht(1)
            Ht(0)=Ht(nx)
            At(nx+1)=At(1)
            At(0)=At(nx)
            Pt(nx+1)=Pt(1)
            Pt(0)=Pt(nx)   

C---------> update state variables

            do i2=0,nx+1
               L(i2)=Lt(i2)
               H(i2)=Ht(i2)
               A(i2)=At(i2)
               P(i2)=Pt(i2)
            enddo

C******************* end time iteration loop *************
         enddo
C*****end program 
      end


C*****SUBROUTINES

C*****random seed generator
        SUBROUTINE init_random_seed()
          INTEGER :: i, n, clock
          INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
          CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
          END SUBROUTINE
        

