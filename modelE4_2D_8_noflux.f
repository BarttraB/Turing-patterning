C    program modelE4_1.f    by Bartlomiej Borek , Sept. 2013

C    **this is basic simulation of model E4 in 2D

      implicit none
      integer, parameter :: nx=261,ny=261
      integer, parameter :: ntime_end=30000000
      double precision, parameter :: pi=3.14159265358979323846
      double precision L(0:nx+1,0:nx+1), H(0:nx+1,0:nx+1)
      double precision A(0:nx+1,0:nx+1), P(0:nx+1,0:nx+1)
      double precision Lt(0:nx+1,0:nx+1), Ht(0:nx+1,0:nx+1)
      double precision At(0:nx+1,0:nx+1), Pt(0:nx+1,0:nx+1)
      double precision lapL, lapH, dL, dH, dA, dP
      double precision rspace, rtime, r1, pert
      double precision D,dx,dt
      double precision a1,a2,g1,g2,a3,g3
      double precision d1,d2
      double precision c1,c2,c3,c4,c5
      double precision Pm,k1,k2
      double precision f1,f2,f3
      integer i1,i2,i,ii,j,jj,j2,iq,ntime,fntime

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
     &    2.5d0,1.1d0,0.675d0/
      data Pm, k1,k2/
     &    1.0d0,1.0d2,2.4815d2/

C--------> open files to write to
      open(unit=90,file='L_E4_2D_8_noflux.txt', status='unknown')
c      open(unit=91,file='H_E4_2D_2.txt', status='unknown')
c      open(unit=92,file='A_E4_2D_2.txt', status='unknown')
c      open(unit=93,file='P_E4_2D_2.txt', status='unknown')

C--------> set initial conditions
       do jj=0,nx+1
         do ii=0,nx+1
	      rspace=(ii-1)*(jj-1)*nx*nx
c           L(ii)=0.005445572791697402d0
c            H(ii)=0.005194499371405287d0
C % the stable fixed point of the homogenous system is numerically found for this set of parameters (see modelE1_simplest1_4.nb code). For the dx and dt chosen here the convergence is to L=0.005213230849236,H=0.005080131526230
C % for 3V model A=0.00521323084940534

c------> random noise w.r.t. space IC
            call random_number(r1)
c		  2*Pi*rspace/(dx*(nx-1))
		  write(99,*) pert,ii,jj
		  pert=0.001*r1
c 	       instead of r1 consider SIN(rspace*0.1/2*pi) 	
c		  write(99,*) pert, rspace, pi
            L(ii,jj)=0.0099537412903138415697630222462534+pert
            H(ii,jj)=0.0073825689719981407350966878264083
            A(ii,jj)=0.019670567928752395134518190011574
            P(ii,jj)=0.0039951539638696398637446589771259
c     +0.1+0.0002*(r1-0.5)
c           L(ii)=0.0099537412903138415697630222462534
         enddo
	  enddo

c         L(100)=0.05
c         r(135:165)=0.5


C************ Start main iteration loop *************
         do ntime=1,ntime_end
            fntime=ntime-1
            rtime=(fntime-1.0)*dt

C----->write spacetime files

            if(mod(fntime,100000)==0)then
              do i=0,nx+1
c               write(99,*) fntime
                write(90,*) (L(i,j),j=0,nx+1)
c               write(91,*) H(1:nx,1:nx)
c               write(92,*) A(1:nx,1:nx)
c               write(93,*) P(1:nx,1:nx)
		    enddo
            endif
         
C----------> integrate the equations (forward Euler and first order Laplacian)
          do j=1,nx
            do i=1,nx
       
              lapL=(L(i+1,j)+L(i-1,j)+L(i,j-1)+L(i,j+1)-4*L(i,j))
c               lapL=0.0             
	         lapH=(H(i+1,j)+H(i-1,j)+H(i,j-1)+H(i,j+1)-4*H(i,j))
c               lapH=0.0
	       
	         f1=a1*(d1+P(i,j))/(P(i,j)+c1)
	         f2=g1*A(i,j)*(L(i,j)/(L(i,j)+c3))
		    f3=k1*(Pm-P(i,j))*L(i,j)-k2*P(i,j)

	         dL=f1-f2-f3
	         dH=a2*(d2+P(i,j))/(P(i,j)+c4)-g2*H(i,j)
	         dA=(a3*H(i,j))/(H(i,j)+c5)-g3*A(i,j)
	         dP=f3

              Lt(i,j)=L(i,j)+dt*dL+dt*lapL/(dx**2)
              Ht(i,j)=H(i,j)+dt*dH+dt*D*lapH/(dx**2)
              At(i,j)=A(i,j)+dt*dA
              Pt(i,j)=P(i,j)+dt*dP

            enddo 
          enddo 

C---------> No-flux boundary conditions for L,H,A,P
	   do i=0,nx+1
	      Lt(nx+1,i)=Lt(nx,i)
	      Lt(0,i)=Lt(1,i)
	      Ht(nx+1,i)=Ht(nx,i)
	      Ht(0,i)=Ht(1,i)
	      At(nx+1,i)=At(nx,i)
	      At(0,i)=At(1,i)
           Pt(nx+1,i)=Pt(nx,i)
           Pt(0,i)=Pt(1,i)   
	   enddo  
	   
	   do j=0,ny+1
	      Lt(j,ny+1)=Lt(j,ny)
	      Lt(j,0)=Lt(j,1)
	      Ht(j,ny+1)=Ht(j,nx)
	      Ht(j,0)=Ht(j,1)
	      At(j,ny+1)=At(j,nx)
	      At(j,0)=At(j,1)
	      Pt(j,ny+1)=Pt(j,nx)
	      Pt(j,0)=Pt(j,1)
	   enddo

C---------> update state variables
          do j2=0,nx+1
            do i2=0,nx+1
               L(i2,j2)=Lt(i2,j2)
               H(i2,j2)=Ht(i2,j2)
               A(i2,j2)=At(i2,j2)
               P(i2,j2)=Pt(i2,j2)
            enddo
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
        

