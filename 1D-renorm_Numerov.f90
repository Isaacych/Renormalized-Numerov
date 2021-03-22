!-------------------------------------------------------------------------------------------------------------
!Purpose of the code: 
!1. Implement the renormlized Numerov method from B. R. Johnson, J. Chem. Phys, 67,9,1, 1977
!The coupling between closed and open channel will be modeled by a Lorentzian function : gamma*C^2/[(R-R_c)^2+C^2]
!The potential used is a 12,6 Lennard-Jones potential
!V22== 4.0*eps2*((sig/R)**12-(sig/R)**6)
!eps2 = 40/au_cm
!sig = 6*a0
!-------------------------------------------------------------------------------------------------------------
Program OneD_SE
Implicit none
Integer iR,Rsteps,iE,Esteps,Num_ch,L,scat_bound,midstep,step2,jR,scount
Real*8 mu,au_cm,au_eV,pi,amu
Real*8 nrg,nrg_ini,nrg_fin,dE,dR,R_ini,R_fin,eps,sig,pot,VL,U1,W1
Real*8 K_mat,f,g,fp,gp,kR,kR2,xs,k,f2,g2,norm,D_E,Vmin,D_prev
Complex*16 S_mat,ci
Real*8, allocatable:: Q1(:),Qb(:),R(:),Psi(:),Psib(:)
Character*20 filename
!------------Set up the parameters---------------------
amu=1822.888486209d0
mu=1.0078d0*amu
pi=dacos(-1.d0)
au_cm=219474.63136320d0
au_eV=27.211386245988d0
ci=(0.d0,1.d0)

R_ini=3.d0;R_fin=500.d0;

dR=1.d-2;Rsteps=ceiling((R_fin-R_ini)/dR)

!-----------------------------------------------------------
Num_ch=1
scat_bound=1	!scat_bound=0 or 1 : bound state/scattering calculation
!--------------single channel scattering case--------------
if ((Num_ch==1).and.(scat_bound==1)) then
	allocate(R(Rsteps),Q1(Rsteps),Psi(Rsteps))
	Open(51,file='Xsec.dat')
	open(61,file='wavefn.dat')
	nrg_ini=1.d-1/au_cm;nrg_fin=1.d1/au_cm;
	Esteps=1000;dE=(nrg_fin-nrg_ini)/dble(Esteps)
	eps=40.d0/au_cm;sig=6.d0 !parameter for Lennard-Jones potential
	L=3
	Do iE=1,Esteps
		Q1=0.d0;Psi=0.d0;norm=0.d0;
		nrg=nrg_ini+dble(iE-1)*dE
		!---------Propagate using Renormalized Numerov method------------
		Do iR=2,Rsteps
			R(iR)=R_ini+dble(iR-1)*dR
			call Pot6_12LJ(eps,sig,R(iR),pot)
			call Centrifugal(L,mu,R(iR),VL)
			W1=pot+VL-nrg	
			W1=2.d0*mu*W1
			U1=(2.d0+5.d0/6.d0*dR**2*W1)/(1.d0-dR**2/12.d0*W1)	
			if (iR==2) then
				Q1(iR)=U1
				Psi(iR)=1.d0
			else
				Q1(iR)=U1-1.d0/Q1(iR-1)
				Psi(iR)=Psi(iR-1)*Q1(iR)
			end if
		End do
		
		k=dsqrt(2.d0*mu*nrg)
		kR=k*R(Rsteps)
		kR2=k*R(Rsteps-1)
		
		call sphbes(L,kR,f,g,fp,gp)
		f=kR*dsqrt(2.d0*mu/(pi*k))*f
		g=kR*dsqrt(2.d0*mu/(pi*k))*g
		
		call sphbes(L,kR2,f2,g2,fp,gp)
		f2=kR2*dsqrt(2.d0*mu/(pi*k))*f2
		g2=kR2*dsqrt(2.d0*mu/(pi*k))*g2
		
		K_mat=(f-Q1(Rsteps)*f2)/(g-Q1(Rsteps)*g2)
		S_mat=(1.d0+ci*K_mat)/(1.d0-ci*K_mat)
		xs=pi/k**2*(2.d0*dble(L)+1.d0)*zabs(S_mat-1.d0)**2	
		write(51,*) nrg*au_cm,xs
		
		!Normalize the wave function, and print it.
		Do iR=1,Rsteps
			norm=dR*dabs(Psi(iR))**2+norm
		End do
		norm=norm-0.5d0*dR*dabs(Psi(Rsteps))**2	!Composite Trapezoidal rule is accurate enough.
		Psi=Psi/dsqrt(norm)
		
		if(iE==500) then
			Do iR=2,Rsteps
			write(61,*) R(iR),Psi(iR)
			End do
		end if	
	End do
	close(51);close(61)	
	deallocate(R,Q1,Psi)
end if
!--------------End---single channel scattering case---------------------
!--------------single channel bound state case--------------------------
if ((Num_ch==1).and.(scat_bound==0)) then
!---Setting up the problem-------------------------------------------	
	nrg_ini=-100.d0;nrg_fin=1.5d0;
	Esteps=1000;dE=(nrg_fin-nrg_ini)/dble(Esteps)
	eps=1.d2;sig=1.d0 !parameter for Lennard-Jones potential
	scount=0
	mu=1.d0
	R_ini=0.6d0;R_fin=3.d0;
	Rsteps=5000
	dR=(R_fin-R_ini)/(Rsteps-1)
	allocate(R(Rsteps))
	open(61,file='D_E.dat')
	open(511,file='bound_pot.dat')
	Vmin=0.d0
	midstep=1
	Do iR=2,Rsteps
		R(iR)=R_ini+dble(iR-1)*dR
		call Pot6_12LJ(eps,sig,R(iR),pot)
		write(511,*) R(iR),pot
		if(pot<Vmin) then
			Vmin=pot
			midstep=iR
		end if	
	End do
	step2=Rsteps-midstep-1
	allocate(Q1(midstep),Psi(Rsteps),Qb(step2),Psib(step2+1))
	close(511)
	print *,'Rstep,midstep,Vmin=',Rsteps,midstep,Vmin
!--------------------------------------------------------------------------------	
	Do iE=1,Esteps
		Q1=0.d0;Qb=0.d0;Psi=0.d0;Psib=0.d0;norm=0.d0;D_E=0.d0
		nrg=nrg_ini+dble(iE-1)*dE
		!---------Propagate using Renormalized Numerov method------------
		Do iR=2,midstep
			R(iR)=R_ini+dble(iR-1)*dR
			call Pot6_12LJ(eps,sig,R(iR),pot)
			W1=pot-nrg	
			W1=2.d0*mu*W1
			U1=(2.d0+5.d0/6.d0*dR**2*W1)/(1.d0-dR**2/12.d0*W1)	
			if (iR==2) then
				Q1(iR)=U1
				Psi(iR)=1.d0
			else
				Q1(iR)=U1-1.d0/Q1(iR-1)
				Psi(iR)=Psi(iR-1)*Q1(iR)
			end if
		End do
		!---------backward propagate using Renormalized Numerov method--------
		Do jR=1,step2
			iR=Rsteps-jR
			R(iR)=R_ini+dble(iR-1)*dR
			call Pot6_12LJ(eps,sig,R(iR),pot)
			W1=pot-nrg	
			W1=2.d0*mu*W1
			U1=(2.d0+5.d0/6.d0*dR**2*W1)/(1.d0-dR**2/12.d0*W1)
			if (iR==Rsteps-1) then
				Qb(jR)=U1
				Psib(jR)=1.d0
				Psib(jR+1)=Psib(jR)*Qb(jR)
			else
				Qb(jR)=U1-1.d0/Qb(jR-1)
				Psib(jR+1)=Psib(jR)*Qb(jR)
			end if	
		End do
	!--------------Match the forward and backward propagated results--------------		
		D_E=1.d0/Qb(step2)-Q1(midstep)
		write(61,*)nrg,D_E
		!Detection of eigenvalues
		if((D_E*D_prev<0.d0).and.(iE>1)) then
			if(dabs(D_E)<dE) then
				print *,'Eigenvalue found!'
				print *,nrg,D_E
				!Eigenvalue found
				scount=scount+1
				Write(filename,'("Psi_",I1,".dat")')scount
				open(66,file=filename)
				!Fix the ratio between forward and backward prop. wavefn.
				norm=Psi(midstep)/Psib(step2+1)
				Psib=Psib*norm
				!Copy the backward wf to full wf
				Do iR=midstep+1,Rsteps-1
					jR=Rsteps-iR
					Psi(iR)=Psib(jR)
				End do
				
				!Normalize the wave function, and print it.
				Do iR=1,Rsteps
					norm=dR*dabs(Psi(iR))**2+norm
				End do
				Psi=Psi/dsqrt(norm)
				Do iR=2,Rsteps
					write(66,*) R(iR),Psi(iR)
				End do
				close(66)
			end if
		end if
		D_prev=D_E
	End do	!End loop for energy
	close(61)	
	deallocate(R,Q1,Qb,Psi,Psib)
end if
!--------------End---single channel bound state case--------------------
End program

!*******************************************************************
	SUBROUTINE Pot6_12LJ(eps,sig,x,pot)
	Implicit none
	Real*8 eps,sig,x,pot
	pot=4.d0*eps*((sig/x)**12-(sig/x)**6)	
	END SUBROUTINE Pot6_12LJ
!*******************************************************************
	SUBROUTINE Centrifugal(L,mu,x,pot)
	Implicit none
	Integer L
	Real*8 mu,x,pot
	pot=dble(L*(L+1))/(2.d0*mu*x**2)	
	END SUBROUTINE Centrifugal
!*******************************************************************

	  SUBROUTINE sphbes(n,x,sj,sy,sjp,syp)

      implicit none
      INTEGER n
      REAL*8 sj,sjp,sy,syp,x,factor,order,rj,rjp,ry,ryp,RTPIO2
!U    USES bessjy

      RTPIO2=sqrt(asin(1.d0))
      if(n < 0.d0.or.x<=0.d0)stop 'bad arguments in sphbes'
      order=n+0.5d0
      call bessjy(x,order,rj,ry,rjp,ryp)
      factor=RTPIO2/sqrt(x)
      sj=factor*rj
      sy=factor*ry
      sjp=factor*rjp-sj/(2.d0*x)
      syp=factor*ryp-sy/(2.d0*x)
      return
	  END SUBROUTINE sphbes

! **************************************
	  SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)

      implicit none
      INTEGER::MAXIT=10000,i,isign,l,nl
      REAL*8::rj,rjp,ry,ryp,x,xnu,PI=3.141592653589793d0
      REAL*8::EPS=1.d-16,FPMIN=1.d-30,XMIN=2.d0,q,r,sum1
!U    USES beschb
      REAL*8 a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e, f,fact,fact2,fact3,ff,gam,gam1,gam2,gammi,gampl,&
       rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,   temp,w,x2,xi,xi2,xmu,xmu2,p,pimu,pimu2
      if(x<=0.d0.or.xnu < 0.d0) stop 'bad arguments in bessjy'
      if(x < XMIN)then
        nl=int(xnu+.5d0)
      else
        nl=max(0,int(xnu-x+1.5d0))
      endif
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      w=xi2/PI
      isign=1
      dr=xnu*xi
      if(dr < FPMIN)dr=FPMIN
      b=xi2*xnu
      d=0.d0
      c=dr
      do i=1,MAXIT
        b=b+xi2
        d=b-d
        if(abs(d) < FPMIN)d=FPMIN
        c=b-1.d0/c
        if(abs(c) < FPMIN)c=FPMIN
        d=1.d0/d
        del=c*d
        dr=del*dr
        if(d < 0.d0)isign=-isign
        if(abs(del-1.d0) < EPS)goto 1
      end do
      print*,x
      stop 'x too large in bessjy; try asymptotic expansion'
1     continue
      rjl=isign*FPMIN
      rjpl=dr*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do l=nl,1,-1
        rjtemp=fact*rjl+rjpl
        fact=fact-xi
        rjpl=fact*rjtemp-rjl
        rjl=rjtemp
      end do
      if(rjl==0.d0)rjl=EPS
      f=rjpl/rjl
      if(x < XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu) < EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e) < EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=2.d0/PI*fact*(gam1*cosh(e)+gam2*fact2*d)
        e=exp(e)
        p=e/(gampl*PI)
        q=1.d0/(e*PI*gammi)
        pimu2=0.5d0*pimu
        if(abs(pimu2) < EPS)then
          fact3=1.d0
        else
          fact3=sin(pimu2)/pimu2
        endif
        r=PI*pimu2*fact3*fact3
        c=1.d0
        d=-x2*x2
        sum=ff+r*q
        sum1=p
        do i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          if(abs(del) < (1.d0+abs(sum))*EPS)goto 2
        end do
        stop 'bessy series failed to converge'
2       continue
        rymu=-sum
        ry1=-sum1*xi2
        rymup=xmu*xi*rymu-ry1
        rjmu=w/(rymup-f*rymu)
      else
        a=.25d0-xmu2
        p=-.5d0*xi
        q=1.d0
        br=2.d0*x
        bi=2.d0
        fact=a*xi/(p*p+q*q)
        cr=br+q*fact
        ci=bi+p*fact
        den=br*br+bi*bi
        dr=br/den
        di=-bi/den
        dlr=cr*dr-ci*di
        dli=cr*di+ci*dr
        temp=p*dlr-q*dli
        q=p*dli+q*dlr
        p=temp
        do i=2,MAXIT
          a=a+2*(i-1)
          bi=bi+2.d0
          dr=a*dr+br
          di=a*di+bi
          if(abs(dr)+abs(di) < FPMIN)dr=FPMIN
          fact=a/(cr*cr+ci*ci)
          cr=br+cr*fact
          ci=bi-ci*fact
          if(abs(cr)+abs(ci) < FPMIN)cr=FPMIN
          den=dr*dr+di*di
          dr=dr/den
          di=-di/den
          dlr=cr*dr-ci*di
          dli=cr*di+ci*dr
          temp=p*dlr-q*dli
          q=p*dli+q*dlr
          p=temp
          if(abs(dlr-1.d0)+abs(dli) < EPS)goto 3
        end do
        stop 'cf2 failed in bessjy'
3       continue
        gam=(p-f)/q
        rjmu=sqrt(w/((p-f)*gam+q))
        rjmu=sign(rjmu,rjl)
        rymu=rjmu*gam
        rymup=rymu*(p+q/gam)
        ry1=xmu*xi*rymu-rymup
      endif
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do i=1,nl
        rytemp=(xmu+i)*xi2*ry1-rymu
        rymu=ry1
        ry1=rytemp
      end do
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      return
	  END SUBROUTINE bessjy

!   ****************************************************************
	  SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER::NUSE1=7,NUSE2=8
      REAL*8 gam1,gam2,gammi,gampl,x,xx,c1(7),c2(8),chebev,modin,odin
!U    USES chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,-4.971736704d-6,-3.3126120d-8,2.42310d-10,&
        -1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      modin=-1.d0
      odin=1.d0
      gam1=chebev(modin,odin,c1,NUSE1,xx)
      gam2=chebev(modin,odin,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
	  END SUBROUTINE beschb

 !   ********************************************
	  Real*8 FUNCTION chebev(a,b,c,m,x)
      INTEGER m,j
      REAL*8 a,b,x,c(m),d,dd,sv,y,y2
      if ((x-a)*(x-b) > 0.d0) stop 'x not in range in chebev'
      d=0.d0
      dd=0.d0
      y=(2.d0*x-a-b)/(b-a)
      y2=2.d0*y
      do j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
      end do
      chebev=y*d-dd+0.5d0*c(1)
      return
	  END FUNCTION chebev
! ***************************************************************
 
