      include  '../ynogk.f90'
      PROGRAM  MAIN
      use constants
      IMPLICIT NONE 

      Double precision inte,r11,r22,eps,a1,a2,a3,a4,alpha,a,B,lz,q,a_spin,cobs,cobs2,g3,z,b1,c1,d1,&
                  sn,u,emmc,cn,dn,muobs,sinobs,p,deltax,deltay,rin,rout,y,x,f1,g1,&
                  h1,f2,g2,h2,a5,b5,integ,robs,lambda,ac,lambda1,lambda2,lambdac,scal,&
                  tinf,tp,hp,b4,p4,p1J1p,gama0 
      Double precision,external :: Fp      
      complex*16 r1,r2,r3,r4
      integer del,i,j,NN,reals,DP,p5
      PARAMETER (DP=kind(1.D0))
      character (len=80) :: outfilename       
      read(unit=5,fmt=*)a_spin,cobs,robs,scal,gama0

      If(cobs.ne.90.D0)then
          cobs2=cobs*dtor
          muobs=dcos(cobs2)
          sinobs=sin(cobs2)
      else
          muobs=0.D0
          sinobs=1.D0
      endif
      If(cobs.eq.180.D0)then
          muobs=-one
          sinobs=zero
      endif

      call redshttori(sinobs,muobs,cobs,a_spin,robs,scal,Fp,gama0)
      END
!*****************************************************************************************************
      Function Fp(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,paras)
!*****************************************************************************************************
      use blcoordinate
      implicit none
      Double precision Fp,a,B,sinobs,muobs,mudisk,a_spin,robs,scal,p,re,mua,rdisk,z,rout,rin,&
                   beta,gamma0,phya,n3,lambda,q,paras(10),f1234(4),mup,rp,&
                  mu,sigma,timec,timea,sigmaa   
      character varble
      integer t1,t2,tr1,tr2,tm1,tm2
       
      n3=paras(1)!0.95D0      
      rin=paras(2)!4.23D0
      rout=paras(3)!50.D0       
      !re=radius(p,f1234(1),lambda,q,a_spin,robs,scal)
      !mua=mucos(p,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
      !phya=phi(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal)   
      call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                       re,mua,phya,timea,sigmaa)  
      !call GEOKERR(p,rp,mup,varble,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,t1,t2,re,mua,timec,phya,sigma)
! (126) of Yang & Wang (2012). 
      beta=n3*sin(PI/two*(re-rin)/(rout-rin))
! (127) of Yang & Wang (2012).
      gamma0=paras(4)*dtors+paras(5)*PI*exp(paras(6)*(rin-re)/(rout-rin))
! (128) of Yang & Wang (2012).
      Fp=tan(beta)*cos(phya-gamma0)+mua/sqrt(one-mua**two)  
      !write(unit=6,fmt=*)re,mua,phya,Fp
      return
  End Function Fp
!*****************************************************************************************************
 subroutine redshttori(sinobs,muobs,cobs,a_spin,robs,scal,Fp,gama0)
!*****************************************************************************************************
      use blcoordinate
      use pemfinding
      implicit none
      Double precision a,B,sinobs,muobs,a_spin,p,deltax,deltay,ra,robs,lambda,q,scal,&
                  Delta,sigma,bigA,somiga,expmu,pem,re,gama0,timea,sigmaa,&
                  bomiga,ut,g,rdisk_out,bomigak,beta,gamma0,theta_dot,phi_dot,&
                  mua,phya,V_phi,V_theta,V,deltap,dmu,sq_Theta,ut_em,V_the_bar,V_phi_bar,&
                  rin,rout,muup,mudown,phy1,phy2,phy0,n3,n1,n2,sinp,V_r,&
                  somiga_obs,expnu_obs,exppsi_obs,expmu1_obs,expmu2_obs,&
                  somiga_em,expnu_em,exppsi_em,expmu1_em,expmu2_em,sin_phy0,cos_phy0,&
                  paras(10),f1234(4),alphac,Betac,alpha1,Beta1,velocity(3),cobs 
      Double precision,external :: Fp!,unifyp
      integer :: i,j,m,caserange,ten=10
      PARAMETER(deltap=1.D-5)
      logical :: bisection
      character(len=100) :: fname, ifname(10), fileimage,IMGEXT = 'txt'


      call metricg(robs,sinobs,muobs,a_spin,somiga_obs,&
                   expnu_obs,exppsi_obs,expmu1_obs,expmu2_obs)

      velocity(1)=-0.D0!expmu1_obs/expnu_obs*robs/(robs**(three/two)+a_spin)!*zero
      velocity(2)=0.D0!expmu2_obs/expnu_obs/(robs**(three/two)+a_spin)*cos(45.D0*dtor)!*zero
      velocity(3)=0.D0!exppsi_obs/expnu_obs*(one/(robs**(three/two)+a_spin)*sin(45.D0*dtor)-somiga_obs)!*zero      
      call center_of_image(robs,sinobs,muobs,a_spin,scal,velocity,alphac,Betac)
      write(*,*)alphac,Betac

      m=400
      deltax=110.D0/m
      deltay=110.D0/m
! n1, n2, n3 are parameters in (126) of Yang & Wang (2012).
      n1=4.D0
      n2=4.D0
      n3=0.95D0 
      rin=rms(a_spin)
      rout=50.D0 
! caserange=4, means muup, mudown and phy1, phy2 are not provided.
      caserange=4   
! parameters to describe the curved surface of warped disk and sended to
! function Fp as dummy variable.   
      paras(1)=n3
      paras(2)=rin
      paras(3)=rout
      paras(4)=gama0
      paras(5)=n1
      paras(6)=n2
!bisection = .false., we use Newton-Raphson method to search the root of equation
! f(p)=0.        
      bisection = .false.
! The name of data file. 
      !write(unit = fname,fmt = '("warpdisk1_",i4.4,".",A3)')floor(abs(paras(4))),IMGEXT
! open data file.
      open(unit=15,file="warpdisk.txt",status="replace")        
!begin to compute the images, from corner (60.D0, 60.D0) of a retangular on 
! the secreen of observer.                  
      DO i=0,m
          Beta1=55.D0-i*deltax+Betac
          write(unit=6,fmt=*)i
          Do j=0,m
              alpha1=55.D0-j*deltay+alphac
! To get the motion constants from subroutine lambdaq.
              call lambdaq(alpha1,Beta1,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)
! Using function pemfind to search the root pem of equation f(p)=0.         
              pem=pemfind(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&                              
                        rin,rout,muup,mudown,phy1,phy2,caserange,Fp,paras,bisection)   
              !write(unit=6,fmt=*)pem
!If the photon do not go to infinity (pem=-1) and fall into the event horizon (pem=-2), use formula 
! (110) to compute the redshfit of photon emitted from the surface of the warped disk.
              If(pem.ne.-one.and.pem.ne.-two)then         
                  call YNOGK(pem,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                                 re,mua,phya,timea,sigmaa) 
                  !re=radius(pem,f1234(1),lambda,q,a_spin,robs,scal) 
                  !mua=mucos(pem,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)
                  sinp=sqrt(one-mua**two)
                  !phya=phi(pem,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal) 
                  !write(unit=6,fmt=*)re,mua,phya
                  beta=n3*sin(PI/two*(re-rin)/(rout-rin))
                  gamma0=paras(4)*dtor+paras(5)*PI*exp(paras(6)*(rin-re)/(rout-rin)) 

                  phy0=atan(tan(phya-gamma0)*cos(beta))
                  sin_phy0=sign(sin(phy0),sin(phya-gamma0))
                  cos_phy0=sign(cos(phy0),cos(phya-gamma0))
                  V=re/(re**(three/two)+a_spin)
                  V_theta=V*(-cos(phya-gamma0)*cos(beta)*sin_phy0*mua+sin(phya-gamma0)*cos_phy0*mua&
                        -sin(beta)*sin_phy0*sinp)
                  V_phi=V*(sin(phya-gamma0)*cos(beta)*sin_phy0+cos(phya-gamma0)*cos_phy0) 

                  theta_dot=V_theta/re
                  phi_dot=V_phi/re/sinp

                  call metricg(re,sinp,mua,a_spin,somiga_em,expnu_em,exppsi_em,expmu1_em,expmu2_em)
      
                  V_the_bar=theta_dot*expmu2_em/expnu_em
                  V_phi_bar=exppsi_em/expnu_em*(phi_dot-somiga_em)
                  ut_em=one/expnu_em/sqrt(one-V_the_bar**two-V_phi_bar**two)

                  dmu=mucos(pem+deltap,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal)-&
                              mucos(pem,f1234(3),f1234(2),lambda,q,sinobs,muobs,a_spin,scal) 
                  sq_Theta=sqrt(q+a_spin**two*mua**two-lambda**two*(mua/sinp)**two)    
! (110) of Yang & Wang (2012).  
                  g=one/expnu_obs*(one-lambda*somiga_obs)/f1234(4)/ut_em/(one+&
                               sign(sq_Theta,-dmu)*V_theta/re-lambda*V_phi/sinp/re) 
                  write(unit=15,fmt=*)g
              else
                  If(pem.eq.-one)then
! Photon goes to infinity.
                     write(unit=15,fmt=*)zero!sqrt(-re)
                  else
! Photon falls into the event horizon.
                     write(unit=15,fmt=*)0.1D0!sqrt(-re)
                  endif
              endif            
          Enddo
      Enddo
      close(unit=15)       
      end subroutine redshttori
!*****************************************************************************************************

