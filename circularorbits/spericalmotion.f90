      include  '../ynogk.f90'  
      PROGRAM  MAIN
      use constants      
      IMPLICIT NONE 

      Double precision  inte,r11,r22,eps,a1,a2,a3,a4,alpha,a,B,lz,q,a_spin,cobs,g3,z,b1,c1,d1,&
                  sn,u,emmc,cn,dn,muobs,sinobs,p,deltax,deltay,rin,rout,y,x,f1,g1,&
                  h1,f2,g2,h2,a5,b5,integ,robs,lambda,ac,lambda1,lambda2,lambdac,scal,&
                  tinf,tp,hp,b4,p4,p1J1p,theta_min 
      complex r1,r2,r3,r4
      integer del,i,j,NN,reals,DP,p5
      parameter (DP=kind(1.D0) )
      character (len=80) :: outfilename       

      read(unit=5,fmt=*)a_spin,cobs,robs 

      If(cobs.ne.90.D0)then
          cobs=cobs*PI/180.D0
          muobs=cos(cobs)
          sinobs=sin(cobs)
      else
          muobs=0.D0
          sinobs=1.D0
      endif   
      
      IF(a_spin*a_spin.ge.one)Then
          write(*,*)'e^2+a^2>=1, there can be no horizon. wrong and stoped.'
          stop
      ENDIF
      call reashift2(sinobs,muobs,a_spin,robs,scal)
      END
!***************************************************************************************************** 
      subroutine reashift2(sinobs,muobs,a_spin,robs,scal)
!*****************************************************************************************************
      use blcoordinate
      use sphericalmotion
      implicit none
      Double precision  a,B,sinobs,muobs,a_spin,delta_theta,deltap,robs,lambda,q,scal,&
             Delta,sigma,bigA,somiga,expmu,pem,re,vptl(3),velo,mve,signs,&
             bomiga,ut,g,p,rhorizon,r_tp,velocity(3),kvec(4),moment(4),&
             rin,rout,muup,mudown,phy1,phy2,dphy,alpha,Beta,alphac,rmss,&
             betac,theta,phi_i,pr,pt,pp,PE,sint,sinp,cost,cosp,ra,mua,phya,&
             x,y,z,delta_phi,expnu,exppsi,expmu1,expmu2,r_in,timea,sigmaa,&
             theta_min,ep,cobs,ptotal,xcc,ycc,zcc   
      integer  i,j,m,caserange,n,k,t1,t2 
      logical :: cir_orbt

      cir_orbt = .true.
      call metricg(robs,sinobs,muobs,a_spin,somiga,expnu,exppsi,expmu1,expmu2)
      velocity(1)=0.0D0!expmu1/expnu*robs/(robs**(2.05D0/two)+a_spin)!*zero
      velocity(2)=0.0D0!expmu2/expnu/(robs**(three/two)+a_spin)*zero
      velocity(3)=exppsi/expnu*(-somiga)  !exppsi/expnu*(one/(robs**(three/two)+a_spin)-somiga) 
      velo=0.35D0 
      m=20
      n=20 
      rhorizon=one+sqrt(one-a_spin**two) 
      delta_theta=180.D0/n
      delta_phi=360.D0/m
      deltap=10.D0/200.D0
      t1=0 
      rmss = rms(a_spin) 
      signs=one

      open(unit=15,file='rayx1.txt',status="replace")
      open(unit=16,file='rayy1.txt',status="replace")        
      open(unit=17,file='rayz1.txt',status="replace")                            
      DO i=m,m
            !write(*,*)i      
            phi_i=i*delta_phi*dtor
            !write(*,*)'phi=',phi
            Do j=n,n
                  theta=90.D0*dtor!j*delta_theta*dtor
                  !write(*,*)'theta=',theta
                  vptl(2)= velo*cos(theta)
                  vptl(1)= velo*sin(theta)*sin(phi_i)
                  vptl(3)= velo*sin(theta)*cos(phi_i)   
                  !call leq(r_circle,theta_min,signs,a_spin,lambda,q) 
                  !call leq(robs,theta_min,signs,a_spin,lambda,q) 
                  call lambdaq_sphericalm(robs,a_spin,lambda,q,theta_min)
                  Do k=0,20000
                      p=k*0.006!deltap  
                      call SPHERICALMOTION_BL(p,one,one,lambda,q,sinobs,muobs,a_spin,&
                               robs,theta_min,phya,timea,sigmaa,mua,t1,t2)   
                      ra = robs    
                      If(ra .ge. rhorizon*0.99)then 
                          xcc=sqrt(ra**two+a_spin**two)*sqrt(one-mua**two)*cos(phya) 
                          ycc=sqrt(ra**two+a_spin**two)*sqrt(one-mua**two)*sin(phya)
                          zcc=sqrt(ra**two+a_spin**two)*mua 
                          write(unit=15,fmt=*)xcc 
                          write(unit=16,fmt=*)ycc
                          write(unit=17,fmt=*)zcc  
                      endif
                  Enddo
            Enddo
      Enddo
      close(unit=15)      
      close(unit=16)
      close(unit=17) 
      return
      end subroutine reashift2










