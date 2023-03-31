      include  '../ynogk.f90' 
      PROGRAM  MAIN
      use constants      
      IMPLICIT NONE 

      double precision :: inte,r11,r22,eps,a1,a2,a3,a4,alpha,a,B,lz,q,a_spin,cobs,g3,z,b1,c1,d1,&
             sn,u,emmc,cn,dn,muobs,sinobs,p,deltax,deltay,rin,rout,y,x,f1,g1,&
             h1,f2,g2,h2,a5,b5,integ,robs,lambda,ac,lambda1,lambda2,lambdac,scal,&
             tinf,tp,hp,b4,p4,p1J1p  
      integer :: del,i,j,NN,reals,DP,p5
      parameter (DP=kind(1.D0))
      character (len=80) :: outfilename 
      
      read(unit=5,fmt=*)a_spin,cobs,robs,scal 

      If(cobs.ne.90.D0)then
          cobs=cobs*dtor
          muobs=cos(cobs)
          sinobs=sin(cobs)
      else
          muobs=0.D0
          sinobs=1.D0
      endif 
      
      call reashift2(sinobs,muobs,a_spin,robs,scal)
      END
!*****************************************************************************************************
      subroutine reashift2(sinobs,muobs,a_spin,robs,scal)
!*****************************************************************************************************
      use blcoordinate
      implicit none
      double precision :: a,B,sinobs,muobs,a_spin,deltax,deltay,ra,robs,lambda,q,scal,&
             Delta,sigma,bigA,somiga,expmu,pem,re,&
             bomiga,ut,g,p,rhorizon,r_tp,velocity(3),f1234(4),moment(4),&
             rin,rout,muup,mudown,phy1,phy2,dphy,alpha,Beta,alphac,betac,x,y,&
             r_tp1,r_tp2,mua,phia,timea,sigmaa
      complex*16  bb(1:4)      
      integer :: i,j,m,caserange,t1,t2,reals,cases_of_tp 
      logical :: robs_eq_rtp,indrhorizon

      velocity(1)=zero
      velocity(2)=zero
      velocity(3)=zero
      call center_of_image(robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac)
      deltax=12.D0/800.D0
      deltay=12.D0/800.D0
      t1=0
      p=one
      rhorizon=one+sqrt(one-a_spin**two)
      x=43.D0
      y=1.D70
      open(unit=15,file='shadow.txt',status="replace")                    
      DO i=0,800
            write(*,*)i      
            Beta=6.D0-i*deltay+betac
            Do j=0,800
                  alpha=4.D0-j*deltax+alphac 
                  call lambdaq(alpha,beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q)  
                  call radiustp(f1234(1),a_spin,robs,lambda,q,r_tp1,&
                                r_tp2,reals,robs_eq_rtp,indrhorizon,cases_of_tp,bb)    
                  If(r_tp1.le.rhorizon)then
                      t1=0
                      t2=0
                      p = r2p(f1234(1),rhorizon,lambda,q,a_spin,robs,scal,t1,t2) 
                      call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                                                 ra,mua,phia,timea,sigmaa) 
                      write(unit=15,fmt=*)sigmaa 
                  else
                      t1=1
                      t2=0
                      p = r2p(f1234(1),robs,lambda,q,a_spin,robs,scal,t1,t2) 
                      call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                                                 ra,mua,phia,timea,sigmaa)  
                      write(unit=15,fmt=*)sigmaa*half 
                  endif      
            Enddo
      Enddo
      close(unit=15)    
      return
      end subroutine reashift2
