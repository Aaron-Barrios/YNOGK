      include  '../ynogk.f90'
      PROGRAM  MAIN 
      use constants
      IMPLICIT NONE 

      double precision :: inte,r11,r22,eps,a1,a2,a3,a4,alpha,a,B,lz,q,a_spin,cobs,g3,z,b1,c1,d1,&
                  sn,u,emmc,cn,dn,muobs,sinobs,p,deltax,deltay,rin,rout,y,x,f1,g1,&
                  h1,f2,g2,h2,a5,b5,integ,robs,lambda,ac,lambda1,lambda2,lambdac,scal,&
                  tinf,tp,hp,b4,p4,p1J1p,mudisk,rdisk_out,theta_disk 
      complex*16 :: r1,r2,r3,r4
      integer :: del,i,j,NN,reals,DP,p5
      parameter (DP=kind(1.D0))
      character (len=80) :: outfilename       
      read(unit=5,fmt=*)a_spin,cobs,robs,scal,theta_disk,rdisk_out 

      If(cobs.ne.90.D0)then
          cobs=cobs*dtors
          muobs=cos(cobs)
          sinobs=sin(cobs)
      else
          muobs=0.D0
          sinobs=1.D0
      endif
      If(theta_disk.ne.zero)then
          mudisk=cos((90.D0-theta_disk)*dtors)
      Else
          mudisk=zero
      ENdif

      call reashiftofthindisk(sinobs,muobs,mudisk,rdisk_out,a_spin,robs,scal)
      END
!*****************************************************************************************************
      subroutine reashiftofthindisk(sinobs,muobs,mudisk,rdisk_out,a_spin,robs,scal)
!*****************************************************************************************************
      use blcoordinate
      implicit none
      double precision :: a,B,sinobs,muobs,a_spin,p,deltax,deltay,ra,robs,lambda,q,scal,&
                  Delta,sigma,bigA,somiga,expmu,pem,re,somiga1,rms1,&
                  bomiga,ut,g,mudisk,rdisk_out,bomigak,expnu,exppsi,expmu1,expmu2,&
                  alpha,beta,velocity(3),alphac,Betac,moment(4),f1234(4),&
                  somiga_obs,expnu_obs,exppsi_obs,expmu1_obs,expmu2_obs,&
                  somiga_em,expnu_em,exppsi_em,expmu1_em,expmu2_em,ut_em       
      integer :: i,j,m

      call metricg(robs,sinobs,muobs,a_spin,somiga_obs,expnu_obs,exppsi_obs,expmu1_obs,expmu2_obs)  
      rms1 = rms(a_spin)
      m=800
! (111) in Yang & Wang (2012).
      velocity(1)= robs/(robs**(three/two)+a_spin)*zero
      velocity(2)= robs/(robs**(three/two)+a_spin)*cos(45.D0*dtors)*zero
      velocity(3)= robs/(robs**(three/two)+a_spin)*sin(45.D0*dtors)*zero 
! (104) and (105) in Yang & Wang (2012).     
      call center_of_image(robs,sinobs,muobs,a_spin,scal,velocity,alphac,betac)
      write(*,*)alphac,betac
      deltax=60.D0/m
      deltay=28.D0/m
      open(unit=15,file='tdiskg.txt',status="replace")                    
      DO i=0,m
          Beta=betac-i*deltay+13.D0 
          Do j=0,m
              alpha=alphac-j*deltax+30.D0  
              call lambdaq(alpha,Beta,robs,sinobs,muobs,a_spin,scal,velocity,f1234,lambda,q) 
! call pemdisk just gives the direct image.            
              !pem = Pemdisk(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mudisk,rdisk_out,rms1)  
! while call pemdisk_all will gives the direct and high-order images.            
              pem = Pemdisk_all(f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,mudisk,rdisk_out,rms1) 
              !write(*,*)'pem == ',pem  
              If(pem.ne.-one.and.pem.ne.-two)then
                  re = radius(pem,f1234(1),lambda,q,a_spin,robs,scal)    
!Keplerian velocity of the particle.    
                  bomiga=one/(a_spin+re**(three/two))
                  call metricg(re,one,zero,a_spin,somiga_em,expnu_em,exppsi_em,expmu1_em,expmu2_em) 
! time component of four momentum of particle.
                  ut_em=one/expnu_em/sqrt(one-(exppsi_em/expnu_em*(bomiga-somiga_em))**two) 
! redshift formula of (113) of Yang and Wang (2012).
                  g=(one-somiga_obs*lambda)/expnu_obs/f1234(4)/(one-bomiga*lambda)/ut_em
                  write(unit=15,fmt=*)g 
              else  
                  If(pem .eq. -one)then
                      write(unit=15,fmt=*)zero!sqrt(-re)
                  else
                      write(unit=15,fmt=*)0.1D0!sqrt(-re)
                  endif 
              endif           
          Enddo
      Enddo
      close(unit=15)      
      end subroutine reashiftofthindisk
!*****************************************************************************************************
