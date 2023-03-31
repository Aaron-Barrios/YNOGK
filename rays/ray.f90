     include  '../ynogk.f90' 
     PROGRAM  MAIN
     use blcoordinate
     IMPLICIT NONE 

     Double precision  inte,r11,r22,eps,a1,a2,a3,a4,alpha,a,B,lz,q,a_spin,cobs,g3,z,b1,c1,d1,&
            sn,u,emmc,cn,dn,muobs,sinobs,p,deltax,deltay,rin,rout,y,x,f1,g1,&
            h1,f2,g2,h2,a5,b5,integ,robs,lambda,ac,lambda1,lambda2,lambdac ,scal,&
            tinf,tp,hp,b4,p4,p1J1p         
      complex r1,r2,r3,r4
      integer del,i,j,NN,reals,DP,p5
      parameter (DP=kind(1.D0) )
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
      a=4.D0
      B=5.D0        
      robs = rms(a_spin)
        
      call reashift2(sinobs,muobs,a_spin,robs,scal)
      END 
!*****************************************************************************************************
      subroutine reashift2(sinobs,muobs,a_spin,robs,scal)
!*****************************************************************************************************
        use blcoordinate
        implicit none
        Double precision  a,B,sinobs,muobs,a_spin,delta_theta,deltap,robs,lambda,q,scal,&
                        Delta,sigma,bigA,somiga,expmu,pem,re,ptotal,&
                        bomiga,ut,g,p,rhorizon,r_tp,velocity(3),f1234(4),moment(4),&
                        rin,rout,muup,mudown,phy1,phy2,dphy,alpha,Beta,alphac,&
                        betac,theta,phi_i,pr,pt,pp,PE,sint,sinp,cost,cosp,ra,mua,phya,&
                        x,y,z,delta_phi,expnu,exppsi,expmu1,expmu2,r_in,timea ,sigmaa 
        integer  i,j,m,caserange,t1,n,k 

        call metricg(robs,sinobs,muobs,a_spin,somiga,expnu,exppsi,expmu1,expmu2)
        velocity(1)=0.0D0!expmu1/expnu*robs/(robs**(2.05D0/two)+a_spin)!*zero
        velocity(2)=0.0D0!expmu2/expnu/(robs**(three/two)+a_spin)*zero
        velocity(3)=exppsi/expnu*(one/(robs**(three/two)+a_spin)-somiga)  
        m=20
        n=20 
        r_in = rms(a_spin)
        rhorizon=one+sqrt(one-a_spin**two) 
        delta_theta=180.D0/n
        delta_phi=360.D0/m
        deltap=1.5D0/200.D0 
        open(unit=15,file='rayx2.txt',status="replace")
        open(unit=16,file='rayy2.txt',status="replace")          
        open(unit=17,file='rayz2.txt',status="replace")                                    
        DO i=0,m
                !write(*,*)i        
                phi_i=i*delta_phi*dtor
                !write(*,*)'phi=',phi
                Do j=0,n
                        theta=(j*delta_theta )*dtor
                        !write(*,*)'theta=',theta
                        pp= cos(theta)
                        pr= sin(theta)*sin(phi_i)
                        pt= sin(theta)*cos(phi_i)
                        call initialdirection(pr,pt,pp,sinobs,muobs,a_spin,robs,velocity,lambda,q,f1234)

                        ptotal = p_total(f1234(1),lambda,q,sinobs,muobs,a_spin,robs,scal)
                        Do k=0,200 
                            p=k*ptotal/200.D0 
                            call YNOGK(p,f1234,lambda,q,sinobs,muobs,a_spin,robs,scal,&
                            ra,mua,phya,timea,sigmaa) 
                            If(ra .ge. 1.2D0*rhorizon)then 
                                write(unit=15,fmt=*)sqrt(ra**two+a_spin**two)*sqrt(one-mua**two)*cos(phya)
                                write(unit=16,fmt=*)-sqrt(ra**two+a_spin**two)*sqrt(one-mua**two)*sin(phya)
                                write(unit=17,fmt=*)ra*mua 
                            Else
                                write(unit=15,fmt=*)sqrt(ra-100.D0) 
                                write(unit=16,fmt=*)sqrt(ra-100.D0) 
                                write(unit=17,fmt=*)sqrt(ra-100.D0)          
                            ENDIF  
                        Enddo
                Enddo
        Enddo
        close(unit=15)        
        close(unit=16)
        close(unit=17)
        return
      end subroutine reashift2










