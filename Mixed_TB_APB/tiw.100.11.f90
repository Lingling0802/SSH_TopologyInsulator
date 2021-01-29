!! modified from mph_ifort.f90
!! To make it "good" or "gd", following mistake was corrected.
!! Mistake in rmbadsxsy subroutine is corrected. See comments.
!! gfortran -o tiw.out  tiw.f90
!! "cotan(x)" is changed to "cos(x)/sin(x)" 
!! IMPORTANT: Cttx or Ctty are -i/sqrt(2.d0) times Ctx or Cty

PROGRAM tiw
  
  IMPLICIT NONE  
  
  INTEGER i,j,k,ki,kr,seed,isign,NDIM,N,imode,iset,ixy(6)
  INTEGER ix,iy,ixx,iyy
  REAL*8  PI,EPS
  REAL*8  A1,A2,A3,B,C1,C3,G1,G2,H1,H2,G1p,G2p
  
  REAL*8  dt,is
  INTEGER nt,ntfin,ntnoise,ntfnoise,ntprint
  
  PARAMETER(N=64)
  PARAMETER(NDIM=2)
  PARAMETER(EPS=1.d-10)
  INTEGER nn(NDIM)
  
  REAL*8 sx(N,N),sx_dat(2*N*N),sxt(N,N),sx_ini(N,N)
  REAL*8 sy(N,N),sy_dat(2*N*N),syt(N,N)
  REAL*8 p1(N,N),p1_dat(2*N*N),p1k0 ! p1=sx^2+sy^2
  REAL*8 p3(N,N),p3_dat(2*N*N),p3k0 ! p3=sx^2-sy^2
  REAL*8 s0,s_var,dnorm,e30,e10,simax,cshift,sdw
  REAL*8 sxmin(N,N),symin(N,N)
  REAL*8 sx_noise(N,N),sy_noise(N,N)

  REAL*8 e1(N,N),e2(N,N),e3(N,N),e1k0,e3k0
  REAL*8 dx(N,N),dy(N,N),d(N,N),th(N,N),Edist(N,N),Ddx(N,N),Ddy(N,N)

  REAL*8 sxtest(N,N),sytest(N,N),p1test(N,N),p3test(N,N)

  REAL*8 B11t_dat(2*N*N),B12t_dat(2*N*N),B22t_dat(2*N*N)
  REAL*8 B11tsx_dat(2*N*N),B12tsy_dat(2*N*N),B12tsx_dat(2*N*N),B22tsy_dat(2*N*N)
  REAL*8 B11tsx(N,N),B12tsy(N,N),B12tsx(N,N),B22tsy(N,N)
  REAL*8 Cttx_dat(2*N*N),Ctty_dat(2*N*N)
  REAL*8 Cttysx_dat(2*N*N),Cttxsy_dat(2*N*N)
  REAL*8 Cs1_dat(2*N*N),Cs1(N,N),Cs3_dat(2*N*N),Cs3(N,N)
  REAL*8 Qx1(N,N),Qy1(N,N),Qx1_dat(2*N*N),Qy1_dat(2*N*N)
  REAL*8 Qx3(N,N),Qy3(N,N),Qx3_dat(2*N*N),Qy3_dat(2*N*N)
  REAL*8 Cttyp1(N,N),Cttyp1_dat(2*N*N),Cttxp1(N,N),Cttxp1_dat(2*N*N)
  REAL*8 Cttyp3(N,N),Cttyp3_dat(2*N*N),Cttxp3(N,N),Cttxp3_dat(2*N*N)

  REAL*8 En1xy,En3xy,En1xy0,En3xy0,En13x0,En130y,En1300
  REAL*8 En1,En2,En3,Entot,Etemp,Entmin

  REAL*8 GE13xysx(N,N),GE13xysy(N,N),GE1300sx(N,N),GE1300sy(N,N)
  REAL*8 GE13x0sx(N,N),GE13x0sy(N,N),GE130ysx(N,N),GE130ysy(N,N)
  REAL*8 GE2sx(N,N),GE2sy(N,N)
  REAL*8 GEtotsx(N,N),GEtotsy(N,N)
  REAL*8 sxpipi,sypipi

  CHARACTER*20 fname

  PI=4.d0*atan(1.d0)
  nn(1)=N
  nn(2)=N

!===== set parameters, initial sx, sy, B11t, B12t, B22t
  open(50,file='tiw.in',status="unknown")
  read(50,*) A1
  read(50,*) A2
  read(50,*) A3
  read(50,*) B
  read(50,*) C1
  read(50,*) C3
  read(50,*) G1
  read(50,*) G2
  read(50,*) H1
  read(50,*) H2
  read(50,*) dt
  read(50,*) ntfin
  read(50,*) ntnoise
  read(50,*) ntfnoise
  read(50,*) ntprint
  read(50,*) s_var
  read(50,*) seed
  read(50,*) is
  read(50,*) dnorm
  read(50,*) simax
  read(50,*) cshift
  read(50,*) iset
  read(50,*) ixy(1),ixy(2),ixy(3),ixy(4)
!!  read(50,*) ixy(1),ixy(2),ixy(3),ixy(4),ixy(5),ixy(6)
  read(50,*) imode
  close(50)  
  
  G1p=G1-2.d0*C3**2/A3-2.d0*C1**2/A1
  G2p=G2+2.d0*C3**2/A3-2.d0*C1**2/A1

  if(G1p.lt.0.d0.and.(G1p**2-4.d0*B*H1).gt.0.d0) then
     s0=sqrt((-G1p+sqrt(G1p**2-4.d0*B*H1))/2.d0/H1)  
     e10=-C1/A1*s0**2
     e30=C3/A3*s0**2
  else
     s0=0.1d0
     e10=0.1d0
     e30=0.1d0
  end if

  sdw=0.327 !! PUT FORMULA HERE HERE HERE HERE XXXXXX !!!!!

  s_var=s_var*s0
  simax=simax*s0
  
!===== set initial sx,sy and read in Bt's and Ctt's
  open(511,file='B11t.dat',status="unknown")
  open(512,file='B12t.dat',status="unknown")
  open(522,file='B22t.dat',status="unknown")
  open(531,file='Cttx.dat',status="unknown")
  open(532,file='Ctty.dat',status="unknown")
  do k=1,2*N*N
     read(511,*)  B11t_dat(k) !  B11t_dat in k-space
     read(512,*)  B12t_dat(k) !  B12t_dat in k-space
     read(522,*)  B22t_dat(k) !  B22t_dat in k-space
     read(531,*)  Cttx_dat(k) !  Cttx_dat in k-space
     read(532,*)  Ctty_dat(k) !  Ctty_dat in k-space
  end do
  close(511)
  close(512)
  close(522)
  close(531)
  close(532)

  call sxsy_ini(sx,sy,s0,simax,N,seed,is,cshift,imode,iset,ixy) 

  if(imode.eq.19) then 
     simax=sdw
     call sxsy_ini(sx,sy,s0,simax,N,seed,is,cshift,imode,iset,ixy)
  end if

  sx_ini=sx

  do iy=N,1,-1
     do ix=1,N
        if(iy.eq.N.and.(ix.eq.1.or.ix.eq.2)) then
           write(14,*) (-1)**ix*simax  ! write these false data as a standard for color coding in xpm file
           write(15,*) (-1)**ix*simax  ! this does not affect sx,sy themselves in the simulation
        else
           write(14,*) sx(ix,iy)
           write(15,*) sy(ix,iy)
        end if
     end do
  end do 

!===== find and write initial configuration
  call find_all(sx,sy,e1,e2,e3,p1,p3,dx,dy,e1k0,e3k0,d,th,Edist,A1,A2,A3,B,C1,C3,G1,G2,H1,H2,N)

  do iy=N,1,-1
     do ix=1,N
        sxt(ix,iy)=sx(ix,iy)*(-1)**(ix+iy)
        syt(ix,iy)=sy(ix,iy)*(-1)**(ix+iy)
     end do
  end do

  e1(1,N)=simax    ! set these in output to help understanding color coding in xpm file
  e1(2,N)=-simax   ! these false data do not affect anything in the simulation
  e2(1,N)=simax    ! But, don't make this change in sx and sy
  e2(2,N)=-simax
  e3(1,N)=simax
  e3(2,N)=-simax
  sxt(1,N)=simax
  sxt(2,N)=-simax
  syt(1,N)=simax
  syt(2,N)=-simax
  p1(1,N)=simax**2
  p1(2,N)=-simax**2
  p3(1,N)=simax**2
  p3(2,N)=-simax**2

  do iy=N,1,-1
     do ix=1,N
        write(11,*) e1(ix,iy)
        write(12,*) e2(ix,iy)
        write(13,*) e3(ix,iy)
        write(16,*) sxt(ix,iy)
        write(17,*) syt(ix,iy)
        write(18,*) p1(ix,iy)
        write(19,*) p3(ix,iy)
        write(20,*) Edist(ix,iy)
        write(61,*) dx(ix,iy)
        write(62,*) dy(ix,iy)
     end do
  end do


!********************** START LOOP ************************  
  Entmin=1.d0/EPS
  DO 999 nt=1,ntfin

!==== Exit if ntfin is zero
     if (ntfin.eq.0) exit


!===== add Gaussian noise
     if ((mod(nt,ntnoise).eq.0).and.(nt.lt.ntfnoise)) then
        call gauss2(seed,sx_noise,N,s_var)
        call gauss2(seed,sy_noise,N,s_var)
	call gauss2(seed,sx_noise,N,s_var)  ! For unknown reason, sx_noise should be called again. 
        sx=sx+sx_noise
        sy=sy+sy_noise
        call rmbadsxsy(sx,sy,N)
     end if

     if(nt.eq.ntfnoise) then
        sx=sxmin ! After last noise, use min conf. as initial      
        sy=symin
     end if

!===== Find p1 & p3
     p1=sx*sx+sy*sy
     p3=sx*sx-sy*sy

!===== Fourier transform sx, sy, p1, p3
     do iy=1,N
        do ix=1,N          
           k=2*(ix+N*(iy-1))-1
           sx_dat(k)=sx(ix,iy)/dble(N*N)   ! sx_dat in r-space
           sx_dat(k+1)=0.d0
           sy_dat(k)=sy(ix,iy)/dble(N*N)   ! sy_dat in r-space
           sy_dat(k+1)=0.d0
           p1_dat(k)=p1(ix,iy)/dble(N*N)   ! p1_dat in r-space
           p1_dat(k+1)=0.d0
           p3_dat(k)=p3(ix,iy)/dble(N*N)   ! p3_dat in r-space
           p3_dat(k+1)=0.d0
        end do
     end do
     
     isign=+1
     call fourn(sx_dat,nn,NDIM,isign) ! gives sx_dat in k-space
     call fourn(sy_dat,nn,NDIM,isign) ! gives sy_dat in k-space
     call fourn(p1_dat,nn,NDIM,isign) ! gives p1_dat in k-space
     call fourn(p3_dat,nn,NDIM,isign) ! gives p3_dat in k-space

     p1k0=p1_dat(1)  
     p3k0=p3_dat(1)


!*** Below this line, NEVER change sx_dat,sy_dat,p1_dat,p3_dat into r-space *******

!===== calculate B11tsx=B11t*sx,B12tsy=B12t*sy,B12tsx=B12t*sx,B22tsy=B22t*sy
     call mulink(N,sx_dat,B11t_dat,B11tsx_dat)   ! B11tsx_dat in k-space  
     call mulink(N,sy_dat,B12t_dat,B12tsy_dat)   ! B12tsy_dat in k-space  
     call mulink(N,sx_dat,B12t_dat,B12tsx_dat)   ! B12tsx_dat in k-space  
     call mulink(N,sy_dat,B22t_dat,B22tsy_dat)   ! B22tsy_dat in k-space  
  
     isign=-1     
     call fourn(B11tsx_dat,nn,NDIM,isign) ! gives B11tsx_dat in r-space
     call fourn(B12tsy_dat,nn,NDIM,isign) ! gives B12tsy_dat in r-space
     call fourn(B12tsx_dat,nn,NDIM,isign) ! gives B12tsx_dat in r-space
     call fourn(B22tsy_dat,nn,NDIM,isign) ! gives B22tsy_dat in r-space
  
     do iy=1,N
        do ix=1,N          
           k=2*(ix+N*(iy-1))-1
           B11tsx(ix,iy)=B11tsx_dat(k)  ! gives B11tsx
           B12tsy(ix,iy)=B12tsy_dat(k)  ! gives B12tsy
           B12tsx(ix,iy)=B12tsx_dat(k)  ! gives B12tsx
           B22tsy(ix,iy)=B22tsy_dat(k)  ! gives B22tsy
        end do
     end do
     
!===== calculate Cttysx=Ctty*sx,Cttxsy=Cttx*sy
     call mulink(N,sx_dat,Ctty_dat,Cttysx_dat)   ! Cttysx_dat in k-space  
     call mulink(N,sy_dat,Cttx_dat,Cttxsy_dat)   ! Cttxsy_dat in k-space  

     Cs1_dat=Cttysx_dat+Cttxsy_dat    ! Cs1_dat in k-space
     Cs3_dat=Cttysx_dat-Cttxsy_dat    ! Cs3_dat in k-space
     
     isign=-1     
     call fourn(Cs1_dat,nn,NDIM,isign) ! gives Cs1_dat in r-space
     call fourn(Cs3_dat,nn,NDIM,isign) ! gives Cs3_dat in r-space

     do iy=1,N
        do ix=1,N          
           k=2*(ix+N*(iy-1))-1
           Cs1(ix,iy)=Cs1_dat(k)  ! gives Cs1
           Cs3(ix,iy)=Cs3_dat(k)  ! gives Cs3
        end do
     end do

!===== calculate Qx1, Qx3, Qy1 and Qy3
     do iy=0,N-1     !  ky=PI*2.d0*dble(iy)/dble(N)     
        do ix=0,N-1  !  kx=PI*2.d0*dble(ix)/dble(N)
           k=(iy*N+ix)*2+1
           if(ix.ne.0.and.iy.eq.0) then
              Qx1_dat(k)=p1_dat(k)         ! Qx1_dat in k-space
              Qx1_dat(k+1)=p1_dat(k+1)
              Qx3_dat(k)=p3_dat(k)         ! Qx3_dat in k-space
              Qx3_dat(k+1)=p3_dat(k+1)
           else
              Qx1_dat(k)=0.d0
              Qx1_dat(k+1)=0.d0
              Qx3_dat(k)=0.d0
              Qx3_dat(k+1)=0.d0
           end if
        end do
     end do
     
     do iy=0,N-1     !  ky=PI*2.d0*dble(iy)/dble(N)     
        do ix=0,N-1  !  kx=PI*2.d0*dble(ix)/dble(N)
           k=(iy*N+ix)*2+1
           if(ix.eq.0.and.iy.ne.0) then
              Qy1_dat(k)=p1_dat(k)        ! Qy1_dat in k-space
              Qy1_dat(k+1)=p1_dat(k+1)
              Qy3_dat(k)=p3_dat(k)        ! Qy3_dat in k-space
              Qy3_dat(k+1)=p3_dat(k+1)
           else
              Qy1_dat(k)=0.d0
              Qy1_dat(k+1)=0.d0
              Qy3_dat(k)=0.d0
              Qy3_dat(k+1)=0.d0
           end if
        end do
     end do
     
     isign=-1     
     call fourn(Qx1_dat,nn,NDIM,isign) ! gives Qx1_dat in r-space
     call fourn(Qy1_dat,nn,NDIM,isign) ! gives Qy1_dat in r-space
     call fourn(Qx3_dat,nn,NDIM,isign) ! gives Qx3_dat in r-space
     call fourn(Qy3_dat,nn,NDIM,isign) ! gives Qy3_dat in r-space
     
     do iy=1,N
        do ix=1,N          
           k=2*(ix+N*(iy-1))-1
           Qx1(ix,iy)=Qx1_dat(k)  ! gives Qx1
           Qy1(ix,iy)=Qy1_dat(k)  ! gives Qy1
           Qx3(ix,iy)=Qx3_dat(k)  ! gives Qx3
           Qy3(ix,iy)=Qy3_dat(k)  ! gives Qy3
        end do
     end do
     
!===== calculate Cttyp1=Ctty*p1, Cttxp1=Cttx*p1, Cttyp3=Ctty*p3 and Cttxp3=Cttx*p3
     call mulink(N,p1_dat,Ctty_dat,Cttyp1_dat)   ! Cttyp1_dat in k-space  
     call mulink(N,p1_dat,Cttx_dat,Cttxp1_dat)   ! Cttxp1_dat in k-space  
     call mulink(N,p3_dat,Ctty_dat,Cttyp3_dat)   ! Cttyp3_dat in k-space  
     call mulink(N,p3_dat,Cttx_dat,Cttxp3_dat)   ! Cttxp3_dat in k-space  
     
     isign=-1     
     call fourn(Cttyp1_dat,nn,NDIM,isign) ! gives Cttyp1_dat in r-space
     call fourn(Cttxp1_dat,nn,NDIM,isign) ! gives Cttxp1_dat in r-space
     call fourn(Cttyp3_dat,nn,NDIM,isign) ! gives Cttyp3_dat in r-space
     call fourn(Cttxp3_dat,nn,NDIM,isign) ! gives Cttxp3_dat in r-space
     
     do iy=1,N
        do ix=1,N          
           k=2*(ix+N*(iy-1))-1
           Cttyp1(ix,iy)=Cttyp1_dat(k)  ! gives Cttyp1
           Cttxp1(ix,iy)=Cttxp1_dat(k)  ! gives Cttxp1
           Cttyp3(ix,iy)=Cttyp3_dat(k)  ! gives Cttyp3
           Cttxp3(ix,iy)=Cttxp3_dat(k)  ! gives Cttxp3
        end do
     end do
     
!===== calculate energy per site, i.e., total energy devided by # of sites
     En1xy=sum((sx*B11tsx+sx*B12tsy+sy*B12tsx+sy*B22tsy)/2.d0)/dble(N*N)
     En3xy=sum(C1*p1*Cs1+C3*p3*Cs3)/dble(N*N) 
     En13x0=sum(-0.5d0/(A1+A3)*(C3*Qx3+C1*Qx1)**2)/dble(N*N)
     En130y=sum(-0.5d0/(A1+A3)*(C3*Qy3-C1*Qy1)**2)/dble(N*N)
     En1300=-0.5d0*C3**2/A3*p3k0*p3k0-0.5d0*C1**2/A1*p1k0*p1k0     
     En2=sum(0.5d0*B*(sx**2+sy**2)+0.25d0*G1*(sx**4+sy**4)+0.5d0*G2*sx**2*sy**2   &
         +1.d0/6.d0*H1*(sx**6+sy**6)+1.d0/6.d0*H2*(sx**2+sy**2)*sx**2*sy**2)/dble(N*N)
     
     Entot=En1xy+En3xy+En13x0+En130y+En1300+En2
     
     if(mod(nt,ntprint).eq.0) then
        write(*,*) nt,Entot 
	write(40,*) nt,sx(1,iset-1)
     end if

!===== save minimum total energy and sx
     if(Entot.lt.Entmin) then
        Entmin=Entot
        sxmin=sx
        symin=sy
     end if

!===== calculate gradient of Etot about sx and sy
     GE13xysx=B11tsx+B12tsy+2.d0*C1*sx*Cs1-C1*Cttyp1+2.d0*C3*sx*Cs3-C3*Cttyp3  !----- grad_sx(E1+E3)_(kx.ne.0 and ky.ne.0)     
     GE13xysy=B12tsx+B22tsy+2.d0*C1*sy*Cs1-C1*Cttxp1-2.d0*C3*sy*Cs3+C3*Cttxp3  !----- grad_sy(E1+E3)_(kx.ne.0 and ky.ne.0)     
     
     GE13x0sx=-2.d0*(C1+C3)/(A1+A3)*sx*(C1*Qx1+C3*Qx3)  !----- grad_sx(E1+E3)_(kx.ne.0 and ky.eq.0)
     GE13x0sy=-2.d0*(C1-C3)/(A1+A3)*sy*(C1*Qx1+C3*Qx3)  !----- grad_sy(E1+E3)_(kx.ne.0 and ky.eq.0)
     
     GE130ysx=-2.d0*(C1-C3)/(A1+A3)*sx*(C1*Qy1-C3*Qy3)  !----- grad_sx(E1+E3)_(kx.eq.0 and ky.ne.0)
     GE130ysy=-2.d0*(C1+C3)/(A1+A3)*sy*(C1*Qy1-C3*Qy3)  !----- grad_sy(E1+E3)_(kx.eq.0 and ky.ne.0)
     
     GE1300sx=-2.d0*C1**2/A1*sx*p1k0-2.d0*C3**2/A3*sx*p3k0   !----- grad_sx(E1+E3)_(kx.eq.0 and ky.eq.0)
     GE1300sy=-2.d0*C1**2/A1*sy*p1k0+2.d0*C3**2/A3*sy*p3k0   !----- grad_sy(E1+E3)_(kx.eq.0 and ky.eq.0)
     
     GE2sx=B*sx+G1*sx**3+G2*sx*sy**2+H1*sx**5+1.d0/3.d0*H2*(2.d0*sx**2+sy**2)*sx*sy**2 !---- grad_sx(E2)
     GE2sy=B*sy+G1*sy**3+G2*sy*sx**2+H1*sy**5+1.d0/3.d0*H2*(2.d0*sy**2+sx**2)*sy*sx**2 !---- grad_sy(E2)
     
     GEtotsx=GE13xysx+GE13x0sx+GE130ysx+GE1300sx+GE2sx  !---  gives GEtotsx=grad_sx(Etot)
     GEtotsy=GE13xysy+GE13x0sy+GE130ysy+GE1300sy+GE2sy  !---  gives GEtotsy=grad_sy(Etot)

!===== find sx and sy for next step
 
     sx=sx-dt*GEtotsx
     sy=sy-dt*GEtotsy
  
!--- This is the extra step to find constrained minimum energy configuartations.
!--- Selected sx(ix,iy)'s are unchanged from initial values, except by the rmbadsxsy subroutine.
!--- This extra step can be turned off by choosing iset > N+N
!--- This step is used, particularly, for imode=13.
     do iy=N,1,-1
	do ix=1,N
	   if(mod(ix+iy,N).eq.iset) then
              sx(ix,iy)=sx_ini(ix,iy)
           end if
	end do
     end do

     call rmbadsxsy(sx,sy,N) 
!--- Due to the above rmbadsxsy step, sx(ix,iy)'s along a line set by iset are 
!--- still changed from sx_ini(ix,iy) values.
!--- But they asymptotically approach to the sx_ini(ix,iy) value. 

999  CONTINUE

!********************** END LOOP ************************  

  sx=sxmin
  sy=symin

!===== find e1,e2,e3,dx,dy
  call find_all(sx,sy,e1,e2,e3,p1,p3,dx,dy,e1k0,e3k0,d,th,Edist,A1,A2,A3,B,C1,C3,G1,G2,H1,H2,N)

  do iy=N,1,-1
     do ix=1,N
        sxt(ix,iy)=sx(ix,iy)*(-1)**(ix+iy)
        syt(ix,iy)=sy(ix,iy)*(-1)**(ix+iy)
     end do
  end do

!===== calculate and write inputs for fle.f90 
     do ix=1,N
        do iy=1,N
           
           if(ix.eq.N) then
              ixx=0
           else
              ixx=ix
           end if
           if(iy.eq.N) then
              iyy=0
           else
              iyy=iy
           end if
           
           ! We can use the following formula since e2(k=0)=0
           Ddx(ix,iy)=dx(ixx+1,iy)-dx(ix,iy)+(e1k0+e3k0)/sqrt(2.d0)
           Ddy(ix,iy)=dy(ix,iyy+1)-dy(ix,iy)+(e1k0-e3k0)/sqrt(2.d0)
           
        end do
     end do

     do iy=1,N
        do ix=1,N
           write(99,203) Ddx(ix,iy),Ddy(ix,iy)  ! write results for fle.f90
203        format(2f8.3)
        end do
     end do

!===== Write the final results for sx and sy without color standards.
!===== If imode=10, then these outputs are saved as sx_ini.in and sy_ini.in
!===== and are used for the next run.
     do iy=N,1,-1
      	do ix=1,N
         write(34,*) sx(ix,iy) 
         write(35,*) sy(ix,iy)  
      	end do
     end do

 if(imode.eq.10) then
    open(250,file='sx_ini.gd.in',status="unknown")
    open(251,file='sy_ini.gd.in',status="unknown")
    do iy=N,1,-1
      do ix=1,N
         write(250,*) sx(ix,iy)
         write(251,*) sy(ix,iy)
      end do
    end do
    close(250)
    close(251)
  end if

!===== write final results without color standards
  do iy=N,1,-1
     do ix=1,N
        write(44,*) sx(ix,iy)
        write(45,*) sy(ix,iy)  
        write(46,*) dx(ix,iy) ! Write displacement output for displacement-based simulations
        write(47,*) dy(ix,iy) ! Write displacement output for displacement-based simulations
        write(30,*) Edist(ix,iy)
     end do
  end do

  do ix=1,N
     do iy=1,N
        write(31,202) ix,iy,dx(ix,iy),dy(ix,iy)
202     format(2i4,2f8.3)
     end do
  end do

!===== write the final results with color standards
  e1(1,N)=s0
  e1(2,N)=-s0
  e2(1,N)=s0
  e2(2,N)=-s0
  e3(1,N)=s0
  e3(2,N)=-s0
  sx(1,N)=s0  ! set these in output to help understanding color coding in xpm file
  sx(2,N)=-s0 ! these false data do not affect anything in the simulation
  sy(1,N)=s0
  sy(2,N)=-s0
  sxt(1,N)=s0
  sxt(2,N)=-s0
  syt(1,N)=s0
  syt(2,N)=-s0
  p1(1,N)=s0**2
  p1(2,N)=-s0**2
  p3(1,N)=s0**2
  p3(2,N)=-s0**2

  do iy=N,1,-1
     do ix=1,N
        write(21,*) e1(ix,iy)
        write(22,*) e2(ix,iy)
        write(23,*) e3(ix,iy)
        write(24,*) sx(ix,iy)
        write(25,*) sy(ix,iy)  
        write(26,*) sxt(ix,iy)
        write(27,*) syt(ix,iy)
        write(28,*) p1(ix,iy)
        write(29,*) p3(ix,iy)
     end do
  end do

!----- test whether dx and dy were calculated correctly

  do ix=1,N
     do iy=1,N
        if(ix.ne.N.and.iy.ne.N) then
           sxtest(ix,iy)=(dx(ix,iy)-dx(ix+1,iy)-dx(ix,iy+1)+dx(ix+1,iy+1))/2.d0
           sytest(ix,iy)=(dy(ix,iy)-dy(ix+1,iy)-dy(ix,iy+1)+dy(ix+1,iy+1))/2.d0
        else if(ix.eq.N.and.iy.ne.N) then
           sxtest(ix,iy)=(dx(ix,iy)-dx(1,iy)-dx(ix,iy+1)+dx(1,iy+1))/2.d0
           sytest(ix,iy)=(dy(ix,iy)-dy(1,iy)-dy(ix,iy+1)+dy(1,iy+1))/2.d0
        else if(ix.ne.N.and.iy.eq.N) then
           sxtest(ix,iy)=(dx(ix,iy)-dx(ix+1,iy)-dx(ix,1)+dx(ix+1,1))/2.d0
           sytest(ix,iy)=(dy(ix,iy)-dy(ix+1,iy)-dy(ix,1)+dy(ix+1,1))/2.d0
        else 
           sxtest(ix,iy)=(dx(ix,iy)-dx(1,iy)-dx(ix,1)+dx(1,1))/2.d0
           sytest(ix,iy)=(dy(ix,iy)-dy(1,iy)-dy(ix,1)+dy(1,1))/2.d0
        end if
     end do
  end do

  p1test=sxtest**2+sytest**2
  p3test=sxtest**2-sytest**2

  sxtest(1,N)=1.0*s0  ! set these in output to help understanding color coding in xpm file
  sxtest(2,N)=-1.0*s0 
  sytest(1,N)=1.0*s0
  sytest(2,N)=-1.0*s0
  p1test(1,N)=1.0*s0**2
  p1test(2,N)=-1.0*s0**2
  p3test(1,N)=1.0*s0**2
  p3test(2,N)=-1.0*s0**2

  do iy=N,1,-1
     do ix=1,N
        write(117,*) sxtest(ix,iy)
        write(118,*) sytest(ix,iy)  
        write(119,*) p1test(ix,iy)
        write(120,*) p3test(ix,iy)
     end do
  end do

!----- test ends here

  fname='vec.d.ps'
  call pswrite(fname,N,d,th,dnorm)

!===== write input parameters  
  write(10,*) '##Input parameters:'
  write(10,*) 'A1=',A1,':A2=',A2,':A3=',A3,':B=',B,':C1=',C1,':C3=',C3
  write(10,*) 'G1=',G1,':G2=',G2,':H1=',H1,':H2=',H2
  write(10,204) G1p,G2p,G1p**2-4.d0*B*H1
204 format(' G1p=',f9.3,':G2p=',f9.3,':G1p**2-4.d0*B*H1=',f9.3)
  write(10,*) 'N=',N
  write(10,*) 'imode=',imode
  write(10,*) 'seed=',seed
  write(10,*) 'dt=',dt
  write(10,*) 'ntfin=',ntfin
  write(10,*) 'ntnoise=',ntnoise
  write(10,*) 'ntfnoise=',ntfnoise
  write(10,*) 'ntprint=',ntprint
  write(10,*) 's_var=',s_var 
  write(10,*) 'simax/s0=',simax/s0
  write(10,*) 'cshift=',cshift
  write(10,*) 'is=',is
  write(10,*) 'iset=',iset,'(un-set, if > N+N)'
  write(10,*) 'ixy=',ixy(1),',',ixy(2),',',ixy(3),',',ixy(4)
  write(10,*) ' '
  write(10,*) '## In homogenious state:'
  write(10,*) '(sx_min) =',s0
  write(10,*) '(e1_min) =',e10
  write(10,*) '(e3_min) =',e30
  write(10,*) ' '
  write(10,*) '## Final results:'
  write(10,*) 'Emin=',Entot
!-- The above energy, Entot, is normalized per site.
  write(10,*) 'Emin*N*N=',Entot*dble(N*N)
  write(10,*) 'e1(k=0)=',e1k0
  write(10,*) 'e3(k=0)=',e3k0

  STOP      
  
END PROGRAM tiw

!******** Main program ends here *************!

!******** Subroutines ************************!
SUBROUTINE sxsy_ini(sx,sy,s0,simax,N,seed,is,cshift,imode,iset,ixy)
  INTEGER N,seed,ixy(6)
  REAL*8 is,cshift,temp
  REAL*8 PI
  REAL*8 sx(N,N),sy(N,N),simax,s0
  REAL*8 ran1

  PI=4.d0*atan(1.d0)

!*** imode=0: homogeneous ground state
  if(imode.eq.0) then
     do iy=1,N
        do ix=1,N
           sx(ix,iy)=simax*(-1)**(ix+iy)
           sy(ix,iy)=0
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=1: random
  if(imode.eq.1) then
     do iy=1,N
        do ix=1,N
           sx(ix,iy)=(ran1(seed)-0.5)*2.0*simax
           sy(ix,iy)=(ran1(seed)-0.5)*2.0*simax
           temp=abs(sx(ix,iy)**2-sy(ix,iy)**2)
           if(temp.gt.simax**2) write(*,*) temp,ix,iy,imode
        end do
     end do
     call rmbadsxsy(sx,sy,N)
     do iy=1,N
        do ix=1,N
           temp=abs(sx(ix,iy)**2-sy(ix,iy)**2)
           if(temp.gt.simax**2) write(*,*) temp,ix,iy
        end do
     end do

  end if

!*** imode=2: start with TB along 135 degree with no phase shift
  if(imode.eq.2) then
     do iy=1,N
        do ix=1,N
           if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).gt.0.d0) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).lt.0.d0) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=simax*(-1)**(ix+iy)
           else
              sx(ix,iy)=0.d0
              sy(ix,iy)=0.d0   
           end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=3: start with TB along 135 degree with phase shift
  if(imode.eq.3) then
     do iy=1,N
        do ix=1,N
           if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).gt.0.d0) then
              sx(ix,iy)=simax*(-1)**(ix+iy+1)
              sy(ix,iy)=0.d0
           else if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).lt.0.d0) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=simax*(-1)**(ix+iy)
           else
              sx(ix,iy)=0.d0
              sy(ix,iy)=0.d0   
           end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=4: start with cos TB along 135 degree with no phase shift
  if(imode.eq.4) then
     do iy=1,N
        do ix=1,N
           if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).gt.0.d0) then
              sx(ix,iy)=simax*(-1)**(ix+iy)*cos((dble(ix+iy)-is)/dble(N)*2.d0*PI)
              sy(ix,iy)=0.d0
           else if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).lt.0.d0) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)*cos((dble(ix+iy)-is)/dble(N)*2.d0*PI)
           else
              sx(ix,iy)=0.d0
              sy(ix,iy)=0.d0   
           end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=5: start with cos TB along 135 degree with phase shift
  if(imode.eq.5) then
     do iy=1,N
        do ix=1,N
           if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).gt.0.d0) then
              sx(ix,iy)=simax*(-1)**(ix+iy+1)*cos((dble(ix+iy)-is)/dble(N)*2.d0*PI)
              sy(ix,iy)=0.d0
           else if(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI).lt.0.d0) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)*cos((dble(ix+iy)-is)/dble(N)*2.d0*PI)
           else
              sx(ix,iy)=0.d0
              sy(ix,iy)=0.d0   
           end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=6: start with cos APB along 90 degree with phase shift
  if(imode.eq.6) then
     do iy=1,N
        do ix=1,N
           sx(ix,iy)=simax*(-1)**(ix+iy+1)*cos((dble(ix)-is)/dble(N)*2.d0*PI)
           sy(ix,iy)=0.d0
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=7: random sx, zero sy
  if(imode.eq.7) then
     do iy=1,N
        do ix=1,N
           sx(ix,iy)=(ran1(seed)-0.5)*2.0*simax
           sy(ix,iy)=0.d0
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=8: random+emboss
  if(imode.eq.8) then
     do iy=1,N
        do ix=1,N
           sx(ix,iy)=(ran1(seed)-0.5)*2.0*simax
           sy(ix,iy)=(ran1(seed)-0.5)*2.0*simax
        end do
     end do

     do i=0,2
        sx(i+1,i+1)=+simax*1.d0
        sy(i+1,i+1)=0.d0
        sx(i+2,i+1)=-simax*1.d0
        sy(i+2,i+1)=0.d0
        sx(i+3,i+1)=0.d0
        sy(i+3,i+1)=-simax*1.d0
        sx(i+4,i+1)=0.d0
        sy(i+4,i+1)=+simax*1.d0

        sx(i+5,i+1)=+simax*1.d0
        sy(i+5,i+1)=0.d0
        sx(i+6,i+1)=-simax*1.d0
        sy(i+6,i+1)=0.d0
        sx(i+7,i+1)=0.d0
        sy(i+7,i+1)=-simax*1.d0
        sx(i+8,i+1)=0.d0
        sy(i+8,i+1)=+simax*1.d0
     end do

     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=9: start with cos of sx along 135 deg
  if(imode.eq.9) then
     do iy=1,N
        do ix=1,N
              sx(ix,iy)=simax*(-1)**(ix+iy)*(cos((dble(ix+iy)-is)/dble(N)*2.d0*PI)+cshift)
              sy(ix,iy)=0.d0
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=10: read input for sx and sy from files named sx_ini.in and sy_ini.in,
!************* and then change sx along a line set by mod(ix+iy,N)=iset by 0.1*s0
!************* Note that this sx is changed further by rmbadsxsy, of course.
  if(imode.eq.10) then
     open(520,file='sx_ini.in',status="unknown")
     open(521,file='sy_ini.in',status="unknown")	
        do iy=N,1,-1
        do ix=1,N
           read(520,*) sx(ix,iy)
           read(521,*) sy(ix,iy)
        end do
        end do
     close(520)  
     close(521)
     call rmbadsxsy(sx,sy,N)

     do iy=N,1,-1
	do ix=1,N
	   if(mod(ix+iy,N).eq.iset) then
              sx(ix,iy)=sx(ix,iy)+0.1*s0
           end if
	end do
     end do
    call rmbadsxsy(sx,sy,N)
  end if

!*** imode=11: two step-wise sx=s0 strpes  along 135 degree
  if(imode.eq.11) then
     do iy=1,N
        do ix=1,N
	  if(mod(ix+iy,64).ge.ixy(1).and.mod(ix+iy,64).le.ixy(2)) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
	  else if(mod(ix+iy,64).ge.ixy(3).and.mod(ix+iy,64).le.ixy(4)) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
          else
              sx(ix,iy)=0.d0
              sy(ix,iy)=0.d0
	  end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=12: two step-wise sx=s0 strpes with different phase along 135 degree
  if(imode.eq.12) then
     do iy=1,N
        do ix=1,N
	  if(mod(ix+iy,64).ge.ixy(1).and.mod(ix+iy,64).le.ixy(2)) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
	  else if(mod(ix+iy,64).ge.ixy(3).and.mod(ix+iy,64).le.ixy(4)) then
              sx(ix,iy)=simax*(-1)**(ix+iy+1)
              sy(ix,iy)=0.d0
          else
              sx(ix,iy)=0.d0
              sy(ix,iy)=0.d0
	  end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=13: read input for sx and sy from files named sx_stable.gd.in and sy_stable.gd.in,
!************* and then set sx distortion in a region 
  if(imode.eq.13) then
     open(520,file='sx_stable.gd.in',status="unknown")
     open(521,file='sy_stable.gd.in',status="unknown")	
        do iy=N,1,-1
        do ix=1,N
           read(520,*) sx(ix,iy)
           read(521,*) sy(ix,iy)
        end do
        end do
     close(520)  
     close(521)
     call rmbadsxsy(sx,sy,N)

     do iy=N,1,-1
     do ix=1,N
        if(ix+iy.ge.ixy(1).and.ix+iy.le.ixy(2).and.ix-iy.ge.ixy(3).and.ix-iy.le.ixy(4)) then
           sx(ix,iy)=s0*(-1)**(ix+iy)
        end if
     end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=14: read input for sx and sy from files named sx_stable.in and sy_stable.in,
!************* and then remove sx distortion in a region 
  if(imode.eq.14) then
     open(520,file='sx_stable.gd.in',status="unknown")
     open(521,file='sy_stable.gd.in',status="unknown")	
        do iy=N,1,-1
        do ix=1,N
           read(520,*) sx(ix,iy)
           read(521,*) sy(ix,iy)
        end do
        end do
     close(520)  
     close(521)
     call rmbadsxsy(sx,sy,N)

     do iy=N,1,-1
     do ix=1,N
        if(ix+iy.ge.ixy(1).and.ix+iy.le.ixy(2).and.ix-iy.ge.ixy(3).and.ix-iy.le.ixy(4)) then
           sx(ix,iy)=0.d0
        end if
     end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=15: read input for sx and sy from files named sx_stable.gd.in and sy_stable.gd.in,
!************* and then set sx distortion in two regions 
  if(imode.eq.15) then
     open(520,file='sx_stable.gd.in',status="unknown")
     open(521,file='sy_stable.gd.in',status="unknown")	
        do iy=N,1,-1
        do ix=1,N
           read(520,*) sx(ix,iy)
           read(521,*) sy(ix,iy)
        end do
        end do
     close(520)  
     close(521)
     call rmbadsxsy(sx,sy,N)

     do iy=N,1,-1
     do ix=1,N
        if(ix+iy.ge.ixy(1).and.ix+iy.le.ixy(2).and.ix-iy.ge.ixy(3).and.ix-iy.le.ixy(4)) then
           sx(ix,iy)=s0*(-1)**(ix+iy)
        end if
        if(ix+iy.ge.ixy(1).and.ix+iy.le.ixy(2).and.ix-iy.ge.ixy(5).and.ix-iy.le.ixy(6)) then
           sx(ix,iy)=s0*(-1)**(ix+iy)
        end if
     end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=16: start with cos APB along 0 degree with phase shift
  if(imode.eq.16) then
     do iy=1,N
        do ix=1,N
           sx(ix,iy)=simax*(-1)**(ix+iy+1)*cos((dble(iy)-is)/dble(N)*2.d0*PI)
           sy(ix,iy)=0.d0
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=17: read input for sx and sy from files named sx.input and sy.input,
!************* and then change signs of sx and sy 
  if(imode.eq.17) then
     open(520,file='sx.input',status="unknown")
     open(521,file='sy.input',status="unknown")	
        do iy=N,1,-1
        do ix=1,N
           read(520,*) sx(ix,iy)
           read(521,*) sy(ix,iy)
        end do
        end do
     close(520)  
     close(521)

     do iy=N,1,-1
     do ix=1,N
        sx(ix,iy)=-sx(ix,iy)
        sy(ix,iy)=-sy(ix,iy)
     end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=18: read input for sx and sy from files named sx.input and sy.input,
!************* and then change signs of pi,pi component of sx and sy 
  if(imode.eq.18) then
     open(520,file='sx.input',status="unknown")
     open(521,file='sy.input',status="unknown")	
        do iy=N,1,-1
        do ix=1,N
           read(520,*) sx(ix,iy)
           read(521,*) sy(ix,iy)
        end do
        end do
     close(520)  
     close(521)

     sxpipi=0.d0
     sypipi=0.d0

     do iy=N,1,-1
     do ix=1,N
        sxpipi=sxpipi+(-1)**(ix+iy)*sx(ix,iy)/dble(N*N)
        sypipi=sypipi+(-1)**(ix+iy)*sy(ix,iy)/dble(N*N)
     end do
     end do

     do iy=N,1,-1
     do ix=1,N
        sx(ix,iy)=sx(ix,iy)-1.d0*sxpipi*(-1)**(ix+iy)
        sy(ix,iy)=sy(ix,iy)-1.d0*sypipi*(-1)**(ix+iy)
     end do
     end do

     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=19: TB along 135 degree and APB of sx within e3 < 0 domain: Pattern #1
  if(imode.eq.19) then
     do iy=1,N
        do ix=1,N
           if(ix+iy.le.N/4+1) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)
           else if(ix+iy.le.3*N/4+1.and.iy-ix.lt.-N/4) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.3*N/4+1.and.iy-ix.ge.-N/4.and.iy-ix.le.N/4) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.3*N/4+1.and.iy-ix.gt.N/4) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.5*N/4) then 
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)
   	   else if(ix+iy.le.7*N/4+1.and.iy-ix.lt.-N/4) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.7*N/4+1.and.iy-ix.ge.-N/4.and.iy-ix.le.N/4) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.7*N/4+1.and.iy-ix.gt.N/4) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.2*N) then 
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)
           end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

!*** imode=20: TB along 135 degree and APB of sx within e3 < 0 domain: Pattern #2
  if(imode.eq.20) then
     do iy=1,N
        do ix=1,N
           if(ix+iy.le.N/4+1) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.3*N/4+1) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)
   	   else if(ix+iy.le.5*N/4+1.and.iy-ix.lt.-N/4) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.5*N/4+1.and.iy-ix.ge.-N/4.and.iy-ix.le.N/4) then
              sx(ix,iy)=simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.5*N/4+1.and.iy-ix.gt.N/4) then
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           else if(ix+iy.le.7*N/4+1) then
              sx(ix,iy)=0.d0
              sy(ix,iy)=-simax*(-1)**(ix+iy)
           else if(ix+iy.le.2*N) then 
              sx(ix,iy)=-simax*(-1)**(ix+iy)
              sy(ix,iy)=0.d0
           end if
        end do
     end do
     call rmbadsxsy(sx,sy,N)
  end if

END SUBROUTINE sxsy_ini

SUBROUTINE rmbadsxsy(sx,sy,N)
  INTEGER N,ix,iy,ixp,iyp
  REAL*8 sx(N,N),sx_bad(N,N),badsx_x(N),badsx_y(N),badsx_0
  REAL*8 sy(N,N),sy_bad(N,N),badsy_x(N),badsy_y(N),badsy_0

  do ix=1,N
     badsx_x(ix)=0.d0
     badsy_x(ix)=0.d0
     do iyp=1,N
        badsx_x(ix)=badsx_x(ix)+sx(ix,iyp) !-- ky=0 part
        badsy_x(ix)=badsy_x(ix)+sy(ix,iyp) !-- ky=0 part
     end do
     badsx_x(ix)=badsx_x(ix)/dble(N)
     badsy_x(ix)=badsy_x(ix)/dble(N)
  end do
  
  do iy=1,N
     badsx_y(iy)=0.d0
     badsy_y(iy)=0.d0
     do ixp=1,N
        badsx_y(iy)=badsx_y(iy)+sx(ixp,iy) !-- kx=0 part
        badsy_y(iy)=badsy_y(iy)+sy(ixp,iy) !-- kx=0 part
     end do
     badsx_y(iy)=badsx_y(iy)/dble(N)
     badsy_y(iy)=badsy_y(iy)/dble(N)
  end do

  badsx_0=0.d0
  badsy_0=0.d0
  do iy=1,N
     do ix=1,N
        badsx_0=badsx_0+sx(ix,iy)  !-- kx=0, ky=0 part
        badsy_0=badsy_0+sy(ix,iy)  !-- kx=0, ky=0 part
     end do
  end do
  badsx_0=badsx_0/dble(N*N)
  badsy_0=badsy_0/dble(N*N)

  do iy=1,N
     do ix=1,N
        sx_bad(ix,iy)=badsx_x(ix)+badsx_y(iy)-badsx_0  !-- substract kx=0 part, ky=0 part, 
        sx(ix,iy)=sx(ix,iy)-sx_bad(ix,iy)     !-- and add doubly substracted kx=0,ky=0 part
!        sy_bad(ix,iy)=badsy_x(ix)+badsy_y(iy)-badsx_0  !-- substract kx=0 part, ky=0 part, MISTAKE HERE "badsx_0"!!
        sy_bad(ix,iy)=badsy_x(ix)+badsy_y(iy)-badsy_0  !-- substract kx=0 part, ky=0 part, MISTAKE CORRECTED!
        sy(ix,iy)=sy(ix,iy)-sy_bad(ix,iy)     !-- and add doubly substracted kx=0,ky=0 part 
     end do
  end do

END SUBROUTINE rmbadsxsy

!**** "mulink" multiply two 2*N*N k-space complex arrays ******
SUBROUTINE mulink(N,dat_1,dat_2,dat_12)
  
  INTEGER N,i,k
  REAL*8 dat_1(2*N*N),dat_2(2*N*N),dat_12(2*N*N)
  
  do i=1,2*N*N-1,2
     k=i+1   
     dat_12(i)=dat_1(i)*dat_2(i)-dat_1(k)*dat_2(k)
     dat_12(k)=dat_1(i)*dat_2(k)+dat_1(k)*dat_2(i)   
  end do
  
END SUBROUTINE mulink


!*** find e1,e2,e3,dx,dy from sx,sy
!*** NOTE that dx(ix,iy),dy(ix,iy) from this do NOT include uniform part.
SUBROUTINE find_all(sx,sy,e1,e2,e3,p1,p3,dx,dy,e1k0,e3k0,d,th,Edist,A1,A2,A3,B,C1,C3,G1,G2,H1,H2,N)

  INTEGER ix,iy,k,N,NDIM,nn(2),isign
  REAL*8 A1,A2,A3,B,C1,C3,G1,G2,H1,H2
  REAL*8 sx(N,N),sy(N,N),p1(N,N),p3(N,N)
  REAL*8 e1(N,N),e2(N,N),e3(N,N)
  REAL*8 dx(N,N),dy(N,N),e1k0,e3k0,d(N,N),th(N,N),Edist(N,N)
  REAL*8 e1xy(N,N),e2xy(N,N),e3xy(N,N)
  REAL*8 e1x0(N,N),e2x0(N,N),e3x0(N,N)
  REAL*8 e10y(N,N),e20y(N,N),e30y(N,N)
  REAL*8 e100(N,N),e200(N,N),e300(N,N)
  REAL*8 Fe1x_dat(2*N*N),Fe2x_dat(2*N*N),Fe3x_dat(2*N*N)
  REAL*8 Fe1y_dat(2*N*N),Fe2y_dat(2*N*N),Fe3y_dat(2*N*N)
  REAL*8 e1x_dat(2*N*N),e2x_dat(2*N*N),e3x_dat(2*N*N)
  REAL*8 e1y_dat(2*N*N),e2y_dat(2*N*N),e3y_dat(2*N*N)
  REAL*8 e1_dat(2*N*N),e2_dat(2*N*N),e3_dat(2*N*N)
  REAL*8 sx_dat(2*N*N),sy_dat(2*N*N),p1_dat(2*N*N),p3_dat(2*N*N)
  REAL*8 dxt_dat(2*N*N),dyt_dat(2*N*N),dx_dat(2*N*N),dy_dat(2*N*N)
  REAL*8 kx,ky,PI

  PI=4.d0*atan(1.d0)
  NDIM=2
  nn(1)=N
  nn(2)=N

  p1=sx**2+sy**2
  p3=sx**2-sy**2

! Fourier transform sx,sy,p1,p3
  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1
        sx_dat(k)=sx(ix,iy)/dble(N*N)   ! sx_dat in r-space
        sx_dat(k+1)=0.d0
        sy_dat(k)=sy(ix,iy)/dble(N*N)   ! sy_dat in r-space
        sy_dat(k+1)=0.d0
        p1_dat(k)=p1(ix,iy)/dble(N*N)   ! p1_dat in r-space
        p1_dat(k+1)=0.d0
        p3_dat(k)=p3(ix,iy)/dble(N*N)   ! p3_dat in r-space
        p3_dat(k+1)=0.d0
     end do
  end do
  isign=+1
  call fourn(sx_dat,nn,NDIM,isign) ! gives sx_dat in k-space
  call fourn(sy_dat,nn,NDIM,isign) ! gives sy_dat in k-space
  call fourn(p1_dat,nn,NDIM,isign) ! gives p1_dat in k-space
  call fourn(p3_dat,nn,NDIM,isign) ! gives p3_dat in k-space

!=== calculate e1,e2,e3
! calculate e1,e2,e3 for kx.ne.0 and ky.ne.0 from constraint : e1xy,e2xy,e3xy
  do iy=0,N-1
     ky=PI*2.d0*dble(iy)/dble(N)
     do ix=0,N-1
        kx=PI*2.d0*dble(ix)/dble(N)
        
        k=2*(ix+1+N*iy)-1
        
        if(ix.eq.0.or.iy.eq.0) then
           Fe1x_dat(k)  =0.d0
           Fe1x_dat(k+1)=0.d0
           Fe2x_dat(k)  =0.d0
           Fe2x_dat(k+1)=0.d0
           Fe3x_dat(k)  =0.d0  
           Fe3x_dat(k+1)=0.d0
           Fe1y_dat(k)  =0.d0
           Fe1y_dat(k+1)=0.d0
           Fe2y_dat(k)  =0.d0
           Fe2y_dat(k+1)=0.d0
           Fe3y_dat(k)  =0.d0  
           Fe3y_dat(k+1)=0.d0
        else
           ! Note the sign difference from my note for Fe1_dat,Fe2_dat,Fe3_dat.
           ! It's coming from different sign at e^(ikr) in Fourier transf
           Fe1x_dat(k)  =0.d0                                        
           Fe1x_dat(k+1)=cos(ky/2.d0)/sin(ky/2.d0)/sqrt(2.d0)  ! Fe1x in k-space
           Fe2x_dat(k)  =0.d0         
           Fe2x_dat(k+1)=cos(kx/2.d0)/sin(kx/2.d0)/sqrt(2.d0)  ! Fe2x in k-space
           Fe3x_dat(k)  =0.d0                                               
           Fe3x_dat(k+1)=cos(ky/2.d0)/sin(ky/2.d0)/sqrt(2.d0)  ! Fe3x in k-space
           Fe1y_dat(k)  =0.d0                                        
           Fe1y_dat(k+1)=cos(kx/2.d0)/sin(kx/2.d0)/sqrt(2.d0)  ! Fe1y in k-space
           Fe2y_dat(k)  =0.d0         
           Fe2y_dat(k+1)=cos(ky/2.d0)/sin(ky/2.d0)/sqrt(2.d0)  ! Fe2y in k-space
           Fe3y_dat(k)  =0.d0                                               
           Fe3y_dat(k+1)=-cos(kx/2.d0)/sin(kx/2.d0)/sqrt(2.d0) ! Fe3y in k-space 
        end if
     end do
  end do

  call mulink(N,sx_dat,Fe1x_dat,e1x_dat)  ! e1x_dat in k-space
  call mulink(N,sx_dat,Fe2x_dat,e2x_dat)  ! e2x_dat in k-space
  call mulink(N,sx_dat,Fe3x_dat,e3x_dat)  ! e3x_dat in k-space
  call mulink(N,sy_dat,Fe1y_dat,e1y_dat)  ! e1y_dat in k-space
  call mulink(N,sy_dat,Fe2y_dat,e2y_dat)  ! e2y_dat in k-space
  call mulink(N,sy_dat,Fe3y_dat,e3y_dat)  ! e3y_dat in k-space
  
  e1_dat=e1x_dat+e1y_dat
  e2_dat=e2x_dat+e2y_dat
  e3_dat=e3x_dat+e3y_dat

  isign=-1     
  call fourn(e1_dat,nn,NDIM,isign)  ! gives e1_dat in r-space
  call fourn(e2_dat,nn,NDIM,isign)  ! gives e2_dat in r-space
  call fourn(e3_dat,nn,NDIM,isign)  ! gives e3_dat in r-space
  
  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1
        e1xy(ix,iy)=e1_dat(k)   ! gives e1xy
        e2xy(ix,iy)=e2_dat(k)   ! gives e2xy
        e3xy(ix,iy)=e3_dat(k)   ! gives e3xy
     end do
  end do

! calculate e1,e2,e3 for kx.ne.0 and ky.eq.0 : e1x0,e2x0,e3x0
  do iy=0,N-1     !  ky=PI*2.d0*dble(iy)/dble(N)     
     do ix=0,N-1  !  kx=PI*2.d0*dble(ix)/dble(N)
        k=(iy*N+ix)*2+1
        if(ix.ne.0.and.iy.eq.0) then
           e3_dat(k)=-1.d0/(A1+A3)*(C1*p1_dat(k)+C3*p3_dat(k))         ! e3_dat in k-space
           e3_dat(k+1)=-1.d0/(A1+A3)*(C1*p1_dat(k+1)+C3*p3_dat(k+1))    
        else
           e3_dat(k)=0.d0
           e3_dat(k+1)=0.d0
        end if
     end do
  end do
  
  e1_dat= e3_dat          ! e1_dat in k-space
  e2_dat=0.d0             ! e2_dat in k-space

  isign=-1     
  call fourn(e1_dat,nn,NDIM,isign)  ! gives e1_dat in r-space
  call fourn(e2_dat,nn,NDIM,isign)  ! gives e2_dat in r-space
  call fourn(e3_dat,nn,NDIM,isign)  ! gives e3_dat in r-space
  
  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1
        e1x0(ix,iy)=e1_dat(k)   ! gives e1x0
        e2x0(ix,iy)=e2_dat(k)   ! gives e2x0
        e3x0(ix,iy)=e3_dat(k)   ! gives e3x0
     end do
  end do
  
! calculate e1,e2,e3 for kx.eq.0 and ky.ne.0 : e10y,e20y,e30y
  do iy=0,N-1     !  ky=PI*2.d0*dble(iy)/dble(N)     
     do ix=0,N-1  !  kx=PI*2.d0*dble(ix)/dble(N)
        k=(iy*N+ix)*2+1
        if(ix.eq.0.and.iy.ne.0) then
           e3_dat(k)=-1.d0/(A1+A3)*(-C1*p1_dat(k)+C3*p3_dat(k))         ! e3_dat in k-space
           e3_dat(k+1)=-1.d0/(A1+A3)*(-C1*p1_dat(k+1)+C3*p3_dat(k+1))    
        else
           e3_dat(k)=0.d0
           e3_dat(k+1)=0.d0
        end if
     end do
  end do
  
  e1_dat=-e3_dat          ! e1_dat in k-space
  e2_dat=0.d0             ! e2_dat in k-space

  isign=-1     
  call fourn(e1_dat,nn,NDIM,isign)  ! gives e1_dat in r-space
  call fourn(e2_dat,nn,NDIM,isign)  ! gives e2_dat in r-space
  call fourn(e3_dat,nn,NDIM,isign)  ! gives e3_dat in r-space
  
  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1
        e10y(ix,iy)=e1_dat(k)   ! gives e10y
        e20y(ix,iy)=e2_dat(k)   ! gives e20y
        e30y(ix,iy)=e3_dat(k)   ! gives e30y
     end do
  end do

! calculate e1,e2,e3 for kx.eq.0 and ky.eq.0 : e100,e200,e300

  e1_dat=0.d0
  e1_dat(1)=-C1/A1*p1_dat(1)
  e1_dat(2)=-C1/A1*p1_dat(2)     ! e1_dat in k-space
  e1k0=e1_dat(1)           ! e1k0=e1(k=0)

  e2_dat=0.d0             ! e2_dat in k-space

  e3_dat=0.d0
  e3_dat(1)=-C3/A3*p3_dat(1)
  e3_dat(2)=-C3/A3*p3_dat(2)     ! e3_dat in k-space
  e3k0=e3_dat(1)           ! e3k0=e3(k=0)

  isign=-1     
  call fourn(e1_dat,nn,NDIM,isign)  ! gives e1_dat in r-space
  call fourn(e2_dat,nn,NDIM,isign)  ! gives e2_dat in r-space
  call fourn(e3_dat,nn,NDIM,isign)  ! gives e3_dat in r-space
  
  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1
        e100(ix,iy)=e1_dat(k)   ! gives e100
        e200(ix,iy)=e2_dat(k)   ! gives e200
        e300(ix,iy)=e3_dat(k)   ! gives e300
     end do
  end do

! calculate e1,e2,e3 for all k components
  e1=e1xy+e1x0+e10y+e100        ! gives e1
  e2=e2xy+e2x0+e20y+e200        ! gives e2     
  e3=e3xy+e3x0+e30y+e300        ! gives e3

!=== calculate dx,dy,d,th
! Fourier transform e1,e2 
  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1
        e1_dat(k)=e1(ix,iy)/dble(N*N)   ! e1_dat in r-space
        e1_dat(k+1)=0.d0
        e2_dat(k)=e2(ix,iy)/dble(N*N)   ! e2_dat in r-space
        e2_dat(k+1)=0.d0
     end do
  end do
  isign=+1
  call fourn(e1_dat,nn,NDIM,isign) ! gives e1_dat in k-space
  call fourn(e2_dat,nn,NDIM,isign) ! gives e2_dat in k-space

! Find dx and dy
  do iy=0,N-1
     ky=PI*2.d0*dble(iy)/dble(N)
     do ix=0,N-1
        kx=PI*2.d0*dble(ix)/dble(N)
        
        k=2*(ix+1+N*iy)-1

        if(ix.ne.0.and.iy.ne.0) then
           dxt_dat(k)=sx_dat(k)/(-2.d0*sin(kx/2.d0)*sin(ky/2.d0))        ! gives dxt_dat in k-space
           dxt_dat(k+1)=sx_dat(k+1)/(-2.d0*sin(kx/2.d0)*sin(ky/2.d0))
           dyt_dat(k)=sy_dat(k)/(-2.d0*sin(kx/2.d0)*sin(ky/2.d0))        ! gives dyt_dat in k-space
           dyt_dat(k+1)=sy_dat(k+1)/(-2.d0*sin(kx/2.d0)*sin(ky/2.d0))
        else if(ix.ne.0.and.iy.eq.0) then
           dxt_dat(k)=e1_dat(k+1)/(-sqrt(2.d0)*sin(kx/2.d0))  ! note the sign change in i compared to my note
           dxt_dat(k+1)=e1_dat(k)/(sqrt(2.d0)*sin(kx/2.d0))
           dyt_dat(k)=e2_dat(k+1)/(-sqrt(2.d0)*sin(kx/2.d0))
           dyt_dat(k+1)=e2_dat(k)/(sqrt(2.d0)*sin(kx/2.d0))
        else if(ix.eq.0.and.iy.ne.0) then
           dxt_dat(k)=e2_dat(k+1)/(-sqrt(2.d0)*sin(ky/2.d0))  ! note the sign change in i compared to my note
           dxt_dat(k+1)=e2_dat(k)/(sqrt(2.d0)*sin(ky/2.d0))
           dyt_dat(k)=e1_dat(k+1)/(-sqrt(2.d0)*sin(ky/2.d0))
           dyt_dat(k+1)=e1_dat(k)/(sqrt(2.d0)*sin(ky/2.d0))
        else
           dxt_dat(k)=0.d0
           dxt_dat(k+1)=0.d0
           dyt_dat(k)=0.d0
           dyt_dat(k+1)=0.d0
        end if

     end do
  end do

     ! Correct the factor in dxt_dat and dyt_dat due to (1/2,1/2) shift to find dx_dat,dy_dat
     ! Also note the sign difference due to the sign at e^(ikr) in Fourier transf
  do iy=0,N-1
     ky=PI*2.d0*dble(iy)/dble(N)
     do ix=0,N-1
        kx=PI*2.d0*dble(ix)/dble(N)          
        k=2*(ix+1+N*iy)-1 
        dx_dat(k)=dxt_dat(k)*cos((kx+ky)/2.d0)-dxt_dat(k+1)*sin((kx+ky)/2.d0)    ! dx_dat in k-space
        dx_dat(k+1)=dxt_dat(k)*sin((kx+ky)/2.d0)+dxt_dat(k+1)*cos((kx+ky)/2.d0)
        dy_dat(k)=dyt_dat(k)*cos((kx+ky)/2.d0)-dyt_dat(k+1)*sin((kx+ky)/2.d0)    ! dy_dat in k-space
        dy_dat(k+1)=dyt_dat(k)*sin((kx+ky)/2.d0)+dyt_dat(k+1)*cos((kx+ky)/2.d0)
     end do
  end do
  
  isign=-1     
  call fourn(dx_dat,nn,NDIM,isign)  ! gives dx_dat in r-space
  call fourn(dy_dat,nn,NDIM,isign)  ! gives dy_dat in r-space

  do iy=1,N
     do ix=1,N          
        k=2*(ix+N*(iy-1))-1  
        dx(ix,iy)=dx_dat(k)   ! gives dx
        dy(ix,iy)=dy_dat(k)   ! gives dy
        
        d(ix,iy) =sqrt(dx(ix,iy)**2+dy(ix,iy)**2)
        
        if(dx(ix,iy).eq.0.d0.and.dy(ix,iy).eq.0.d0) then
           th(ix,iy)=0.d0
        else if(dx(ix,iy).eq.0.d0.and.dy(ix,iy).gt.0.d0) then
           th(ix,iy)=90.d0
        else if(dx(ix,iy).eq.0.d0.and.dy(ix,iy).lt.0.d0) then
           th(ix,iy)=-90.d0
        else if(dx(ix,iy).gt.0.d0) then
           th(ix,iy)=atan(dy(ix,iy)/dx(ix,iy))*180.d0/PI
        else 
           th(ix,iy)=atan(dy(ix,iy)/dx(ix,iy))*180.d0/PI+180
        end if

     end do
  end do

  Edist=0.5d0*A1*e1**2+0.5d0*A2*e2**2+0.5d0*A3*e3**2+0.5d0*B*(sx**2+sy**2)  &
       + C1*(sx**2+sy**2)*e1 + C3*(sx**2-sy**2)*e3  &
       +0.25d0*G1*(sx**4+sy**4)+0.5d0*G2*sx**2*sy**2    &
       +1.d0/6.d0*H1*(sx**6+sy**6)+1.d0/6.d0*H2*(sx**2+sy**2)*sx**2*sy**2
  
END SUBROUTINE find_all

!**** Write ps file
SUBROUTINE pswrite(fname,N,d,th,dnorm)

  CHARACTER*20 fname
  INTEGER N,ix,iy
  REAL*8 d(N,N),th(N,N),dnorm
  
  open(51,file=fname,status="unknown")
  write(51,*) '%!PS-Adobe-2.0'
  write(51,*) '%%Creator: Ilan Schnell'
  write(51,*) '%%Pages: 1'
  write(51,*) '%%BoundingBox: 20 20 650 650'
  write(51,*) '%%EndComments'
  write(51,*) '%%Page: 1 1'
  write(51,*) ' '
  write(51,*) '/bdf { bind def } def'
  write(51,*) '/Cshow { dup stringwidth pop -2 div 0 rmoveto show } bdf'
  write(51,*) '/ci { newpath 0 360 arc closepath fill } bdf'
  write(51,*) ' '
  write(51,*) '/draw { % draws one arrow'
  write(51,*) '  /angle exch def'
  write(51,*) '  /r exch def'
  write(51,*) '  /y exch def'
  write(51,*) '  /x exch def'
  write(51,*) '  gsave'
  write(51,*) '  x y translate angle rotate'
  write(51,*) '  0 0 .05 ci'
  write(51,*) '  .05 setlinewidth newpath 0 0 moveto r .09 sub 0 lineto stroke'
  write(51,*) '  newpath r .2 sub .19 moveto .2 -.19 rlineto -.2 -.19 rlineto closepath fill'
  write(51,*) '  grestore'
  write(51,*) '} bdf'
  write(51,*) ' '
  write(51,*) '% main program starts here'
  write(51,*) ' '
  write(51,*) '80 80 translate'
  write(51,*) '15 dup scale'
  write(51,*) ' '
  write(51,*) '.05 setlinewidth .4 1 .4 setrgbcolor'
  write(51,*) '1 1 33 { dup'
  write(51,*) '  newpath 1 moveto 0 32 rlineto stroke'
  write(51,*) '  newpath 1 exch moveto 32 0 rlineto stroke'
  write(51,*) '} for'
  write(51,*) '0 0 0 setrgbcolor'
  write(51,*) ' '
  do ix=1,N
     do iy=1,N
        write(51,300) ix, iy, d(ix,iy)*dnorm, th(ix,iy)
300     format(i4,i4,f9.3,f7.1,' draw')
     end do
  end do
  do ix=1,N
     write(51,300) ix, N+1, d(ix,1)*dnorm, th(ix,1)
  end do
  do iy=1,N
     write(51,300) N+1, iy, d(1,iy)*dnorm, th(1,iy)
  end do
  write(51,300)   N+1, N+1, d(1,1)*dnorm, th(1,1)
  write(51,*) ' '
  write(51,*) 'showpage'
  close(51)

END SUBROUTINE pswrite


!***** Fast Fourier Transformation
SUBROUTINE fourn(DATA,nn,ndim,isign)
  INTEGER isign,ndim,nn(ndim)
  REAL*8 DATA(*)
  INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,&
       k2,n,nprev,nrem,ntot
  REAL*8 tempi,tempr
  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
  ntot=1
  DO idim=1,ndim
     ntot=ntot*nn(idim)
  END DO
  nprev=1
  DO idim=1,ndim
     n=nn(idim)
     nrem=ntot/(n*nprev)
     ip1=2*nprev
     ip2=ip1*n
     ip3=ip2*nrem
     i2rev=1
     DO i2=1,ip2,ip1
        IF(i2.LT.i2rev)THEN
           DO i1=i2,i2+ip1-2,2
              DO  i3=i1,ip3,ip2
                 i3rev=i2rev+i3-i2
                 tempr=DATA(i3)
                 tempi=DATA(i3+1)
                 DATA(i3)=DATA(i3rev)
                 DATA(i3+1)=DATA(i3rev+1)
                 DATA(i3rev)=tempr
                 DATA(i3rev+1)=tempi
              END DO
           END DO
        END IF
        ibit=ip2/2
1       IF ((ibit.GE.ip1).AND.(i2rev.GT.ibit)) THEN
           i2rev=i2rev-ibit
           ibit=ibit/2
           GOTO 1
        END IF
        i2rev=i2rev+ibit
     END DO
     ifp1=ip1
2    IF(ifp1.LT.ip2)THEN
        ifp2=2*ifp1
        theta=isign*6.28318530717959d0/(ifp2/ip1)
        wpr=-2.d0*SIN(0.5d0*theta)**2
        wpi=SIN(theta)
        wr=1.d0
        wi=0.d0
        DO  i3=1,ifp1,ip1
           DO  i1=i3,i3+ip1-2,2
              DO  i2=i1,ip3,ifp2
                 k1=i2
                 k2=k1+ifp1
                 tempr=sngl(wr)*DATA(k2)-sngl(wi)*DATA(k2+1)
                 tempi=sngl(wr)*DATA(k2+1)+sngl(wi)*DATA(k2)
                 DATA(k2)=DATA(k1)-tempr
                 DATA(k2+1)=DATA(k1+1)-tempi
                 DATA(k1)=DATA(k1)+tempr
                 DATA(k1+1)=DATA(k1+1)+tempi
              END DO
           END DO
           wtemp=wr
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
        END DO
        ifp1=ifp2
        GOTO 2
     END IF
     nprev=n*nprev
  END DO
END SUBROUTINE fourn
 

FUNCTION RAN1(IDUM)
  REAL*8 RAN1
  DIMENSION R(97)
  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=1./M1)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=1./M2)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF/0/
  IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
     IFF=1
     IX1=MOD(IC1-IDUM,M1)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX2=MOD(IX1,M2)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX3=MOD(IX1,M3)
     DO  J=1,97
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IA2*IX2+IC2,M2)
        R(J)=(DBLE(IX1)+DBLE(IX2)*RM2)*RM1
     END DO
     IDUM=1
  END IF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC3,M3)
  J=1+(97*IX3)/M3
!  IF(J.GT.97.OR.J.LT.1)PAUSE
  RAN1=R(J)
  R(J)=(DBLE(IX1)+DBLE(IX2)*RM2)*RM1
END FUNCTION RAN1


!----------This is as Sharon had written this--------------
SUBROUTINE gauss2(seed,gauu,gNDAT,var)
  IMPLICIT REAL*8 (a-h,o-z) 
  INTEGER gNDAT

  REAL*8 gauu(gNDAT,gNDAT),  var,  gasdev(gNDAT,gNDAT)
  REAL*8 V1(gNDAT,gNDAT),V2(gNDAT,gNDAT),R(gNDAT,gNDAT), FAC
  
  INTEGER seed
  
  DO i=1,gNDAT
     DO j=1,gNDAT
1       V1(i,j)=2.*ran1(seed)-1.
        V2(i,j)=2.*ran1(seed)-1.
        R(i,j)=V1(i,j)**2+V2(i,j)**2
        IF (R(i,j).GE.1.) GO TO 1
        FAC=SQRT(-2.*LOG(R(i,j))/R(i,j))
        GASDEV(i,j)=V2(i,j)*FAC
     END DO
  END DO
  
  gauu = gasdev*SQRT(var)

END SUBROUTINE gauss2

