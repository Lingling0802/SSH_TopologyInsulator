!! To complie the code: f90 -o fle.out  fle.f90
!!                      f90 -fast -xO5 -ftrap=%none -o fle.out fle.f90
!! To run: fle.out
!! Input file: fle.in

PROGRAM Hamiltion

  IMPLICIT NONE

  INTEGER N,NT,m,NEMAX,SITEMAX,EFMAX,My,N1,N2
  
  REAL*8 EPS,PI
  
  PARAMETER(N1=16)
  PARAMETER(N2=16)
  PARAMETER(N=N1*2*N2)
  PARAMETER(NT=N*N) 
  PARAMETER(m=N)
  PARAMETER(My=32)

  PARAMETER(PI=3.1415926)
  PARAMETER(NEMAX=1000)
  PARAMETER(SITEMAX=10)  ! If this is increased, then change corresponding part with "write".
  PARAMETER(EFMAX=5)     ! If this is increased, then change corresponding part with "write".
  
  INTEGER ix,iy,it,itx,ity,ixx,iyy,im,imm,j,i,itt,j0,i0,j00,i00
  INTEGER iwork(N),ierr,modeD,matz
  INTEGER m1,m2,q1,q2,nx,ny
  INTEGER xsA,ysA,xsB,ysB,xoA,yoA,xoB,yoB,xA,yA,xB,yB
  REAL*8 alp,betaR,betaB,betaG,betaY
  REAL*8 H(N,N),ini,w(N),z(N,m),fwork(8*N)


!---- Set alp and bat
   alp=750.0d0
   betaR=250.0d0
   betaB=150.0d0
   betaG=120.0d0
   betaY=120.0d0

!---- Construct the Hamiltion matrix in real space

  ini = 0.d0
  do i=1,N
     do j=1,N
        H(i,j)=ini
     end do
  end do

  do i=1,N
     do j=1,N

        if(i.eq.j) then
          H(i,j)=alp
        end if

        if(MOD(i,2).eq.1) then
            if((i+1).eq.j) then
                H(i,j)=-betaB
                H(j,i)=-betaB
            else if((i+2*N1+1).eq.j) then
                H(i,j)=-betaG
                H(j,i)=-betaG
            end if
        end if

        if(MOD(i,2).eq.0) then
            if((i+1).eq.j) then
               if(MOD(i,2*N1).ne.0) then
                  H(i,j)=-betaR
                  H(j,i)=-betaR
               end if
            end if
        end if

        if(MOD(i,2).eq.0) then
          if((i+2*N1+1).eq.j) then
            if(MOD(i,2*N1).ne.0) then
               H(i,j)=-betaY
               H(j,i)=-betaY
            end if
          end if
        end if

      end do
   end do


!---- Find eigensystems of the matrix
  call rsm(N,N,H,w,m,z,fwork,iwork,ierr)

!---- Write eigenvector magnitude
   do i=1,N
      do j=1,N
         write(10,*) i,j,H(i,j)
      end do
   end do

   do it=1,N
      write(16,*)  w(it)
   end do

   do it=1,N
      write(18,*) sqrt(w(it))
   end do

!------------------------------------------------------------------

END PROGRAM Hamiltion

!========================================================
!           subroutines
!========================================================

subroutine rsm(nm,n,a,w,m,z,fwork,iwork,ierr)
!
integer n,nm,m,iwork(n),ierr
integer k1,k2,k3,k4,k5,k6,k7
double precision a(nm,n),w(n),z(nm,m),fwork(8*n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find all of the eigenvalues and some of the eigenvectors
!     of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        m  the eigenvectors corresponding to the first m eigenvalues
!           are to be computed.
!           if m = 0 then no eigenvectors are computed.
!           if m = n then all of the eigenvectors are computed.
!
!     on output
!
!        w  contains all n eigenvalues in ascending order.
!
!        z  contains the orthonormal eigenvectors associated with
!           the first m eigenvalues.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat,
!           imtqlv and tinvit.  the normal completion code is zero.
!
!        fwork  is a temporary storage array of dimension 8*n.
!
!        iwork  is an integer temporary storage array of dimension n.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 10 * n
if (n .gt. nm .or. m .gt. nm) go to 50
k1 = 1
k2 = k1 + n
k3 = k2 + n
k4 = k3 + n
k5 = k4 + n
k6 = k5 + n
k7 = k6 + n
k8 = k7 + n
if (m .gt. 0) go to 10
!     .......... find eigenvalues only ..........
call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
call  tqlrat(n,w,fwork(k2),ierr)
go to 50
!     .......... find all eigenvalues and m eigenvectors ..........
10 call  tred1(nm,n,a,fwork(k1),fwork(k2),fwork(k3))
call  imtqlv(n,fwork(k1),fwork(k2),fwork(k3),w,iwork,  &
ierr,fwork(k4))
call  tinvit(nm,n,fwork(k1),fwork(k2),fwork(k3),m,w,iwork,z,ierr, &
fwork(k4),fwork(k5),fwork(k6),fwork(k7),fwork(k8))
call  trbak1(nm,n,a,fwork(k2),m,z)
50 return
end


!==================================================================

double precision function epslon (x)
double precision x
!
!     estimate unit roundoff in quantities of size x.
!
double precision a,b,c,eps
!
!     this program should function properly on all systems
!     satisfying the following two assumptions,
!        1.  the base used in representing floating point
!            numbers is not a power of three.
!        2.  the quantity  a  in statement 10 is represented to
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying
!     assumption 2.
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of 1.0 from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/83.
!
a = 4.0d0/3.0d0
10 b = a - 1.0d0
c = b + b + b
eps = dabs(c-1.0d0)
if (eps .eq. 0.0d0) go to 10
epslon = eps*dabs(x)
return
end

subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1)
!
integer i,j,k,l,m,n,ii,mml,tag,ierr
double precision d(n),e(n),e2(n),w(n),rv1(n)
double precision b,c,f,g,p,r,s,tst1,tst2,pythag
integer ind(n)
!
!     this subroutine is a variant of  imtql1  which is a translation of
!     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
!     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues of a symmetric tridiagonal
!     matrix by the implicit ql method and associates with them
!     their corresponding submatrix indices.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2(1) is arbitrary.
!
!     on output
!
!        d and e are unaltered.
!
!        elements of e2, corresponding to elements of e regarded
!          as negligible, have been replaced by zero causing the
!          matrix to split into a direct sum of submatrices.
!          e2(1) is also set to zero.
!
!        w contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        ind contains the submatrix indices associated with the
!          corresponding eigenvalues in w -- 1 for eigenvalues
!          belonging to the first submatrix from the top,
!          2 for those belonging to the second submatrix, etc..
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 0
k = 0
tag = 0
!
do 100 i = 1, n
w(i) = d(i)
if (i .ne. 1) rv1(i-1) = e(i)
100 continue
!
e2(1) = 0.0d0
rv1(n) = 0.0d0
!
do 290 l = 1, n
j = 0
!     .......... look for small sub-diagonal element ..........
105    do 110 m = l, n
if (m .eq. n) go to 120
tst1 = dabs(w(m)) + dabs(w(m+1))
tst2 = tst1 + dabs(rv1(m))
if (tst2 .eq. tst1) go to 120
!     .......... guard against underflowed element of e2 ..........
if (e2(m+1) .eq. 0.0d0) go to 125
110    continue
!
120    if (m .le. k) go to 130
if (m .ne. n) e2(m+1) = 0.0d0
125    k = m
tag = tag + 1
130    p = w(l)
if (m .eq. l) go to 215
if (j .eq. 30) go to 1000
j = j + 1
!     .......... form shift ..........
g = (w(l+1) - p) / (2.0d0 * rv1(l))
r = pythag(g,1.0d0)
g = w(m) - p + rv1(l) / (g + dsign(r,g))
s = 1.0d0
c = 1.0d0
p = 0.0d0
mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
do 200 ii = 1, mml
i = m - ii
f = s * rv1(i)
b = c * rv1(i)
r = pythag(f,g)
rv1(i+1) = r
if (r .eq. 0.0d0) go to 210
s = f / r
c = g / r
g = w(i+1) - p
r = (w(i) - g) * s + 2.0d0 * c * b
p = s * r
w(i+1) = g + p
g = c * r - b
200    continue
!
w(l) = w(l) - p
rv1(l) = g
rv1(m) = 0.0d0
go to 105
!     .......... recover from underflow ..........
210    w(i+1) = w(i+1) - p
rv1(m) = 0.0d0
go to 105
!     .......... order eigenvalues ..........
215    if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
do 230 ii = 2, l
i = l + 2 - ii
if (p .ge. w(i-1)) go to 270
w(i) = w(i-1)
ind(i) = ind(i-1)
230    continue
!
250    i = 1
270    w(i) = p
ind(i) = tag
290 continue
!
go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end

subroutine tqlrat(n,d,e2,ierr)
!
integer i,j,l,m,n,ii,l1,mml,ierr
double precision d(n),e2(n)
double precision b,c,f,g,h,p,r,s,t,epslon,pythag
!
!     this subroutine is a translation of the algol procedure tqlrat,
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 0
if (n .eq. 1) go to 1001
!
do 100 i = 2, n
100 e2(i-1) = e2(i)
!
f = 0.0d0
t = 0.0d0
e2(n) = 0.0d0
!
do 290 l = 1, n
j = 0
h = dabs(d(l)) + dsqrt(e2(l))
if (t .gt. h) go to 105
t = h
b = epslon(t)
c = b * b
!     .......... look for small squared sub-diagonal element ..........
105    do 110 m = l, n
if (e2(m) .le. c) go to 120
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
110    continue
!
120    if (m .eq. l) go to 210
130    if (j .eq. 30) go to 1000
j = j + 1
!     .......... form shift ..........
l1 = l + 1
s = dsqrt(e2(l))
g = d(l)
p = (d(l1) - g) / (2.0d0 * s)
r = pythag(p,1.0d0)
d(l) = s / (p + dsign(r,p))
h = g - d(l)
!
do 140 i = l1, n
140    d(i) = d(i) - h
!
f = f + h
!     .......... rational ql transformation ..........
g = d(m)
if (g .eq. 0.0d0) g = b
h = g
s = 0.0d0
mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
do 200 ii = 1, mml
i = m - ii
p = g * h
r = p + e2(i)
e2(i+1) = s * r
s = e2(i) / r
d(i+1) = h + s * (h + d(i))
g = d(i) - e2(i) / g
if (g .eq. 0.0d0) g = b
h = g * p / r
200    continue
!
e2(l) = s * g
d(l) = h
!     .......... guard against underflow in convergence test ..........
if (h .eq. 0.0d0) go to 210
if (dabs(e2(l)) .le. dabs(c/h)) go to 210
e2(l) = h * e2(l)
if (e2(l) .ne. 0.0d0) go to 130
210    p = d(l) + f
!     .......... order eigenvalues ..........
if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
do 230 ii = 2, l
i = l + 2 - ii
if (p .ge. d(i-1)) go to 270
d(i) = d(i-1)
230    continue
!
250    i = 1
270    d(i) = p
290 continue
!
go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end

double precision function pythag(a,b)
double precision a,b
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
double precision p,r,s,t,u
p = dmax1(dabs(a),dabs(b))
if (p .eq. 0.0d0) go to 20
r = (dmin1(dabs(a),dabs(b))/p)**2
10 continue
t = 4.0d0 + r
if (t .eq. 4.0d0) go to 20
s = r/t
u = 1.0d0 + 2.0d0*s
p = u*p
r = (s/u)**2 * r
go to 10
20 pythag = p
return
end



subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,ierr,rv1,rv2,rv3,rv4,rv6)
!
integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
double precision d(n),e(n),e2(n),w(m),z(nm,m),rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
double precision u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon,pythag
integer ind(m)
!
!     this subroutine is a translation of the inverse iteration tech-
!     nique in the algol procedure tristurm by peters and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
!
!     this subroutine finds those eigenvectors of a tridiagonal
!     symmetric matrix corresponding to specified eigenvalues,
!     using inverse iteration.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        e2 contains the squares of the corresponding elements of e,
!          with zeros corresponding to negligible elements of e.
!          e(i) is considered negligible if it is not larger than
!          the product of the relative machine precision and the sum
!          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
!          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
!          if the eigenvalues are in descending order.  if  bisect,
!          tridib, or  imtqlv  has been used to find the eigenvalues,
!          their output e2 array is exactly what is expected here.
!
!        m is the number of specified eigenvalues.
!
!        w contains the m eigenvalues in ascending or descending order.
!
!        ind contains in its first m positions the submatrix indices
!          associated with the corresponding eigenvalues in w --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!
!     on output
!
!        all input arrays are unaltered.
!
!        z contains the associated set of orthonormal eigenvectors.
!          any vector which fails to converge is set to zero.
!
!        ierr is set to
!          zero       for normal return,
!          -r         if the eigenvector corresponding to the r-th
!                     eigenvalue fails to converge in 5 iterations.
!
!        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 0
if (m .eq. 0) go to 1001
tag = 0
order = 1.0d0 - e2(1)
q = 0
!     .......... establish and process next submatrix ..........
100 p = q + 1
!
do 120 q = p, n
if (q .eq. n) go to 140
if (e2(q+1) .eq. 0.0d0) go to 140
120 continue
!     .......... find vectors by inverse iteration ..........
140 tag = tag + 1
s = 0
!
do 920 r = 1, m
if (ind(r) .ne. tag) go to 920
its = 1
x1 = w(r)
if (s .ne. 0) go to 510
!     .......... check for isolated root ..........
xu = 1.0d0
if (p .ne. q) go to 490
rv6(p) = 1.0d0
go to 870
490    norm = dabs(d(p))
ip = p + 1
!
do 500 i = ip, q
500    norm = dmax1(norm, dabs(d(i))+dabs(e(i)))
!     .......... eps2 is the criterion for grouping,
!                eps3 replaces zero pivots and equal
!                roots are modified by eps3,
!                eps4 is taken very small to avoid overflow ..........
eps2 = 1.0d-3 * norm
eps3 = epslon(norm)
uk = q - p + 1
eps4 = uk * eps3
uk = eps4 / dsqrt(uk)
s = p
505    group = 0
go to 520
!     .......... look for close or coincident roots ..........
510    if (dabs(x1-x0) .ge. eps2) go to 505
group = group + 1
if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
!     .......... elimination with interchanges and
!                initialization of vector ..........
520    v = 0.0d0
!
do 580 i = p, q
rv6(i) = uk
if (i .eq. p) go to 560
if (dabs(e(i)) .lt. dabs(u)) go to 540
!     .......... warning -- a divide check may occur here if
!                e2 array has not been specified correctly ..........
xu = u / e(i)
rv4(i) = xu
rv1(i-1) = e(i)
rv2(i-1) = d(i) - x1
rv3(i-1) = 0.0d0
if (i .ne. q) rv3(i-1) = e(i+1)
u = v - xu * rv2(i-1)
v = -xu * rv3(i-1)
go to 580
540       xu = e(i) / u
rv4(i) = xu
rv1(i-1) = u
rv2(i-1) = v
rv3(i-1) = 0.0d0
560       u = d(i) - x1 - xu * v
if (i .ne. q) v = e(i+1)
580    continue
!
if (u .eq. 0.0d0) u = eps3
rv1(q) = u
rv2(q) = 0.0d0
rv3(q) = 0.0d0
!     .......... back substitution
!                for i=q step -1 until p do -- ..........
600    do 620 ii = p, q
i = p + q - ii
rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
v = u
u = rv6(i)
620    continue
!     .......... orthogonalize with respect to previous
!                members of group ..........
if (group .eq. 0) go to 700
j = r
!
do 680 jj = 1, group
630       j = j - 1
if (ind(j) .ne. tag) go to 630
xu = 0.0d0
!
do 640 i = p, q
640       xu = xu + rv6(i) * z(i,j)
!
do 660 i = p, q
660       rv6(i) = rv6(i) - xu * z(i,j)
!
680    continue
!
700    norm = 0.0d0
!
do 720 i = p, q
720    norm = norm + dabs(rv6(i))
!
if (norm .ge. 1.0d0) go to 840
!     .......... forward substitution ..........
if (its .eq. 5) go to 830
if (norm .ne. 0.0d0) go to 740
rv6(s) = eps4
s = s + 1
if (s .gt. q) s = p
go to 780
740    xu = eps4 / norm
!
do 760 i = p, q
760    rv6(i) = rv6(i) * xu
!     .......... elimination operations on next vector
!                iterate ..........
780    do 820 i = ip, q
u = rv6(i)
!     .......... if rv1(i-1) .eq. e(i), a row interchange
!                was performed earlier in the
!                triangularization process ..........
if (rv1(i-1) .ne. e(i)) go to 800
u = rv6(i-1)
rv6(i-1) = rv6(i)
800       rv6(i) = u - rv4(i) * rv6(i-1)
820    continue
!
its = its + 1
go to 600
!     .......... set error -- non-converged eigenvector ..........
830    ierr = -r
xu = 0.0d0
go to 870
!     .......... normalize so that sum of squares is
!                1 and expand to full order ..........
840    u = 0.0d0
!
do 860 i = p, q
860    u = pythag(u,rv6(i))
!
xu = 1.0d0 / u
!
870    do 880 i = 1, n
880    z(i,r) = 0.0d0
!
do 900 i = p, q
900    z(i,r) = rv6(i) * xu
!
x0 = x1
920 continue
!
if (q .lt. n) go to 100
1001 return
end

subroutine trbak1(nm,n,a,e,m,z)
!
integer i,j,k,l,m,n,nm
double precision a(nm,n),e(n),z(nm,m)
double precision s
!
!     this subroutine is a translation of the algol procedure trbak1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine forms the eigenvectors of a real symmetric
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  tred1.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction by  tred1
!          in its strict lower triangle.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is arbitrary.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     note that trbak1 preserves vector euclidean norms.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
if (m .eq. 0) go to 200
if (n .eq. 1) go to 200
!
do 140 i = 2, n
l = i - 1
if (e(i) .eq. 0.0d0) go to 140
!
do 130 j = 1, m
s = 0.0d0
!
do 110 k = 1, l
110       s = s + a(i,k) * z(k,j)
!     .......... divisor below is negative of h formed in tred1.
!                double division avoids possible underflow ..........
s = (s / a(i,l)) / e(i)
!
do 120 k = 1, l
120       z(k,j) = z(k,j) + s * a(i,k)
!
130    continue
!
140 continue
!
200 return
end

subroutine tred1(nm,n,a,d,e,e2)
!
integer i,j,k,l,n,ii,nm,jp1
double precision a(nm,n),d(n),e(n),e2(n)
double precision f,g,h,scale
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
do 100 i = 1, n
d(i) = a(n,i)
a(n,i) = a(i,i)
100 continue
!     .......... for i=n step -1 until 1 do -- ..........
do 300 ii = 1, n
i = n + 1 - ii
l = i - 1
h = 0.0d0
scale = 0.0d0
if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
do 120 k = 1, l
120    scale = scale + dabs(d(k))
!
if (scale .ne. 0.0d0) go to 140
!
do 125 j = 1, l
d(j) = a(l,j)
a(l,j) = a(i,j)
a(i,j) = 0.0d0
125    continue
!
130    e(i) = 0.0d0
e2(i) = 0.0d0
go to 300
!
140    do 150 k = 1, l
d(k) = d(k) / scale
h = h + d(k) * d(k)
150    continue
!
e2(i) = scale * scale * h
f = d(l)
g = -dsign(dsqrt(h),f)
e(i) = scale * g
h = h - f * g
d(l) = f - g
if (l .eq. 1) go to 285
!     .......... form a*u ..........
do 170 j = 1, l
170    e(j) = 0.0d0
!
do 240 j = 1, l
f = d(j)
g = e(j) + a(j,j) * f
jp1 = j + 1
if (l .lt. jp1) go to 220
!
do 200 k = jp1, l
g = g + a(k,j) * d(k)
e(k) = e(k) + a(k,j) * f
200       continue
!
220       e(j) = g
240    continue
!     .......... form p ..........
f = 0.0d0
!
do 245 j = 1, l
e(j) = e(j) / h
f = f + e(j) * d(j)
245    continue
!
h = f / (h + h)
!     .......... form q ..........
do 250 j = 1, l
250    e(j) = e(j) - h * d(j)
!     .......... form reduced a ..........
do 280 j = 1, l
f = d(j)
g = e(j)
!
do 260 k = j, l
260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
!
280    continue
!
285    do 290 j = 1, l
f = d(j)
d(j) = a(l,j)
a(l,j) = a(i,j)
a(i,j) = f * scale
290    continue
!
300 continue
!
return
end
