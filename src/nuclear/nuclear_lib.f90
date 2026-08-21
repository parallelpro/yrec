!----------------------------------------------------------------------
! nuclear_lib
!----------------------------------------------------------------------
! Aggregates YREC's COMMON-free nuclear-physics utility routines into a
! single module, following the numerics_lib/const_lib convention (see
! GUIDELINES.md, and MESA's own one small public *_lib module per
! domain). Each routine below was originally its own file
! (nuclear/<name>.f90, with per-file conversion notes preserved below);
! bodies are unchanged, only concatenated under one module. Callers use
! this via `use nuclear_lib`. nulosses.f90 contained two subroutines
! (nulosses, azbar); both are included here.
module nuclear_lib
contains

!----------------------------------------------------------------------
! lir
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original lir.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Generic table-lookup interpolation/extrapolation routine, entered
! either as LIR (cubic interpolation/extrapolation, unless num_points
! < 4) or as LIR1 (always linear interpolation/extrapolation).
!
! FOR A SUCH THAT target_z=table_z(A),  SETS result_y(I)=table_y(I,A), I=1,num_y
! table_z(N),table_y(I,N) MUST BE SUPPLIED FOR N=1,num_points AND I=1,num_y
! y_stride IS FIRST DIMENSION OF table_y
! interp_flag IS SET TO 1 FOR INTERPOLATION AND 0 FOR EXTRAPOLATION
! IF continue_search.LE.1, SCAN TO FIND THE table_z(N) WHICH IMMEDIATELY
!   BOUND target_z, STARTING AT N=1
! IF continue_search.GT.1, SCAN STARTS FROM VALUE OF N FROM PREVIOUS
!   CALL OF LIR
! NOTE
! MOST OF THE COMPUTATION IS PERFORMED IN SINGLE PRECISION
subroutine lir(target_z,table_z,result_y,table_y,num_y,y_stride, &
     num_points,continue_search,interp_flag)

      implicit none

      double precision, intent(in) :: target_z
      double precision, intent(in) :: table_z(*)
      double precision, intent(out) :: result_y(*)
      double precision, intent(in) :: table_y(*)
      integer, intent(in) :: num_y, y_stride, num_points, continue_search
      integer, intent(out) :: interp_flag

      double precision :: weight(4)
      integer :: search_idx
      data search_idx/-1/
      save

      integer :: linear_mode
      integer :: stride, stride_m1, y_strided, num_y_strided, table_end
      double precision :: diff
      integer :: pivot, closest
      double precision :: y1, y2, y3, y4, z1, z2, z3, z4, z12, z34
      integer :: y_idx, j, k, base_idx
      double precision :: yy

      linear_mode=0
      go to 1
      entry lir1(target_z,table_z,result_y,table_y,num_y,y_stride, &
           num_points,continue_search,interp_flag)
      linear_mode=1
    1 continue
      stride=1
! CHECK NT AND RESET IL IF NECESSARY
      if(num_points.lt.2) go to 101
      if(num_points.lt.4) linear_mode=1
! ADDRESSING CONSTANTS
      interp_flag=1
      stride_m1=stride-1
      y_strided=stride*y_stride
      num_y_strided=(num_y-1)*stride+1
      table_end=(num_points-1)*stride+1
      diff=table_z(table_end)-table_z(1)
! SET INDEX FOR START OF SEARCH
      search_idx=(search_idx-2)*stride+1
      if(continue_search.le.1.or.search_idx.lt.1) search_idx=1
! DETERMINE POSITION OF target_z WITHIN table_z
    2 if(search_idx.gt.table_end) go to 8
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF(DIFF) 4,102,3
!     3 IF(ZI(N)-Z) 5,6,9
!     4 IF(ZI(N)-Z) 9,6,5
      if (diff .lt. 0.0) then
         goto 4
      else if (diff .eq. 0.0) then
         goto 102
      else
         goto 3
      end if
    3 if (table_z(search_idx) .lt. target_z) then
         goto 5
      else if (table_z(search_idx) .eq. target_z) then
         goto 6
      else
         goto 9
      end if
    4 if (table_z(search_idx) .lt. target_z) then
         goto 9
      else if (table_z(search_idx) .eq. target_z) then
         goto 6
      else
         goto 5
      end if
    5 search_idx=search_idx+stride
      go to 2
! SET Y WHEN Z LIES ON A MESH POINT
    6 base_idx=(search_idx-1)*y_stride
      do 7 y_idx=1,num_y_strided
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       Y(I)=YI(I+J)
!     7 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         result_y(y_idx)=table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
    7 continue
      go to 30
! CONTROL WHEN Z DOES NOT LIE ON A MESH POINT
    8 interp_flag=0
    9 if(search_idx.le.1) interp_flag=0
      if(linear_mode.eq.1) go to 20
! CUBIC INTERPOLATION/EXTRAPOLATION
! PIVOTAL POINT (M) AND POINT (K) CLOSEST TO Z
!    10 M=N
      pivot=search_idx
      closest=3
      if(search_idx.gt.1+stride) go to 11
      pivot=1+stride+stride
      closest=search_idx
   11 if(search_idx.lt.table_end) go to 12
      pivot=table_end-stride
      closest=4
! WEIGHTING FACTORS
   12 y1=table_z(pivot-stride*2)
      y2=table_z(pivot-stride)
      y3=table_z(pivot)
      y4=table_z(pivot+stride)
      z1=target_z-y1
      z2=target_z-y2
      z3=target_z-y3
      z4=target_z-y4
!    13 Z12=Z1*Z2
      z12=z1*z2
      z34=z3*z4
!    14 A(1)=Z2*Z34/((Y1-Y2)*(Y1-Y3)*(Y1-Y4))
      weight(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      weight(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      weight(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      weight(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
! CORRECT A(K)
!    15 DIFF=A(1)+A(2)+A(3)+A(4)
      diff=weight(1)+weight(2)+weight(3)+weight(4)
      weight(closest)=(1.d0+weight(closest))-diff
! COMPUTE Y
!    16 M=(M-1)/IR-3
      pivot=(pivot-1)/stride-3
      pivot=pivot*y_strided
      do 18 y_idx=1,num_y_strided
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       K=I+M
!       YY=0.D0
!       DO 17 J=1,4
!       K=K+IRD
!       DIFF=YI(K)
!    17 YY=YY+A(J)*DIFF
!       Y(I)=YY
!    18 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         k=y_idx+pivot
         yy=0.d0
         do 17 j=1,4
            k=k+y_strided
            diff=table_y(k)
            yy=yy+weight(j)*diff
   17    continue
         result_y(y_idx)=yy
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
   18 continue
      go to 30
! LINEAR INTERPOLATION/EXTRAPOLATION
   20 if(search_idx.eq.1) search_idx=1+stride
      if(search_idx.gt.table_end) search_idx=table_end
      z1=table_z(search_idx)
      y1=(z1-target_z)/(z1-table_z(search_idx-stride))
      y2=1.0d0-y1
      base_idx=(search_idx-1)*y_stride
      pivot=base_idx-y_strided
      do 21 y_idx=1,num_y_strided,stride
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
!       Y(I)=Y1*YI(I+M)+Y2*YI(I+J)
!    21 IF(Y(I).EQ.0.D0) Y(I+IR1)=0.D0
         result_y(y_idx)=y1*table_y(y_idx+pivot)+y2*table_y(y_idx+base_idx)
         if(result_y(y_idx).eq.0.d0) result_y(y_idx+stride_m1)=0.d0
   21 continue
! RESET N
   30 search_idx=(search_idx+stride-1)/stride
      return
! DIAGNOSTICS
  101 continue
      return
  102 continue
      return
!  1001 FORMAT(/1X,10('*'),5X,'THERE ARE FEWER THAN TWO DATA POINTS IN',
!      *      ' LIR     NT =',I4,5X,10('*')/)
!  1002 FORMAT(/1X,10('*'),5X,'EXTREME VALUES OF INDEPENDENT VARIABLE',
!      *      ' EQUAL IN LIR',5X,10('*')/16X,'ZI(1) =',1PE13.5,',   ',
!      *       'ZI(',I4,') =',1PE13.5/)
end subroutine lir

!----------------------------------------------------------------------
! neutrino
!----------------------------------------------------------------------
! Modernized (free-form) 2026 as part of the YREC readability refactor.
! Logic and numerics are unchanged from the original neutrino.f; only
! the source form, explicit typing, and header/comment style were
! updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Added by Grant Newsham 9/06.
!
! Sets up the (currently 4-species: H, He, C, O) composition arrays
! expected by nulosses/sneut and calls through to them. Returns the
! neutrino energy-loss rate (snu) and its derivatives with respect to
! temperature and density (dsnudt, dsnudd).
subroutine neutrino(temp,den,x,y,zc,zo,snu,dsnudt,dsnudd)
      implicit none

      double precision, intent(in) :: temp, den, x, y, zc, zo
      double precision, intent(out) :: snu, dsnudt, dsnudd

      integer, parameter :: ionmax = 4
      double precision :: xmass(ionmax), ymass(ionmax), aion(ionmax), &
           zion(ionmax)

!..set the mass fractions, z's and a's of the composition
!..hydrogen

      aion(1)  = 1.0d0
      zion(1)  = 1.0d0

!..helium

      aion(2)  = 4.0d0
      zion(2)  = 2.0d0

!..carbon 12

      aion(3)  = 12.0d0
      zion(3)  = 6.0d0

!..oxygen 16

      aion(4)  = 16.0d0
      zion(4)  = 8.0d0

      xmass(1) = x
      xmass(2) = y
      xmass(3) = zc
      xmass(4) = zo

      call nulosses(temp,den,snu,xmass,ymass,aion,zion,dsnudt,dsnudd)

      return
end subroutine neutrino

!----------------------------------------------------------------------
! nulosses
!----------------------------------------------------------------------
! Modernized (free-form) 2026 as part of the YREC readability refactor.
! Logic and numerics are unchanged from the original nulosses.f; only
! the source form and header/comment style were updated -- this file's
! variable names were already lowercase and descriptive and were left
! as-is. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Neutrino loss subroutine modified by G. Newsham from codes of
! N. Itoh & F. Timmes (9/06).
!
! Currently uses just the carbon coefficients for the bremsstrahlung
! neutrino losses (the differences for other elements is slight.)
! We also currently set the fermi temperature cuttoff at 0.01 Tfermi.
! This is applicable for helium but will need to be amended for
! carbon and oxygen white dwarfs and cores of AGB stars.

subroutine nulosses(temp,den,snu,xmass,ymass,aion,zion, &
                     dsnudt,dsnudd)

      implicit none
      save

!..tests the neutrino loss rate routine

!..ionmax  = number of isotopes in the network
!..xmass   = mass fractions
!..ymass   = molar fractions
!..aion    = number of nucleons
!..zion    = number of protons

      integer          ionmax
      parameter        (ionmax=4)
      double precision xmass(ionmax),ymass(ionmax), &
                       aion(ionmax),zion(ionmax),temp,den,abar,zbar, &
                       snu,dsnudt,dsnudd,dsnuda,dsnudz



!..get abar and zbar
      call azbar(xmass,aion,zion,ionmax, &
                 ymass,abar,zbar)



!..get the neutrino losses
      call sneut(temp,den,abar,zbar, &
                  snu,dsnudt,dsnudd,dsnuda,dsnudz)



      return
end subroutine nulosses


subroutine azbar(xmass,aion,zion,ionmax, &
                 ymass,abar,zbar)
      implicit none
      save

!..input
!..mass fractions     = xmass(1:ionmax)
!..number of nucleons = aion(1:ionmax)
!..charge of nucleus  = zion(1:ionmax)
!..number of isotopes = ionmax

!..output:
!..molar abundances        = ymass(1:ionmax),
!..mean number of nucleons = abar
!..mean nucleon charge     = zbar



!..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       ymass(ionmax),abar,zbar,zbarxx,ytot1


      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
end subroutine azbar

!----------------------------------------------------------------------
! ratext
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ratext.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Diagonal rational-function extrapolator used to extrapolate a
! sequence of sub-stepped burning-rate estimates (indexed by
! decreasing step size) to the zero-step-size limit. Same algorithm as
! subroutine RZEXTR in Numerical Recipes, p.566. Called repeatedly by
! liburn/liburn2 with increasing est_index as the lithium/beryllium
! burning sub-stepping is refined.
subroutine ratext(est_index, x_est, y_est, y_extrap, y_err, num_vars, &
     max_use)
      implicit none
      integer, parameter :: imax=11, nmax=15, ncol=7

      integer, intent(in) :: est_index, num_vars, max_use
      double precision, intent(in) :: x_est
      double precision, intent(in) :: y_est(num_vars)
      double precision, intent(out) :: y_extrap(num_vars), y_err(num_vars)

      double precision :: x_hist(imax), tableau(nmax,ncol), fx(ncol)
      save

! SAME AS SR RZEXTR FROM NUMERICAL RECIPES, P.566.
!
! SAVE CURRENT INDEPENDENT VARIABLE.
      integer :: var_idx, k, k2, num_use
      double precision :: yy, v, c, b1, b, ddy

      x_hist(est_index) = x_est
      if (est_index.eq.1) then
         do var_idx = 1, num_vars
            y_extrap(var_idx) = y_est(var_idx)
            tableau(var_idx,1) = y_est(var_idx)
            y_err(var_idx) = y_est(var_idx)
         end do
      else
!        USE AT MOST max_use PREVIOUS MEMBERS.
         num_use = min(est_index,max_use)
         do k = 1, num_use-1
            fx(k+1) = x_hist(est_index-k)/x_est
         end do
!        EVALUATE NEXT DIAGONAL IN TABLEAU.
         do var_idx=1,num_vars
            yy = y_est(var_idx)
            v = tableau(var_idx,1)
            c = yy
            tableau(var_idx,1) = yy
            do k2 = 2, num_use
               b1 = fx(k2)*v
               b = b1 - c
!              CARE NEEDED TO AVOID DIVISION BY ZERO.
               if (b.ne.0d0) then
                  b = (c - v)/b
                  ddy = c*b
                  c = b1*b
               else
                  ddy = v
               end if
               v = tableau(var_idx,k2)
               tableau(var_idx,k2) = ddy
               yy = yy + ddy
            end do
            y_err(var_idx) = ddy
            y_extrap(var_idx) = yy
         end do
      end if
      return
end subroutine ratext

!----------------------------------------------------------------------
! sneut
!----------------------------------------------------------------------
! Modernized (free-form) 2026 as part of the YREC readability refactor.
! Logic and numerics are unchanged from the original sneut.f; only the
! source form and header/comment style were updated -- this file's
! variable names were already lowercase and physically descriptive
! (matching the notation of Itoh et al. 1996, ApJS 102, 411) and were
! left as-is. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Computes neutrino energy-loss rates (erg/g/s) and their derivatives
! with respect to temperature, density, mean atomic weight, and mean
! charge, from the analytic fits of Itoh et al. (1996) to the pair,
! plasma, photoneutrino, bremsstrahlung, and recombination processes.
subroutine sneut(temp,den,abar,zbar, &
                  snu,dsnudt,dsnudd,dsnuda,dsnudz)
implicit none
save


!..this routine computes neutrino losses from the analytic fits of
!..itoh et al. apjs 102, 411, 1996, and also returns their derivatives.


!..input:
!..temp = temperature
!..den  = density
!..abar = mean atomic weight
!..zbar = mean charge


!..output:
!..snu    = total neutrino loss rate in erg/g/sec
!..dsnudt = derivative of snu with temperature
!..dsnudd = derivative of snu with density
!..dsnuda = derivative of snu with abar
!..dsnudz = derivative of snu with zbar



!..declare the pass
double precision temp,den,abar,zbar, &
                 snu,dsnudt,dsnudd,dsnuda,dsnudz


!..local variables
!       integer          i
double precision spair,spairdt,spairdd,spairda,spairdz, &
                 splas,splasdt,splasdd,splasda,splasdz, &
                 sphot,sphotdt,sphotdd,sphotda,sphotdz, &
                 sbrem,sbremdt,sbremdd,sbremda,sbremdz, &
                 sreco,srecodt,srecodd,srecoda,srecodz


double precision t9,xl,xldt,xlp5,xl2,xl3,xl4,xl5,xl6,xl7,xl8,xl9, &
                 xlmp5,xlm1,xlm2,xlm3,xlm4,xlnt,cc,den6,tfermi, &
                 a0,a1,a2,a3,b1,b2,c00,c01,c02,c03,c04,c05,c06, &  ! b3
                 c10,c11,c12,c13,c14,c15,c16,c20,c21,c22,c23,c24, &
                 c25,c26,dd00,dd01,dd02,dd03,dd04,dd05,dd11,dd12, &
                 dd13,dd14,dd15,dd21,dd22,dd23,dd24,dd25,b,c,d,f0, &
                 f1,deni,tempi,abari,zbari,f2,f3,z,xmue,ye, &
                 dum,dumdt,dumdd,dumda,dumdz, &  ! rp1,rn1
                 gum,gumdt,gumdd,gumda,gumdz



!..pair production
double precision rm,rmdd,rmda,rmdz,rmi,gl,gldt, &
                 zeta,zetadt,zetadd,zetada,zetadz,zeta2,zeta3, &
                 xnum,xnumdt,xnumdd,xnumda,xnumdz, &
                 xden,xdendt,xdendd,xdenda,xdendz, &
                 fpair,fpairdt,fpairdd,fpairda,fpairdz, &
                 qpair,qpairdt,qpairdd,qpairda,qpairdz


!..plasma
double precision gl2,gl2dt,gl2dd,gl2da,gl2dz,gl12,gl32,gl72,gl6, &
                 ft,ftdt,ftdd,ftda,ftdz,fl,fldt,fldd,flda,fldz, &
                 fxy,fxydt,fxydd,fxyda,fxydz


!..photo
double precision tau,taudt,cos1,cos2,cos3,cos4,cos5,sin1,sin2, &
                 sin3,sin4,sin5,last,xast, &
                 fphot,fphotdt,fphotdd,fphotda,fphotdz, &
                 qphot,qphotdt,qphotdd,qphotda,qphotdz


!..brem
double precision t8,t812,t832,t82,t83,t85,t86,t8m1,t8m2,t8m3,t8m5, &
                 t8m6, &
                 eta,etadt,etadd,etada,etadz,etam1,etam2,etam3, &
                 fbrem,fbremdt,fbremdd,fbremda,fbremdz, &
                 gbrem,gbremdt,gbremdd,gbremda,gbremdz, &
                 u,gm1,gm2,gm13,gm23,gm43,gm53,v,w,fb,gt,gb, &
                 fliq,fliqdt,fliqdd,fliqda,fliqdz, &
                 gliq,gliqdt,gliqdd,gliqda,gliqdz


!..recomb
double precision ifermi12,zfermim12,nu,nudt,nudd,nuda,nudz, &
                 nu2,nu3,bigj,bigjdt,bigjdd,bigjda,bigjdz




!..numerical constants
double precision pi,fac1,fac2,fac3,oneth,twoth,con1,sixth,iln10
parameter        (pi     = 3.1415926535897932384d0, &
                  fac1   = 5.0d0 * pi / 3.0d0, &
                  fac2   = 10.0d0 * pi, &
                  fac3   = pi / 5.0d0, &
                  oneth  = 1.0d0/3.0d0, &
                  twoth  = 2.0d0/3.0d0, &
                  con1   = 1.0d0/5.9302d0, &
                  sixth  = 1.0d0/6.0d0, &
                  iln10  = 4.342944819032518d-1)



!..theta is sin**2(theta_weinberg) = 0.2319 plus/minus 0.00005 (1996)
!..xnufam is the number of neutrino flavors = 3.02 plus/minus 0.005 (1998)
!..change theta and xnufam if need be, and the changes will automatically
!..propagate through the routine. cv and ca are the vector and axial currents.


double precision theta,xnufam,cv,ca,cvp,cap,tfac1,tfac2,tfac3, &
                 tfac4,tfac5,tfac6
parameter        (theta  = 0.2319d0, &
                  xnufam = 3.0d0, &
                  cv     = 0.5d0 + 2.0d0 * theta, &
                  cvp    = 1.0d0 - cv, &
                  ca     = 0.5d0, &
                  cap    = 1.0d0 - ca, &
                  tfac1  = cv*cv + ca*ca + &
                           (xnufam-1.0d0) * (cvp*cvp+cap*cap), &
                  tfac2  = cv*cv - ca*ca + &
                           (xnufam-1.0d0) * (cvp*cvp - cap-cap), &
                  tfac3  = tfac2/tfac1, &
                  tfac4  = 0.5d0 * tfac1, &
                  tfac5  = 0.5d0 * tfac2, &
                  tfac6  = cv*cv + 1.5d0*ca*ca + (xnufam - 1.0d0)* &
                           (cvp*cvp + 1.5d0*cap*cap))




!..initialize

spair   = 0.0d0
spairdt = 0.0d0
spairdd = 0.0d0
spairda = 0.0d0
spairdz = 0.0d0


splas   = 0.0d0
splasdt = 0.0d0
splasdd = 0.0d0
splasda = 0.0d0
splasdz = 0.0d0


sphot   = 0.0d0
sphotdt = 0.0d0
sphotdd = 0.0d0
sphotda = 0.0d0
sphotdz = 0.0d0


sbrem   = 0.0d0
sbremdt = 0.0d0
sbremdd = 0.0d0
sbremda = 0.0d0
sbremdz = 0.0d0


sreco   = 0.0d0
srecodt = 0.0d0
srecodd = 0.0d0
srecoda = 0.0d0
srecodz = 0.0d0


snu     = 0.0d0
dsnudt  = 0.0d0
dsnudd  = 0.0d0
dsnuda  = 0.0d0
dsnudz  = 0.0d0


if (temp .lt. 1.0e7) return



!..to avoid lots of divisions
deni  = 1.0d0/den
tempi = 1.0d0/temp
abari = 1.0d0/abar
zbari = 1.0d0/zbar



!..some composition variables
ye    = zbar*abari
xmue  = abar*zbari





!..some frequent factors
t9     = temp * 1.0d-9
xl     = t9 * con1
xldt   = 1.0d-9 * con1
xlp5   = sqrt(xl)
xl2    = xl*xl
xl3    = xl2*xl
xl4    = xl3*xl
xl5    = xl4*xl
xl6    = xl5*xl
xl7    = xl6*xl
xl8    = xl7*xl
xl9    = xl8*xl
xlmp5  = 1.0d0/xlp5
xlm1   = 1.0d0/xl
xlm2   = xlm1*xlm1
xlm3   = xlm1*xlm2
xlm4   = xlm1*xlm3


rm     = den*ye
rmdd   = ye
rmda   = -rm*abari
rmdz   = den*abari
rmi    = 1.0d0/rm


a0     = rm * 1.0d-9
a1     = a0**oneth
zeta   = a1 * xlm1
zetadt = -a1 * xlm2 * xldt
a2     = oneth * a1*rmi * xlm1
zetadd = a2 * rmdd
zetada = a2 * rmda
zetadz = a2 * rmdz

zeta2 = zeta * zeta
zeta3 = zeta2 * zeta





!..pair neutrino section
!..for reactions like e+ + e- => nu_e + nubar_e


!..equation 2.8
gl   = 1.0d0 - 13.04d0*xl2 +133.5d0*xl4 +1534.0d0*xl6 +918.6d0*xl8
gldt = xldt*(-26.08d0*xl +534.0d0*xl3 +9204.0d0*xl5 +7348.8d0*xl7)


!..equation 2.7


a1     = 6.002d19 + 2.084d20*zeta + 1.872d21*zeta2
a2     = 2.084d20 + 2.0d0*1.872d21*zeta


if (t9 .lt. 10.0) then
 b1     = exp(-5.5924d0*zeta)
 b2     = -b1*5.5924d0
else
 b1     = exp(-4.9924d0*zeta)
 b2     = -b1*4.9924d0
end if

xnum   = a1 * b1
c      = a2*b1 + a1*b2
xnumdt = c*zetadt
xnumdd = c*zetadd
xnumda = c*zetada
xnumdz = c*zetadz


if (t9 .lt. 10.0) then
 a1   = 9.383d-1*xlm1 - 4.141d-1*xlm2 + 5.829d-2*xlm3
 a2   = -9.383d-1*xlm2 + 2.0d0*4.141d-1*xlm3 - 3.0d0*5.829d-2*xlm4
else
 a1   = 1.2383d0*xlm1 - 8.141d-1*xlm2
 a2   = -1.2383d0*xlm2 + 2.0d0*8.141d-1*xlm3
end if


b1   = 3.0d0*zeta2


xden   = zeta3 + a1
xdendt = b1*zetadt + a2*xldt
xdendd = b1*zetadd
xdenda = b1*zetada
xdendz = b1*zetadz


a1      = 1.0d0/xden
fpair   = xnum*a1
fpairdt = (xnumdt - fpair*xdendt)*a1
fpairdd = (xnumdd - fpair*xdendd)*a1
fpairda = (xnumda - fpair*xdenda)*a1
fpairdz = (xnumdz - fpair*xdendz)*a1



!..equation 2.6
a1     = 10.7480d0*xl2 + 0.3967d0*xlp5 + 1.005d0
a2     = xldt*(2.0d0*10.7480d0*xl + 0.5d0*0.3967d0*xlmp5)
xnum   = 1.0d0/a1
xnumdt = -xnum*xnum*a2


a1     = 7.692d7*xl3 + 9.715d6*xlp5
a2     = xldt*(3.0d0*7.692d7*xl2 + 0.5d0*9.715d6*xlmp5)


c      = 1.0d0/a1
b1     = 1.0d0 + rm*c


xden   = b1**(-0.3d0)


d      = -0.3d0*xden/b1
xdendt = -d*rm*c*c*a2
xdendd = d*rmdd*c
xdenda = d*rmda*c
xdendz = d*rmdz*c


qpair   = xnum*xden
qpairdt = xnumdt*xden + xnum*xdendt
qpairdd = xnum*xdendd
qpairda = xnum*xdenda
qpairdz = xnum*xdendz




!..equation 2.5
a1    = exp(-2.0d0*xlm1)
a2    = a1*2.0d0*xlm2*xldt


spair   = a1*fpair
spairdt = a2*fpair + a1*fpairdt
spairdd = a1*fpairdd
spairda = a1*fpairda
spairdz = a1*fpairdz


a1      = spair
spair   = gl*a1
spairdt = gl*spairdt + gldt*a1
spairdd = gl*spairdd
spairda = gl*spairda
spairdz = gl*spairdz


a1      = tfac4*(1.0d0 + tfac3 * qpair)
a2      = tfac4*tfac3


a3      = spair
spair   = a1*a3
spairdt = a1*spairdt + a2*qpairdt*a3
spairdd = a1*spairdd + a2*qpairdd*a3
spairda = a1*spairda + a2*qpairda*a3
spairdz = a1*spairdz + a2*qpairdz*a3





!..plasma neutrino section
!..for collective reactions like gamma_plasmon => nu_e + nubar_e
!..equation 4.6


a1   = 1.019d-6*rm
a2   = a1**twoth
a3   = twoth*a2/a1


b1   =  sqrt(1.0d0 + a2)
b2   = 1.0d0/b1

c00  = 1.0d0/(temp*temp*b1)


gl2   = 1.1095d11 * rm * c00


gl2dt = -2.0d0*gl2*tempi
d     = rm*c00*b2*0.5d0*b2*a3*1.019d-6
gl2dd = 1.1095d11 * (rmdd*c00  - d*rmdd)
gl2da = 1.1095d11 * (rmda*c00  - d*rmda)
gl2dz = 1.1095d11 * (rmdz*c00  - d*rmdz)



gl    = sqrt(gl2)
gl12  = sqrt(gl)
gl32  = gl * gl12
gl72  = gl2 * gl32
gl6   = gl2 * gl2 * gl2



!..equation 4.7
ft   = 2.4d0 + 0.6d0*gl12 + 0.51d0*gl + 1.25d0*gl32
gum  = 1.0d0/gl2
a1   =(0.25d0*0.6d0*gl12 +0.5d0*0.51d0*gl +0.75d0*1.25d0*gl32)*gum
ftdt = a1*gl2dt
ftdd = a1*gl2dd
ftda = a1*gl2da
ftdz = a1*gl2dz



!..equation 4.8
a1   = 8.6d0*gl2 + 1.35d0*gl72
a2   = 8.6d0 + 1.75d0*1.35d0*gl72*gum


b1   = 225.0d0 - 17.0d0*gl + gl2
b2   = -0.5d0*17.0d0*gl*gum + 1.0d0


c    = 1.0d0/b1
fl   = a1*c


d    = (a2 - fl*b2)*c
fldt = d*gl2dt
fldd = d*gl2dd
flda = d*gl2da
fldz = d*gl2dz



!..equation 4.9 and 4.10
cc   = log10(2.0d0*rm)
xlnt = log10(temp)


xnum   = sixth * (17.5d0 + cc - 3.0d0*xlnt)
xnumdt = -iln10*0.5d0*tempi
a2     = iln10*sixth*rmi
xnumdd = a2*rmdd
xnumda = a2*rmda
xnumdz = a2*rmdz


xden   = sixth * (-24.5d0 + cc + 3.0d0*xlnt)
xdendt = iln10*0.5d0*tempi
xdendd = a2*rmdd
xdenda = a2*rmda
xdendz = a2*rmdz



!..equation 4.11
if (abs(xnum) .gt. 0.7d0  .or.  xden .lt. 0.0d0) then
 fxy   = 1.0d0
 fxydt = 0.0d0
 fxydd = 0.0d0
 fxydz = 0.0d0
 fxyda = 0.0d0


else


 a1  = 0.39d0 - 1.25d0*xnum - 0.35d0*sin(4.5d0*xnum)
 a2  = -1.25d0 - 4.5d0*0.35d0*cos(4.5d0*xnum)


 b1  = 0.3d0 * exp(-1.0d0*(4.5d0*xnum + 0.9d0)**2)
 b2  = -b1*2.0d0*(4.5d0*xnum + 0.9d0)*4.5d0


 c   = min(0.0d0, xden - 1.6d0 + 1.25d0*xnum)
 if (c .eq. 0.0) then
  dumdt = 0.0d0
  dumdd = 0.0d0
  dumda = 0.0d0
  dumdz = 0.0d0
 else
  dumdt = xdendt + 1.25d0*xnumdt
  dumdd = xdendd + 1.25d0*xnumdd
  dumda = xdenda + 1.25d0*xnumda
  dumdz = xdendz + 1.25d0*xnumdz
 end if


 d   = 0.57d0 - 0.25d0*xnum
 a3  = c/d
 c00 = exp(-1.0d0*a3**2)


 f1  = -c00*2.0d0*a3/d
 c01 = f1*(dumdt + a3*0.25d0*xnumdt)
 c02 = f1*(dumdd + a3*0.25d0*xnumdd)
 c03 = f1*(dumda + a3*0.25d0*xnumda)
 c04 = f1*(dumdz + a3*0.25d0*xnumdz)


 fxy   = 1.05d0 + (a1 - b1)*c00
 fxydt = (a2*xnumdt -  b2*xnumdt)*c00 + (a1-b1)*c01
 fxydd = (a2*xnumdd -  b2*xnumdd)*c00 + (a1-b1)*c02
 fxyda = (a2*xnumda -  b2*xnumda)*c00 + (a1-b1)*c03
 fxydz = (a2*xnumdz -  b2*xnumdz)*c00 + (a1-b1)*c04


end if




!..equation 4.1 and 4.5
splas   = (ft + fl) * fxy
splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz


a2      = exp(-gl)
a3      = -0.5d0*a2*gl*gum


a1      = splas
splas   = a2*a1
splasdt = a2*splasdt + a3*gl2dt*a1
splasdd = a2*splasdd + a3*gl2dd*a1
splasda = a2*splasda + a3*gl2da*a1
splasdz = a2*splasdz + a3*gl2dz*a1


a2      = gl6
a3      = 3.0d0*gl6*gum


a1      = splas
splas   = a2*a1
splasdt = a2*splasdt + a3*gl2dt*a1
splasdd = a2*splasdd + a3*gl2dd*a1
splasda = a2*splasda + a3*gl2da*a1
splasdz = a2*splasdz + a3*gl2dz*a1



a2      = 0.93153d0 * 3.0d21 * xl9
a3      = 0.93153d0 * 3.0d21 * 9.0d0*xl8*xldt


a1      = splas
splas   = a2*a1
splasdt = a2*splasdt + a3*a1
splasdd = a2*splasdd
splasda = a2*splasda
splasdz = a2*splasdz





!..photoneutrino process section
!..for reactions like e- + gamma => e- + nu_e + nubar_e
!..                   e+ + gamma => e+ + nu_e + nubar_e
!..equation 3.8 for tau, equation 3.6 for cc,
!..and table 2 written out for speed
if (temp .ge. 1.0d7  .and. temp .lt. 1.0d8) then
 tau  =  log10(temp * 1.0d-7)
 cc   =  0.5654d0 + tau
 c00  =  1.008d11
 c01  =  0.0d0
 c02  =  0.0d0
 c03  =  0.0d0
 c04  =  0.0d0
 c05  =  0.0d0
 c06  =  0.0d0
 c10  =  8.156d10
 c11  =  9.728d8
 c12  = -3.806d9
 c13  = -4.384d9
 c14  = -5.774d9
 c15  = -5.249d9
 c16  = -5.153d9
 c20  =  1.067d11
 c21  = -9.782d9
 c22  = -7.193d9
 c23  = -6.936d9
 c24  = -6.893d9
 c25  = -7.041d9
 c26  = -7.193d9
 dd01 =  0.0d0
 dd02 =  0.0d0
 dd03 =  0.0d0
 dd04 =  0.0d0
 dd05 =  0.0d0
 dd11 = -1.879d10
 dd12 = -9.667d9
 dd13 = -5.602d9
 dd14 = -3.370d9
 dd15 = -1.825d9
 dd21 = -2.919d10
 dd22 = -1.185d10
 dd23 = -7.270d9
 dd24 = -4.222d9
 dd25 = -1.560d9


else if (temp .ge. 1.0d8  .and. temp .lt. 1.0d9) then
 tau   =  log10(temp * 1.0d-8)
 cc   =  1.5654d0
 c00  =  9.889d10
 c01  = -4.524d8
 c02  = -6.088d6
 c03  =  4.269d7
 c04  =  5.172d7
 c05  =  4.910d7
 c06  =  4.388d7
 c10  =  1.813d11
 c11  = -7.556d9
 c12  = -3.304d9
 c13  = -1.031d9
 c14  = -1.764d9
 c15  = -1.851d9
 c16  = -1.928d9
 c20  =  9.750d10
 c21  =  3.484d10
 c22  =  5.199d9
 c23  = -1.695d9
 c24  = -2.865d9
 c25  = -3.395d9
 c26  = -3.418d9
 dd01 = -1.135d8
 dd02 =  1.256d8
 dd03 =  5.149d7
 dd04 =  3.436d7
 dd05 =  1.005d7
 dd11 =  1.652d9
 dd12 = -3.119d9
 dd13 = -1.839d9
 dd14 = -1.458d9
 dd15 = -8.956d8
 dd21 = -1.549d10
 dd22 = -9.338d9
 dd23 = -5.899d9
 dd24 = -3.035d9
 dd25 = -1.598d9


else if (temp .ge. 1.0d9) then
 tau  =  log10(t9)
 cc   =  1.5654d0
 c00  =  9.581d10
 c01  =  4.107d8
 c02  =  2.305d8
 c03  =  2.236d8
 c04  =  1.580d8
 c05  =  2.165d8
 c06  =  1.721d8
 c10  =  1.459d12
 c11  =  1.314d11
 c12  = -1.169d11
 c13  = -1.765d11
 c14  = -1.867d11
 c15  = -1.983d11
 c16  = -1.896d11
 c20  =  2.424d11
 c21  = -3.669d9
 c22  = -8.691d9
 c23  = -7.967d9
 c24  = -7.932d9
 c25  = -7.987d9
 c26  = -8.333d9
 dd01 =  4.724d8
 dd02 =  2.976d8
 dd03 =  2.242d8
 dd04 =  7.937d7
 dd05 =  4.859d7
 dd11 = -7.094d11
 dd12 = -3.697d11
 dd13 = -2.189d11
 dd14 = -1.273d11
 dd15 = -5.705d10
 dd21 = -2.254d10
 dd22 = -1.551d10
 dd23 = -7.793d9
 dd24 = -4.489d9
 dd25 = -2.185d9
end if


taudt = iln10*tempi



!..equation 3.7, compute the expensive trig functions only one time
cos1 = cos(fac1*tau)
cos2 = cos(fac1*2.0d0*tau)
cos3 = cos(fac1*3.0d0*tau)
cos4 = cos(fac1*4.0d0*tau)
cos5 = cos(fac1*5.0d0*tau)
last = cos(fac2*tau)


sin1 = sin(fac1*tau)
sin2 = sin(fac1*2.0d0*tau)
sin3 = sin(fac1*3.0d0*tau)
sin4 = sin(fac1*4.0d0*tau)
sin5 = sin(fac1*5.0d0*tau)
xast = sin(fac2*tau)


a0 = 0.5d0*c00 &
     + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2 &
     + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4 &
     + c05*cos5 + dd05*sin5 + 0.5d0*c06*last


f0 =  taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0d0 &
     + dd02*cos2*2.0d0 - c03*sin3*3.0d0 + dd03*cos3*3.0d0 &
     - c04*sin4*4.0d0 + dd04*cos4*4.0d0 &
     - c05*sin5*5.0d0 + dd05*cos5*5.0d0) &
     - 0.5d0*c06*xast*fac2*taudt


a1 = 0.5d0*c10 &
     + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2 &
     + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4 &
     + c15*cos5 + dd15*sin5 + 0.5d0*c16*last


f1 = taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0d0 &
     + dd12*cos2*2.0d0 - c13*sin3*3.0d0 + dd13*cos3*3.0d0 &
     - c14*sin4*4.0d0 + dd14*cos4*4.0d0 - c15*sin5*5.0d0 &
     + dd15*cos5*5.0d0) - 0.5d0*c16*xast*fac2*taudt


a2 = 0.5d0*c20 &
     + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2 &
     + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4 &
     + c25*cos5 + dd25*sin5 + 0.5d0*c26*last


f2 = taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0d0 &
     + dd22*cos2*2.0d0 - c23*sin3*3.0d0 + dd23*cos3*3.0d0 &
     - c24*sin4*4.0d0 + dd24*cos4*4.0d0 - c25*sin5*5.0d0 &
     + dd25*cos5*5.0d0) - 0.5d0*c26*xast*fac2*taudt


!..equation 3.4
dum   = a0 + a1*zeta + a2*zeta2
dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0d0*a2*zeta*zetadt
dumdd = a1*zetadd + 2.0d0*a2*zeta*zetadd
dumda = a1*zetada + 2.0d0*a2*zeta*zetada
dumdz = a1*zetadz + 2.0d0*a2*zeta*zetadz


z      = exp(-cc*zeta)


xnum   = dum*z
xnumdt = dumdt*z - dum*z*cc*zetadt
xnumdd = dumdd*z - dum*z*cc*zetadd
xnumda = dumda*z - dum*z*cc*zetada
xnumdz = dumdz*z - dum*z*cc*zetadz


xden   = zeta3 + 6.290d-3*xlm1 + 7.483d-3*xlm2 + 3.061d-4*xlm3


dum    = 3.0d0*zeta2
xdendt = dum*zetadt - xldt*(6.290d-3*xlm2 &
         + 2.0d0*7.483d-3*xlm3 + 3.0d0*3.061d-4*xlm4)
xdendd = dum*zetadd
xdenda = dum*zetada
xdendz = dum*zetadz


dum      = 1.0d0/xden
fphot   = xnum*dum
fphotdt = (xnumdt - fphot*xdendt)*dum
fphotdd = (xnumdd - fphot*xdendd)*dum
fphotda = (xnumda - fphot*xdenda)*dum
fphotdz = (xnumdz - fphot*xdendz)*dum



!..equation 3.3
a0     = 1.0d0 + 2.045d0 * xl
xnum   = 0.666d0*a0**(-2.066d0)
xnumdt = -2.066d0*xnum/a0 * 2.045d0*xldt


dum    = 1.875d8*xl + 1.653d8*xl2 + 8.449d8*xl3 - 1.604d8*xl4
dumdt  = xldt*(1.875d8 + 2.0d0*1.653d8*xl + 3.0d0*8.449d8*xl2 &
         - 4.0d0*1.604d8*xl3)


z      = 1.0d0/dum
xden   = 1.0d0 + rm*z
xdendt =  -rm*z*z*dumdt
xdendd =  rmdd*z
xdenda =  rmda*z
xdendz =  rmdz*z


z      = 1.0d0/xden
qphot = xnum*z
qphotdt = (xnumdt - qphot*xdendt)*z
dum      = -qphot*z
qphotdd = dum*xdendd
qphotda = dum*xdenda
qphotdz = dum*xdendz


!..equation 3.2
sphot   = xl5 * fphot
sphotdt = 5.0d0*xl4*xldt*fphot + xl5*fphotdt
sphotdd = xl5*fphotdd
sphotda = xl5*fphotda
sphotdz = xl5*fphotdz


a1      = sphot
sphot   = rm*a1
sphotdt = rm*sphotdt
sphotdd = rm*sphotdd + rmdd*a1
sphotda = rm*sphotda + rmda*a1
sphotdz = rm*sphotdz + rmdz*a1


a1      = tfac4*(1.0d0 - tfac3 * qphot)
a2      = -tfac4*tfac3


a3      = sphot
sphot   = a1*a3
sphotdt = a1*sphotdt + a2*qphotdt*a3
sphotdd = a1*sphotdd + a2*qphotdd*a3
sphotda = a1*sphotda + a2*qphotda*a3
sphotdz = a1*sphotdz + a2*qphotdz*a3


if (sphot .le. 0.0) then
 sphot   = 0.0d0
 sphotdt = 0.0d0
 sphotdd = 0.0d0
 sphotda = 0.0d0
 sphotdz = 0.0d0
end if






!..bremsstrahlung neutrino section
!..for reactions like e- + (z,a) => e- + (z,a) + nu + nubar
!..                   n  + n     => n + n + nu + nubar
!..                   n  + p     => n + p + nu + nubar
!..equation 4.3


den6   = den * 1.0d-6
t8     = temp * 1.0d-8
t812   = sqrt(t8)
t832   = t8 * t812
t82    = t8*t8
t83    = t82*t8
t85    = t82*t83
t86    = t85*t8
t8m1   = 1.0d0/t8
t8m2   = t8m1*t8m1
t8m3   = t8m2*t8m1
t8m5   = t8m3*t8m2
t8m6   = t8m5*t8m1



tfermi = 5.9302d9*(sqrt(1.0d0+1.018d0*(den6*ye)**twoth)-1.0d0)


!.."weak" degenerate electrons only



if (temp .gt. 0.01d0 * tfermi) then


!..equation 5.3
 dum   = 7.05d6 * t832 + 5.12d4 * t83
 dumdt = (1.5d0*7.05d6*t812 + 3.0d0*5.12d4*t82)*1.0d-8


 z     = 1.0d0/dum
 eta   = rm*z
 etadt = -rm*z*z*dumdt
 etadd = rmdd*z
 etada = rmda*z
 etadz = rmdz*z


 etam1 = 1.0d0/eta
 etam2 = etam1 * etam1
 etam3 = etam2 * etam1



!..equation 5.2
 a0    = 23.5d0 + 6.83d4*t8m2 + 7.81d8*t8m5
 f0    = (-2.0d0*6.83d4*t8m3 - 5.0d0*7.81d8*t8m6)*1.0d-8
 xnum  = 1.0d0/a0


 dum   = 1.0d0 + 1.47d0*etam1 + 3.29d-2*etam2
 z     = -1.47d0*etam2 - 2.0d0*3.29d-2*etam3
 dumdt = z*etadt
 dumdd = z*etadd
 dumda = z*etada
 dumdz = z*etadz


 c00   = 1.26d0 * (1.0d0+etam1)
 z     = -1.26d0*etam2
 c01   = z*etadt
 c02   = z*etadd
 c03   = z*etada
 c04   = z*etadz

 z      = 1.0d0/dum
 xden   = c00*z
 xdendt = (c01 - xden*dumdt)*z
 xdendd = (c02 - xden*dumdd)*z
 xdenda = (c03 - xden*dumda)*z
 xdendz = (c04 - xden*dumdz)*z


 fbrem   = xnum + xden
 fbremdt = -xnum*xnum*f0 + xdendt
 fbremdd = xdendd
 fbremda = xdenda
 fbremdz = xdendz



!..equation 5.9
 a0    = 230.0d0 + 6.7d5*t8m2 + 7.66d9*t8m5
 f0    = (-2.0d0*6.7d5*t8m3 - 5.0d0*7.66d9*t8m6)*1.0d-8


 z     = 1.0d0 + rm*1.0d-9
 dum   = a0*z
 dumdt = f0*z
 z     = a0*1.0d-9
 dumdd = z*rmdd
 dumda = z*rmda
 dumdz = z*rmdz


 xnum   = 1.0d0/dum
 z      = -xnum*xnum
 xnumdt = z*dumdt
 xnumdd = z*dumdd
 xnumda = z*dumda
 xnumdz = z*dumdz


 c00   = 7.75d5*t832 + 247.0d0*t8**(3.85d0)
 dd00  = (1.5d0*7.75d5*t812 + 3.85d0*247.0d0*t8**(2.85d0))*1.0d-8


 c01   = 4.07d0 + 0.0240d0 * t8**(1.4d0)
 dd01  = 1.4d0*0.0240d0*t8**(0.4d0)*1.0d-8


 c02   = 4.59d-5 * t8**(-0.110d0)
 dd02  = -0.11d0*4.59d-5 * t8**(-1.11d0)*1.0d-8


 z     = den**(0.656d0)
 dum   = c00*rmi  + c01  + c02*z
 dumdt = dd00*rmi + dd01 + dd02*z
 z     = -c00*rmi*rmi
 dumdd = z*rmdd + 0.656d0*c02*den**(-0.454d0)
 dumda = z*rmda
 dumdz = z*rmdz


 xden  = 1.0d0/dum
 z      = -xden*xden
 xdendt = z*dumdt
 xdendd = z*dumdd
 xdenda = z*dumda
 xdendz = z*dumdz


 gbrem   = xnum + xden
 gbremdt = xnumdt + xdendt
 gbremdd = xnumdd + xdendd
 gbremda = xnumda + xdenda
 gbremdz = xnumdz + xdendz



!..equation 5.1
 dum    = 0.5738d0*zbar*ye*t86*den
 dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
 dumdd  = 0.5738d0*zbar*ye*t86
 dumda  = -dum*abari
 dumdz  = 0.5738d0*2.0d0*ye*t86*den


 z       = tfac4*fbrem - tfac5*gbrem
 sbrem   = dum * z
 sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
 sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
 sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
 sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)





!..liquid metal with c12 parameters (not too different for other elements)
!..equation 5.18 and 5.16


else
 u     = fac3 * (log10(den) - 3.0d0)
 a0    = iln10*fac3*deni


!..compute the expensive trig functions of equation 5.21 only once
 cos1 = cos(u)
 cos2 = cos(2.0d0*u)
 cos3 = cos(3.0d0*u)
 cos4 = cos(4.0d0*u)
 cos5 = cos(5.0d0*u)


 sin1 = sin(u)
 sin2 = sin(2.0d0*u)
 sin3 = sin(3.0d0*u)
 sin4 = sin(4.0d0*u)
 sin5 = sin(5.0d0*u)


!..equation 5.21
 fb =  0.5d0 * 0.17946d0  + 0.00945d0*u + 0.34529d0 &
       - 0.05821d0*cos1 - 0.04969d0*sin1 &
       - 0.01089d0*cos2 - 0.01584d0*sin2 &
       - 0.01147d0*cos3 - 0.00504d0*sin3 &
       - 0.00656d0*cos4 - 0.00281d0*sin4 &
       - 0.00519d0*cos5


 c00 =  a0*(0.00945d0 &
       + 0.05821d0*sin1       - 0.04969d0*cos1 &
       + 0.01089d0*sin2*2.0d0 - 0.01584d0*cos2*2.0d0 &
       + 0.01147d0*sin3*3.0d0 - 0.00504d0*cos3*3.0d0 &
       + 0.00656d0*sin4*4.0d0 - 0.00281d0*cos4*4.0d0 &
       + 0.00519d0*sin5*5.0d0)



!..equation 5.22
 ft =  0.5d0 * 0.06781d0 - 0.02342d0*u + 0.24819d0 &
       - 0.00944d0*cos1 - 0.02213d0*sin1 &
       - 0.01289d0*cos2 - 0.01136d0*sin2 &
       - 0.00589d0*cos3 - 0.00467d0*sin3 &
       - 0.00404d0*cos4 - 0.00131d0*sin4 &
       - 0.00330d0*cos5


 c01 = a0*(-0.02342d0 &
       + 0.00944d0*sin1       - 0.02213d0*cos1 &
       + 0.01289d0*sin2*2.0d0 - 0.01136d0*cos2*2.0d0 &
       + 0.00589d0*sin3*3.0d0 - 0.00467d0*cos3*3.0d0 &
       + 0.00404d0*sin4*4.0d0 - 0.00131d0*cos4*4.0d0 &
       + 0.00330d0*sin5*5.0d0)



!..equation 5.23
 gb =  0.5d0 * 0.00766d0 - 0.01259d0*u + 0.07917d0 &
       - 0.00710d0*cos1 + 0.02300d0*sin1 &
       - 0.00028d0*cos2 - 0.01078d0*sin2 &
       + 0.00232d0*cos3 + 0.00118d0*sin3 &
       + 0.00044d0*cos4 - 0.00089d0*sin4 &
       + 0.00158d0*cos5


 c02 = a0*(-0.01259d0 &
       + 0.00710d0*sin1       + 0.02300d0*cos1 &
       + 0.00028d0*sin2*2.0d0 - 0.01078d0*cos2*2.0d0 &
       - 0.00232d0*sin3*3.0d0 + 0.00118d0*cos3*3.0d0 &
       - 0.00044d0*sin4*4.0d0 - 0.00089d0*cos4*4.0d0 &
       - 0.00158d0*sin5*5.0d0)



!..equation 5.24
 gt =  -0.5d0 * 0.00769d0  - 0.00829d0*u + 0.05211d0 &
       + 0.00356d0*cos1 + 0.01052d0*sin1 &
       - 0.00184d0*cos2 - 0.00354d0*sin2 &
       + 0.00146d0*cos3 - 0.00014d0*sin3 &
       + 0.00031d0*cos4 - 0.00018d0*sin4 &
       + 0.00069d0*cos5


 c03 = a0*(-0.00829d0 &
       - 0.00356d0*sin1       + 0.01052d0*cos1 &
       + 0.00184d0*sin2*2.0d0 - 0.00354d0*cos2*2.0d0 &
       - 0.00146d0*sin3*3.0d0 - 0.00014d0*cos3*3.0d0 &
       - 0.00031d0*sin4*4.0d0 - 0.00018d0*cos4*4.0d0 &
       - 0.00069d0*sin5*5.0d0)



 dum   = 2.275d-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
 dumdt = -dum*tempi
 dumdd = oneth*dum*deni
 dumda = -oneth*dum*abari
 dumdz = 2.0d0*dum*zbari

 gm1   = 1.0d0/dum
 gm2   = gm1*gm1
 gm13  = gm1**oneth
 gm23  = gm13 * gm13
 gm43  = gm13*gm1
 gm53  = gm23*gm1



!..equation 5.25 and 5.26
 v  = -0.05483d0 - 0.01946d0*gm13 + 1.86310d0*gm23 - 0.78873d0*gm1
 a0 = oneth*0.01946d0*gm43 - twoth*1.86310d0*gm53 + 0.78873d0*gm2


 w  = -0.06711d0 + 0.06859d0*gm13 + 1.74360d0*gm23 - 0.74498d0*gm1
 a1 = -oneth*0.06859d0*gm43 - twoth*1.74360d0*gm53 + 0.74498d0*gm2



!..equation 5.19 and 5.20
 fliq   = v*fb + (1.0d0 - v)*ft
 fliqdt = a0*dumdt*(fb - ft)
 fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0d0 - v)*c01
 fliqda = a0*dumda*(fb - ft)
 fliqdz = a0*dumdz*(fb - ft)


 gliq   = w*gb + (1.0d0 - w)*gt
 gliqdt = a1*dumdt*(gb - gt)
 gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0d0 - w)*c03
 gliqda = a1*dumda*(gb - gt)
 gliqdz = a1*dumdz*(gb - gt)



!..equation 5.17
 dum    = 0.5738d0*zbar*ye*t86*den
 dumdt  = 0.5738d0*zbar*ye*6.0d0*t85*den*1.0d-8
 dumdd  = 0.5738d0*zbar*ye*t86
 dumda  = -dum*abari
 dumdz  = 0.5738d0*2.0d0*ye*t86*den


 z       = tfac4*fliq - tfac5*gliq
 sbrem   = dum * z
 sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
 sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
 sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
 sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)


end if





!..recombination neutrino section
!..for reactions like e- (continuum) => e- (bound) + nu_e + nubar_e
!..equation 6.11 solved for nu
xnum   = 1.10520d8 * den * ye /(temp*sqrt(temp))
xnumdt = -1.50d0*xnum*tempi
xnumdd = xnum*deni
xnumda = -xnum*abari
xnumdz = xnum*zbari


!..the chemical potential
nu   = ifermi12(xnum)


!..a0 is d(nu)/d(xnum)
a0 = 1.0d0/(0.5d0*zfermim12(nu))
nudt = a0*xnumdt
nudd = a0*xnumdd
nuda = a0*xnumda
nudz = a0*xnumdz


nu2  = nu * nu
nu3  = nu2 * nu


!..table 12
if (nu .ge. -20.0  .and. nu .lt. 0.0) then
 a1 = 1.51d-2
 a2 = 2.42d-1
 a3 = 1.21d0
 b  = 3.71d-2
 c  = 9.06e-1
 d  = 9.28d-1
 f1 = 0.0d0
 f2 = 0.0d0
 f3 = 0.0d0
else if (nu .ge. 0.0  .and. nu .le. 10.0) then
 a1 = 1.23d-2
 a2 = 2.66d-1
 a3 = 1.30d0
 b  = 1.17d-1
 c  = 8.97e-1
 d  = 1.77d-1
 f1 = -1.20d-2
 f2 = 2.29d-2
 f3 = -1.04d-3
end if



!..equation 6.7, 6.13 and 6.14
if (nu .ge. -20.0  .and.  nu .le. 10.0) then


 zeta   = 1.579d5*zbar*zbar*tempi
 zetadt = -zeta*tempi
 zetadd = 0.0d0
 zetada = 0.0d0
 zetadz = 2.0d0*zeta*zbari


 c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
 c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
 dum    = zeta*c00
 dumdt  = zetadt*c00 + zeta*c01*nudt
 dumdd  = zeta*c01*nudd
 dumda  = zeta*c01*nuda
 dumdz  = zetadz*c00 + zeta*c01*nudz



 z      = 1.0d0/dum
 dd00   = dum**(-2.25)
 dd01   = dum**(-4.55)
 c00    = a1*z + a2*dd00 + a3*dd01
 c01    = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z



 z      = exp(c*nu)
 dd00   = b*z*(1.0d0 + d*dum)
 gum    = 1.0d0 + dd00
 gumdt  = dd00*c*nudt + b*z*d*dumdt
 gumdd  = dd00*c*nudd + b*z*d*dumdd
 gumda  = dd00*c*nuda + b*z*d*dumda
 gumdz  = dd00*c*nudz + b*z*d*dumdz



 z   = exp(nu)
 a1  = 1.0d0/gum


 bigj   = c00 * z * a1
 bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
 bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
 bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
 bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz



!..equation 6.5
 z     = exp(zeta + nu)
 dum   = 1.0d0 + z
 a1    = 1.0d0/dum
 a2    = 1.0d0/bigj


 sreco   = tfac6 * 2.649d-18 * ye * zbar**13 * den * bigj*a1
 srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
 srecodd = sreco*(1.0d0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
 srecoda = sreco*(-1.0d0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
 srecodz = sreco*(14.0d0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)


end if



!..convert from erg/cm^3/s to erg/g/s
!..comment these out to duplicate the itoh et al plots


spair   = spair*deni
spairdt = spairdt*deni
spairdd = spairdd*deni - spair*deni
spairda = spairda*deni
spairdz = spairdz*deni


splas   = splas*deni
splasdt = splasdt*deni
splasdd = splasdd*deni - splas*deni
splasda = splasda*deni
splasdz = splasdz*deni


sphot   = sphot*deni
sphotdt = sphotdt*deni
sphotdd = sphotdd*deni - sphot*deni
sphotda = sphotda*deni
sphotdz = sphotdz*deni


sbrem   = sbrem*deni
sbremdt = sbremdt*deni
sbremdd = sbremdd*deni - sbrem*deni
sbremda = sbremda*deni
sbremdz = sbremdz*deni


sreco   = sreco*deni
srecodt = srecodt*deni
srecodd = srecodd*deni - sreco*deni
srecoda = srecoda*deni
srecodz = srecodz*deni



!..the total neutrino loss rate
snu    =  splas + spair + sphot + sbrem + sreco
dsnudt =  splasdt + spairdt + sphotdt + sbremdt + srecodt
dsnudd =  splasdd + spairdd + sphotdd + sbremdd + srecodd
dsnuda =  splasda + spairda + sphotda + sbremda + srecoda
dsnudz =  splasdz + spairdz + sphotdz + sbremdz + srecodz



return
end subroutine sneut

end module nuclear_lib
