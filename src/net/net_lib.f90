!----------------------------------------------------------------------
! net_lib
!----------------------------------------------------------------------
! Aggregates YREC's COMMON-free nuclear-physics utility routines into a
! single module, following the numerics_lib/const_lib convention (see
! GUIDELINES.md, and MESA's own one small public *_lib module per
! domain). Each routine below was originally its own file
! (nuclear/<name>.f90, with per-file conversion notes preserved below);
! bodies are unchanged, only concatenated under one module. Callers use
! this via `use net_lib`. nulosses.f90 contained two subroutines
! (nulosses, azbar); both are included here.
!
! lir and ratext, originally also housed in nuclear/ (and formerly
! part of this module), moved out to numerics_lib (2026): both are
! generic numerical utilities with no nuclear-physics content -- they
! were only filed here because this module's own liburn (below) happens
! to be one of their callers. See numerics_lib.f90 and GUIDELINES.md's rule that
! folder/module placement should track function, not caller.
!
! Completed the merge (2026, phase-two reorg -- see GUIDELINES.md's
! "Physics domains still entangled with the solver"): dburn.f90,
! dburnm.f90 (a genuine near-duplicate of dburn, not collapsed into
! it -- different rate-data source, timestep units, convergence
! threshold, and accretion-weighting formula; investigated and kept
! separate), deutrate.f90, engeb.f90, liburn.f90 (which itself
! bundled a second small subroutine, safedivexp), liburn2.f90,
! lirate88.f90, and the two private Fermi-Dirac integral functions
! ifermi12.f90/zfermim12.f90 (previously relocated here from kap/
! during that domain's sweep; called only from this module's own
! `rates`, so private to it) all moved in unchanged. engeb.f90's own
! `use net_lib` (it called `neutrino`) was removed -- a module
! cannot use itself; host association already gives every contained
! procedure access to every sibling procedure. Callers that didn't
! already `use net_lib` (setup/midmod.f90, misc/coefft.f90,
! core/main.f90, util/ytime.f90, rotation/getw.f90,
! mixing/bursmix.f90) had it added.
! ---------------------------------------------------------------------
! MAP OF THE MODULE (2026; updated at the physics-purity split).
! net_lib now holds ONLY the pure kernels -- functions of
! (logRho, logT, composition, controls), no model state, the surface
! test_nuclear pins:
!     rates      13 reaction rates + branching fractions at one point
!                (consumes cross_section_scale from setup/remap.f90).
!     sneut      Itoh et al. (1996) neutrino losses (pure fits).
!     nulosses   azbar + sneut wrapper for a full mixture.
!     neutrino   legacy interface to the loss rates.
!     azbar      abar/zbar/molar-abundance bookkeeping.
!     safedivexp division-with-exponent guard.
!     ifermi12, zfermim12   Fermi-Dirac integral inverses (rates).
! The BURN DRIVERS (engeb/eqburn, dburn/dburnm/deutrate,
! liburn/liburn2/lirate88) -- everything that reads or writes
! star_info -- moved to core/burn_lib.f90: they advance the model and
! are star-layer code (MESA: struct_burn_mix, not net).
! ---------------------------------------------------------------------
module net_lib
contains

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
! ifermi12/zfermim12 are no longer declared as external doubles here:
! since both moved into this module (2026, see net_lib.f90's own
! header), declaring their names again as plain locals would shadow
! the host-associated module functions rather than resolve to them.
double precision nu,nudt,nudd,nuda,nudz, &
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

!----------------------------------------------------------------------
! rates
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rates.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model). Numeric literals (including
! ones missing the customary D0/D-3/etc suffix, e.g. many of the
! DATA-statement values below) are copied character-for-character from
! the original -- do not "correct" them, see project notes.
!
! JULY 3, 1991 (MHP)
! THIS SUBROUTINE COMPUTES THE NUCLEAR BURNING RATES FOR USE IN KEMCOM.
! IT IS A STRIPPED-DOWN VERSION OF ENGEB; SEE ENGEB FOR DETAILED NOTES
! ON THE REACTIONS.
! ALL PHYSICAL QUANTITIES ARE CGS. THE REACTION RATES ARE IN
!  GM^{-1}S^{-1}.
!  MHP 10/97
!
! On 10/13/97, JNB converted the nuclear masses from neutral nuclear
! masses to bare nuclear masses by subtracting Z(I)*(m_e)*c^2 from
! the neutral nuclear masses. This caused changes in a number of
! places: in ANUC(I), in Q1(I)-Q7(I), in the calculation of the Be7electron
! Be7proton rates, and in the calculation of the N15p branching ratio.
!
! JNB made some purely cosmetic changes on 1/20/96. Revised all of the
! input Q's on 9/23-25/97 to agree with submitted version of Solar
! Fusion Workshop paper. The SStandard are fixed to agree with
! the Workshop paper. JNB recalculated all of the EG(I) to determine
! the best values for the energy generation for all the reactions,
! taking account of my improved calculations of neutrino energy loss.
! The calculations are documented in Vol. 19, 132-141, 1997 of my notes.
! The neutral atom mass differences are taken from Table of Isotopes,
! 8th Ed, 1996 and the neutrino energy losses from Bahcall, Gallium
! solar neutrino experiments, Phys Rev C, in press, 1997.
!
! ALL NUMBERS IN THIS SUBROUTINE HAVE BEEN CHECKED AND REVISED, WHERE
!  NECESSARY, BY JOHN BAHCALL SO THAT THEY AGREE WITH THE MODERN NUMBERS
!  IN NEUTRINO ASTROPHYSICS (1989).
! *****************************************
! CHANGING NUCLEAR REACTION CROSS SECTIONS.
! *****************************************
! NUCLEAR CROSS SECTIONS CAN BE CHANGED SIMPLY BY INSERTING NEW NUMBERS
!  FOR THE DATA VALUES OF SSTANDARD(I), WHICH ARE THE RATIOS OF THE
!  DESIRED CROSS SECTION FACTORS TO THE VALUES GIVEN IN TABLES 3.2 AND
!  3.4 OF NEUTRINO ASTROPHYSICS (1989). IF THE VALUE GIVEN IN NEUTRINO
!  ASTROPHYSICS IS USED, THEN SSTANDARD(I) = 1.0 .  TO INCREASE THE CROSS
!  SECTION FOR REACTION K BY A FACTOR OF TWO COMPARED TO THE STANDARD
!  VALUE, SET SSTANDARD(K) = 2.0 .  TO DETERMINE WHICH VALUE OF I GOES
!  WITH WHICH REACTION, SEE THE SECTION JUST BELOW.
! THE ENERGY DERIVATIVES ENTER IN A FORM IN WHICH THEY ARE DIVIDED BY
!  THE ABSOLUTE VALUES OF THE CROSS SECTIONS AT ZERO ENERGY.  THUS IF
!  THE SHAPE OF THE CROSS SECTION EXTRAPOLATION IS UNCHANGED AND ONLY
!  THE INTERCEPT OF S(E) AT ZERO ENERGY IS CHANGED, THEN NO CORRECTION
!  NEED BE MADE FOR THE DERIVATIVES.  THEY ARE AUTOMATICALLY SCALED
!  CORRECTLY.  THE EXACT WAY THAT THE DERIVATIVES ENETER THE
!  CALCUALTIONS IS DESCRIBED IN THE SECTION LABELED ``DEFINING THE
!  Q(I)' THAT IS PRESENTED BELOW.
! ************************************
! IDENTIFYING THE REACTIONS.
! ************************************
!  THE VALUE OF J DENOTES WHICH OF THE REACTIONS THE COEFFICIENTS
!  REFER TO:
!  J = 1, PP; J = 2, HE3+HE3; J = 3, HE3+ HE4; J =4, P + C12;  J = 5, P+C13;
!  J = 6. P + N14; J = 7, P + O16 ; J = 8, HE4+C13; J = 10, HE4+C12;
!  J = 11, HE4 + N14; J = 12, TRIPLE ALPHA.
!  REACTION 14 IS PEP; REACTION 15 IS BE7 ELECTRON CAPTURE; REACTION 16 IS
!   BE7 PROTON CAPTURE; REACTION 17 IS THE HEP REACTION.
!  do not change sstandard(14) unless you want to change the ratio of
!  pep to pp.
!   REACTIONS 14-17 WERE NOT EXPLICITLY INCLUDED IN THE YALE
!    PREVIOUS VERSION OF THE CODE, BUT THEY ARE MOST OF THE STORY FOR
!    THE SOLAR NEUTRINO PROBLEM.
!   THE BRANCHING OF THE N15 + P REACTIONS IS TREATED IN A SERIES OF
!    SEPARATE STATEMENTS FOLLOWING THE CALCULATION OF THE BE7 + P
!    REACTION. SEE THE DEFINITIONS OF F3 AND F4.  IF THE CROSS-SECTION
!    FACTORS OF THE N15 + P REACTIONS ARE REVISED, THEN THE NUMERICAL
!    COEFFICIENTS MUST BE CHANGED IN THE DEFINITION OF C12ALPHA AND
!    O16GAMMA.
! FOR Q1(I), ...,Q(5(I), I = 8 CORRESPONDS TO THE BE7 +  P REACTION.
!  THIS ASSIGNMENT FOR I = 8 IS ONLY VALID FOR THE LISTED Q'S AND NOT
!  FOR OTHER ARRAYS IN THE PROGRAM.
! IU IS THE SHELL NUMBER.
subroutine rates(log_density,log_temperature,hydrogen_fraction, &
     helium_fraction,he3_fraction,c12_fraction,c13_fraction,n14_fraction, &
     o16_fraction,o18_fraction,zone_idx,rate_pp,rate_he3_he3,rate_he3_he4, &
     rate_c12_p,rate_c13_p,rate_n14_p,rate_o16_p,rate_c13_alpha,rate_zero9, &
     rate_c12_alpha,rate_n14_alpha,rate_triple_alpha,rate_zero13, &
     frac_c12_alpha,frac_be7_electron)

      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: log_density, log_temperature, &
           hydrogen_fraction, helium_fraction, he3_fraction, c12_fraction, &
           c13_fraction, n14_fraction, o16_fraction, o18_fraction
      integer, intent(in) :: zone_idx
      double precision, intent(out) :: rate_pp(json), rate_he3_he3(json), &
           rate_he3_he4(json), rate_c12_p(json), rate_c13_p(json), &
           rate_n14_p(json), rate_o16_p(json), rate_c13_alpha(json), &
           rate_zero9(json), rate_c12_alpha(json), rate_n14_alpha(json), &
           rate_triple_alpha(json), rate_zero13(json)
      double precision, intent(out) :: frac_c12_alpha(json), &
           frac_be7_electron(json)





      double precision :: mass_frac(13), rate(13), screening_factor(13), &
           charge_product(13), z53(13), z43(13), z23(13), z86(13)
      double precision :: q1(8), q2(8), q3(8), q4(8), q5(8), q6(7), q7(7), &
           q8(7)
      double precision :: atomic_mass(13), atomic_charge(13)
      integer :: num_isotopes, num_reactions

! ***************************
! ANUC ARE ATOMIC MASS UNITS.
! ***************************
! On 10/13/97, JNB converted ANUC(I) from neutral atomic masses to
! bare nuclear masses.
!
! The scale is the mass of C12 divided by 12 or 931.49432 MeV,
! which is 1.6605402 times 10^{-24} gm. The values are obtained
! by dividing the mass excess (expressed in MeV) by 931.49432 MeV and
! adding to this the atomic mass number, A.  The value for Be7, which
! is used implicitly in this subroutine was, 7.016930, until 13/13/97.
! We now use the bare nuclear mass of 7.014735 .
!
      data atomic_mass/1.008665,1.007276,2.013553,3.015501,3.014933,4.001506, &
           11.996709,13.000064,13.999233,15.990526,17.994772,19.986954, &
           23.978458/,atomic_charge/0.,1.,1.,1.,2.,2.,6.,6.,7.,8.,8.,10.,12./, &
           num_isotopes/13/
! THE ISOTOPES ARE NEUTRON, H1, D, H3, HE3, HE4, C12, C13, N14, O16, O18,
!  NE20, MG24, RESPECTIVELY. ALL OF THESE NUMBERS WERE CHECKED.
! NELEM IS THE NUMBER OF ISOTOPES INCLUDED.
! **************************************************************************
! THE QUANTITIES Q1(J), Q2(J), ...,Q5(J) ARE THE TERMS IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS AND IN EQUATION 53 OF FOWLER, CAUGHLAN, AND
!  ZIMMERMAN (1967) EQ. 53, IN BOTH CASES MULTIPLIED BY T SUB 9 ^(-2/3).
!  THE REACTIONS CORRESPONDING TO EACH J ARE LISTED ABOVE, UNDER:
!   IDENTIFYING THE RECTIONS.
!  FOR THIS SET OF PARAMETERS, AND ONLY FOR THIS SET OF PARAMETERS,
!  J = 8 CORRESPONDS TO THE BE7 + P REACTION.
! **************************************************************************
! GENERAL EXPRESSION FOR Q'S:
! T_9^(-2/3)[S_EFF/S(0)] =
!   [T_9^(-2/3) + Q1(I)T_9^(-1/3) + Q2(I) + Q3(I)T_9^(1/3) +
!    Q4(I)T_9^(2/3) + Q(5)T_9 ]
!
! BY COMPARISON WITH EQUATION 3.14 OF NEUTRINO ASTROPHYSICS, WE SEE THAT
!
!  Q1 = (5/(12*TAU))*T_9^(-1/3) .
!  Q2 = (S'/S)(E_0))*T_9^(-2/3) .
!  Q3 = (S'/S)(35/36)(K*10^9 K)
!  Q4 = (S''/2S)(E_0^2)(T_9^(-4/3)
!  Q5 = (89/72)(S''/S)(E_0)(KT)(T_9^(-5/3)
!
! EACH OF THE Q'S IS INDEPENDENT OF TEMPERATURE (T), AS CAN BE SEEN FROM
!  EQUATIONS 3.10 AND 3.11 .
!
! ALL OF THE VALUES OF THE Q1, ...., Q5 HAVE BEEN RECALCULATED, USING
!  WHERE NEEDED NUCLEAR CROSS SECTIONS GIVEN IN TABLES 3.2 AND 3.4 OF NEUTRINO
!  ASTROPHYSICS.
! ******************************************************************
! Q6 IS THE COEFFICIENT OF THE TEMPERATURE TERM IN THE DEFINITION OF
!  TAU, EQUATION 3.10 OF NEUTRINO ASTROPHYSICS.
!  TAU = Q6*(T SUB 9 TO THE (-1/3) POWER ).
! ******************************************************************
! SLIGHT CHANGES HAVE BEEN MADE IN THE PREVIOUS VALUES OF Q6 TO MAKE
!  THE DATA MORE ACCURATE.
! NOTE Q6 IS NEGATIVE.
! ********************************************************************
! Q7 IS THE CONSTANT IN FRONT OF THE REACTION RATE. THE LOGARITHMS ARE
!  RELATED TO S FACTORS BY A LOGARITHM AND VARIOUS NUMERICAL FACTORS.
! ********************************************************************
! THE GENERAL RELATION IS:
!  Q7 = 70.62860 + ((LN(Z0*Z1/A))/3) -LN(A0*A1) - LN(S SUB 0)
!       -LN(1 + DELTA_01)
!  HERE S SUB 0 IS THE CROSS SECTION FACTOR IN UNITS OF KEV-BARNS.
!  THE NUMERICAL VALUES USED HERE ARE TAKEN FROM TABLES 3.4 AND 3.2
!   OF NUCLEAR ASTROPHYSICS.
!  THE QUANTITY DELTA_01 IS NON-ZERO (EQUAL TO UNITY) ONLY WHEN THE
!   TWO REACTING NUCLEI, 0 AND 1, ARE IDENTICAL.
! Q8 REFLECTS A TERM IN THE EXPONENTIAL THAT OCCURS IN THE RATES, THE
!  TERM BEING PROPORTIONAL TO E^( CONSTANT*T_9^2). SEE HARRIS, ET AL.
!  1983, ANNUAL REV. ASTRON. ASTROPHYS.21, 165 (1983) FOR THE MEANING
!  OF THIS OBSCURE TERM.
! THE VALUES OF Q1(I),...Q7(I) GIVEN BELOW WERE OBTAINED USING THE DATA
!  IN NEUTRINO ASTROPHYSICS WITH THE AID OF AUXILARY COMPUTER CODES
!  THAT GENERATED ACCURATE EVALUATIONS. THE NUMBERS GIVEN HERE ARE
!  IN MANY CASES COMPLETELY DIFFERENT FROM THE VALUES IN THE ORIGINAL
!  YALE CODE.
! NRXNS IS THE NUMBER OF REACTIONS BEING TRACKED.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! 10/13/97. Changed Q(I) so that are now calculated for the bare nuclear
! masses. The calculates were made cues.f . Last date on the calculations
! in my notes is 10/10/97.  Note Q6 is -tau(T_9 = 1); in cues.f, I calculate
! tau(T_6 = 1). The connection is:  Q6 = 0.1*tau(cues.f).
!
!
! BP00 values
!      DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
!     1 Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, -0.2095/,
!     2 Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/,
!     3 Q4/0.,0.,0.,4.79,7.456,0.0,.69703,.16762/,
!     4 Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,.12114/,
!     5 Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,
!     5    -16.6925/,
!     6 Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7 Q8/0.,0.,0.,0.0,0.0,0.0,0.0/,
!     8NRXNS/13/
! Revised values 7/21/03
!      DATA Q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/,
!     1Q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224, 0.0/,
!     2Q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,0.0/,
!     3Q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.0/,
!     4Q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.0/,
!     5Q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/,
!     6Q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/,
!     7Q8/0.,0.,0.,0.0,0.0,0.0,0./,
!     8NRXNS/13/
! MHP 9/14 RESTORED BE7+P DERIVATIVES THAT WERE ZEROED OUT
      data q1/0.12317,.03392,.0325,.0304,.03035,.0273,.02494,.040572/, &
      q2/1.08749,-.273,-.2085,.7630,-0.4044,-1.60,-1.224,-0.2095/, &
      q3/.93833,-.0648,-.0474,.1626,-.08598,-.3064,-.2139,-0.0595/, &
      q4/0.,0.,0.,4.79,7.456,0.0,.69703,0.16762/, &
      q5/0.,0.,0.,2.595,4.032,0.0,0.3097,0.12114/, &
      q6/-3.3804,-12.2757,-12.826,-13.6899,-13.7173,-15.2281,-16.6925/, &
      q7/20.8964,76.6003,67.8036,69.130,70.3809,69.8517,70.8012/, &
      q8/0.,0.,0.,0.0,0.0,0.0,0./, &
      num_reactions/13/
! For different values of SSTANDARD, check the comments in ENGEB
! *********************************************************************
! THE VALUES OF SSTANDARD(I) ARE TO BE CHANGED FROM UNITY IF THE CROSS
!  SECTION FACTORS ARE NOT THE ONES GIVEN IN NEUTRINO ASTROPHYSICS,
!  TABLE 3.2 AND TABLE 3.4 .
! THE CURRENT VALUES ARE CHANGED TO BE THE PREFERRED VALUES LISTED IN
!  THE LAST COLUMN OF TABLE 1 OF BAHCALL AND PINSONNEAULT (1991),
!   6/14/91.
!  THE VALUE OF SSTANDARD(17) FOR HEP IS TAKEN FROM WOLFS ET AL. , PHYS.
!   REV. LETTERS, 63, 2721 (1989).  IT CORRESPONDS TO AN S SUB 0 =
!   15.3E-20 KEV-BARNS, A FACTOR OF 1.913 LARGER THAN FOR THE OLDER
!   MEASUREMENTS USED IN NEUTRINO ASTROPHYSICS.
! *********************************************************************
! ZPRD IS USED IN SCREENING CALCULATIONS. IT IS THE PRODUCT OF THE
!  CHARGES OF THE INTERACTING IONS. ZPRD WAS CHECKED. Z86 IS USED
!  IN CALCULATING INTERMEDIATE SCREENING AND IS DEFINED BY GRABOSKE ET
!  AL, AP. J. 181, PAGE 465 (1973), IN TABLE 4.  Z86 WAS CHECKED AND
!  SOME NUMERICAL VALUES WERE MADE SLIGHTLY MORE ACCURATE.  Z53, Z43,
!  AND Z23 ARE ALSO DEFINED IN TABLE 4 (SEE ABOVE).  SINCE THEY ARE
!  ONLY USED IN STRONG SCREENING, THE VALUES OF Z53, Z43, AND Z23
!  WERE NOT CHECKED.
      double precision :: c21
      data charge_product/1.,4.,4.,6.,6.,7.,8.,12.,16.,12.,14.,12.,36./,z53/ &
         1.175,3.73,3.73,4.804,4.804,5.385,5.941,9.014,11.24,9.014, &
         10.15,9.104,23.28/,z43/0.52,1.31,1.31,1.488,1.488,1.61,1.721, &
         2.577,3.025,2.577,2.81,2.577,5.668/,z23/-0.413,-0.655,-0.655, &
         -0.643,-0.643,-0.659,-0.673,-0.889,-0.946,-0.889,-0.92,-0.889, &
         -1.36/,z86/1.630,5.917,5.917,8.302,8.302,9.520,10.716,16.192, &
         20.978,16.192,18.606,16.192,45.6635/, &
         c21/5.240358E-8/
! DEFINE NEXT THE FRACTIONAL ABUNDANCES BY MASS OF THE IMPORTANT
!  ISOTOPES.
! X, Y, Z, XHE3,..., XBE9 ARE THE MASS FRACTIONS OF THE ISOTOPES.
!  THE ABUNDANCES OF NEUTRONS, H2, H3, NE20,AND MG24, WHICH ARE,
!  RESPECTIVELY, XFRAC(I) FOR I = 1,3,4,12,13, ARE NO LONGER USED.

      integer :: i, nz
      double precision :: mu_ion_inv, mu_e_inv, xtr, zeta0, trm
      double precision :: rho_over_mu_e, log_rho_local, density, t9, t9p13, &
           t9p23, t9m13, t9m23, t9m1, t9m2, t9m12, t9m32
      double precision :: fermi_mom_sq, fermi_energy_over_kt, f_prime_over_f, &
           log_degeneracy
      double precision :: lambda0, lambda0_23, lambda0_86, z_curl, z_bar, &
           z_curl_58, z_bar_28, z_bar_13, lambda0_zcurl
      double precision :: weak_screening_u, intermediate_screening_u, &
           strong_screening_u
      double precision :: r1, r2, a1, a2, a3, a4, a5
      double precision :: be7_electron_rate, be7_proton_rate, be7_temp_factor, &
           be7_q_factor, be7p_charge_product, be7p_z86, be7p_screening_u, &
           be7_mass_factor
      double precision :: be7_electron_frac, c12_alpha_frac, o16_gamma_frac
      double precision :: o16_gamma_rate, c12_alpha_n15p_rate

      mass_frac(1) = 0.0
      mass_frac(2) = hydrogen_fraction
      mass_frac(3) = 0.0
      mass_frac(4) = 0.0
      mass_frac(5) = he3_fraction
      mass_frac(6) = helium_fraction
      mass_frac(7) = c12_fraction
      mass_frac(8) = c13_fraction
      mass_frac(9) = n14_fraction
      mass_frac(10) = o16_fraction
      mass_frac(11) = o18_fraction
      mass_frac(12) = 0.0
      mass_frac(13) = 0.0
! *******************************************************************
! BEGIN CALCULATION OF SCREENING CORRECTION.
! *******************************************************************
!  THE BASIC REFERENCES ARE SALPETER, AUSTRALIAN JOURNAL OF PHYSICS,
!  VOL. 7, 373 (1954). THE FORMULA FOR WEAK SCREENING THAT IS BEING
!  PROGRAMMED IS EQUATION (25) OF THIS PAPER. THE OTHER IMPORTANT
!  REFERENCES ARE: DEWITT, GRABOSKE, AND COOPER, AP. J. 181, 439 (1973)
!  AND GRABOSKE ET AL., AP. J. 181, 457 (1973). THE VALUES OF EMU AND
!  ZET ARE ESSENTIAL FOR COMPUTING WEAK SCREENING; THE VALUE OF AMU IS
!  USED IN AN NON-ESSENTIAL WAY IN THIS COMPUTATION. XTR IS USED IN
!  COMPUTING INTERMEDIATE SCREENING.
! AMU IS ONE OVER THE MEAN MOLECULAR WEIGHT OF THE IONS, MU SUB I .
! EMU IS ONE OVER THE ELECTRON MEAN MOLECULAR WEIGHT, MU SUB E.
!  EMU IS USED HERE AS THE NAME FOR THE SECOND PART OF THE ZETA FUNCTION
!  IN THE SALPETER EXPRESSION FOR WEAK SCREENING.
! XTR IS USED LATER IN THE INTERMEDIATE SCREENING CALCULATION.  THE AVERAGE
!  OF THE QUANTITY Z**(3B -1) IS EQUAL TO XTR/AMU.
! ZET IS THE FIRST PART OF THE SALPETER SCREENING ZETA VARIABLE.
! MU = SUM OVER I OF [X(I)/A(I)].
! MU SUB E = SUM OVER I [ Z(I)*X(I)/A(I)].
      mu_ion_inv = 0.
      mu_e_inv = 0.
      xtr = 0.
      zeta0 = 0.
      do i = 1,num_isotopes
         trm = mass_frac(i)/atomic_mass(i)
         mu_ion_inv = mu_ion_inv+trm
         mu_e_inv = mu_e_inv+trm*atomic_charge(i)
         xtr = xtr+trm*atomic_charge(i)**1.58
         zeta0 = zeta0+trm*atomic_charge(i)**2
   10 continue
      end do
! DL AND DT ARE THE THE LOG10 OF THE DENSITY AND TEMPERATURE.
!  THE UNIT OF TEMPERATURE IS 10^9 K AND THE UNIT OF DENSITY IS
!  GM PER CM^3 .
! DD = LOG RHO TO THE BASE 10.
! CLN = LN10.  CLN IS CONVERSION BETWEEN LOG10 AND LN.
! CONVERT DENSITY TO UNLOGGED FORM.
! RWE = RHO/(MU SUB E), I. E., THE NUMBER OF ELECTRONS DIVIDED BY
!  AVOGADRO'S NUMBER.
      rho_over_mu_e = ( exp(ln10*log_density) )*mu_e_inv
! THE EXPRESSION FOR RWE WAS INCORRECT IN THE ORIGINAL YALE SUBROUTINE.
!  THE ORIGINAL VERSION HAD ( EXP(CLN*DL) ) DIVIDED BY EMU INSTEAD OF
!  MULTIPLIED BY EMU.  RWE IS USED LATER IN COMPUTING THE SCREENING
!  CORRECTION.
      log_rho_local = log_density
! SET RATES EQUAL TO ZERO FOR THE LOG_10(T) < TCUT(1)
      if(log_temperature.le.tcut(1)) then
         do i = 1,num_reactions
            rate(i) = 0.
   20    continue
         end do
         go to 200
      endif
! T9P13 IS THE TEMPERATURE IN UNITS OF 10^9 DEGREES K TO THE PLUS 1/3
!  POWER.  MINUS IS DENOTED BY M.  HERE T9 IS THE TEMPERATURE IN UNITS
!  OF 10^9 K, CONVERTED FROM THE LOG_10 (T) AND RHO IS THE DENSITY IN
!  CGS UNITS.
      density=exp(ln10*log_rho_local)
      t9 = exp(ln10*(log_temperature - 9.0d0))
      t9p13 = t9**cc13
      t9p23 = t9p13**2
      t9m13=1./t9p13
      t9m23=t9m13**2
      t9m1=1./t9
      t9m2=t9m1**2
      t9m12=1./dsqrt(t9)
      t9m32=t9m1*t9m12
! ***********************************
! F PRIME/F
! ***********************************
! THE NEXT PIECE OF CODE COMPUTES FIRST THE FERMI ENERGY DIVIDED BY KT, WHERE
!  PFMC2 IS THE FERMI MOMENTUM DIVIDED BY MASS OF ELECTRON TIMES C, ALL
!  SQUARED AND EFMKT IS THE FERMI ENERGY DIVIDED BY KT.  THE QUANTITY
!  FPRF IS THE RATIO OF F PRIME TO F IN SALPETER'S SCREENING CORRECTION.
!  THE VALUE OF FPRF IS DETERMINED IN THE INTERMEDIATE CASE BY AN
!  INTERPOLATION FORMULA DEPENDING UPON THE DEGREE OF DEGENERACY, AS
!  MEASURED BY DEGD = LOG10(E_F/KT). THE NUMERICAL VALUES FOR THE FIT
!  WERE TAKEN FROM SALPETER'S ORIGINAL PAPER, FIGURE 1. THE ONLY CHANGES
!  IN THIS PART OF THE SUBROUTINE WERE THE CORRECTION OF THE ERROR IN
!  THE DEFINITION OF RWE (SEE ABOVE) AND REFINEMENTS OF THE COEFFICIENTS
!  IN THE EXPRESSIONS FOR PFMC2 AND EFMKT.
      fermi_mom_sq=1.017677E-4*rho_over_mu_e**0.6666667
      fermi_energy_over_kt=5.92986*t9m1*(dsqrt(1.+fermi_mom_sq)-1.)
      if(fermi_energy_over_kt.le.1.E-2) then
         f_prime_over_f=1.0
      else
         log_degeneracy=dlog10(fermi_energy_over_kt)
         if(log_degeneracy.ge.1.5) then
            f_prime_over_f=0.0
         else
            f_prime_over_f=0.75793-0.54621*log_degeneracy-0.30964*log_degeneracy**2+0.12535*log_degeneracy**3+ &
            0.1203*log_degeneracy**4-0.012857*log_degeneracy**5-0.014768*log_degeneracy**6
         endif
      endif
! END OF CALCULATION OF FERMI ENERGY DIVIDED BY KT AND THE INTERPOLATION
!  FORMULA FOR (F PRIME/F), WHICH APPEARS IN SALPETER'S SCREENING
!  FORMULA.
! NOW WE GET TO THE COMPUTATION OF WEAK SCREENING (SEE ALSO PAGE 61 OF
!  NEUTRINO ASTROPHYSICS, WHICH GIVES ONLY A SIMPLIFIED FORMULA) AND
!  ALSO THE MORE COMPLICATED INTERMEDIATE AND STRONG SCREENING CASES
!  (SEE REFERENCES TO AP. J. 181 ABOVE FOR INTERMEDIATE
!  AND STRONG SCREENING), ESPECIALLY PAGE 465.
! XXL IS USED IN ALL THE SCREENING FORMULAE.  THIS QUANTITY IS THE
!  IS THE FUNCTION CALLED LAMBDA SUB ZERO BY GRABOSKE ET AL. THE
!  QUANTITY XXL**0.86 = XXL8 IS USED IN CALCULATING INTERMEDIATE
!  SCREENING.
! ZCURL IS THE QUANTITY DEFINED BY EQUATION (4) OF DEWITT ET AL; IT
!  IS THEIR Z WITH A CURLY SYMBOL ON ITS TOP. ZCURL IS USED IN WEAK
!  AND IN INTERMEDIATE SCREENING AND WAS FIRST DEFINED BY SALPETER.
!  ZCURL IS THE SAME AS SALPETER'S ZETA EXCEPT FOR THE FACTOR OF 1/AMU.
! (Z SUB 1 TIMES Z SUB 2)*XXL*ZCURL GIVES THE WEAK SCREENING FACTOR,
!  THE SAME AS SALPETER OR AS EQUATION (19) OF GRABOSKE ET AL.
! Z BAR IS THE SAME AS Z BAR OF DE WITT ET AL AND IS THE AVERAGE CHARGE
!  OF THE IONS.  IT IS EQUAL TO EMU/AMU.  THE QUANTITY (Z BAR)**0.28 =
!  Z28 OCCURS IN THE COMPUTATION OF INTERMEDIATE SCREENING.
! XXL6 IS USED FOR COMPUTING STRONG STRONG SCREENING.
! THE NOTATION USED HERE IS EXPLAINED IN LARGE PART BY
!  EQUATION (4) OF DEWITT ET AL., AP. J. 181, PAGE 439.
! THE FINAL EXPRESSION FOR WEAK SCREENING IS EXACTLY EQUAL TO SALPETER'S
!  FORMULA, WHICH INCLUDES A DEGENERACY CORRECTION. THE MORE GENERAL
!  EXPRESSIONS ARE GIVEN IN TABLE 4 AND EQUATION (19) OF GRABOSKE ET AL.
      lambda0=5.9426E-6*t9m32*dsqrt(density*mu_ion_inv)
      lambda0_23=lambda0**0.666667
      lambda0_86=lambda0**0.86
      z_curl=dsqrt((zeta0+f_prime_over_f*mu_e_inv)/mu_ion_inv)
      z_bar=mu_e_inv/mu_ion_inv
      z_curl_58=z_curl**0.58
      z_bar_28=z_bar**0.28
      z_bar_13=z_bar**cc13
      lambda0_zcurl=lambda0*z_curl
! COMPUTE SCREENING FOR EACH OF THE REACTIONS.
      do i=1,num_reactions
         weak_screening_u=lambda0_zcurl*charge_product(i)
         if(weak_screening_u.le.weak_screening_threshold) then
! WEAKSCREENING IS A NUMERICAL PARAMETER PASSED IN THE FLUX COMMON
!  BLOCK. TO OBTAIN THE GRABOSKE ET AL. AND SALPETER STANDARD RESULTS,
!  USE: WEAKSCREENING = 0.03.  FOR THE STANDARD SOLAR MODEL, THIS IS THE
!  VALUE THAT SHOULD BE ADOPTED. TO INVESTIGATE THE EFFECT OF ALWAYS USING
!  WEAK SCREENING, USE A LARGE VALUE FOR WEAKSCREENING, E. G., 30.  AS
!  LONG AS WEAKSCREENING IS ASSUMED TO BE BIGGER THAN ONE, THE PROGRAM
!  WILL ALWAYS CALCULATE FOR THE SUN WITH THE WEAK SCREENING
!  APPROXIMATION.
! UTOT IS THE FINAL SCREENING CORRECTION WHICH APPEARS IN THE
!  RATE EXPRESSION AS: EXP(UTOT) .
            screening_factor(i)=weak_screening_u
         else
            intermediate_screening_u=0.38*lambda0_86*xtr*z86(i)/(mu_ion_inv*z_curl_58*z_bar_28)
            if(weak_screening_u.le.2.) then
               screening_factor(i)=intermediate_screening_u
            else
               strong_screening_u=0.624*z_bar_13*lambda0_23*(z53(i)+0.316*z_bar_13*z43(i)+0.737* &
               z23(i)/(z_bar*lambda0_23))
               if(strong_screening_u.lt.intermediate_screening_u.or.weak_screening_u.ge.5.) then
                  screening_factor(i)=strong_screening_u
               else
                  screening_factor(i)=intermediate_screening_u
               endif
            endif
         endif
   30 continue
      end do
! ****************************************************************
! END OF SCREENING CALCULATION. WEAK AND INTERMEDIATE SCREENING FORMS
!  ARE GIVEN CORRECTLY.  STRONG SCREENING WAS NOT CHECKED BECAUSE IT IS
!  NOT RELEVANT FOR THE SUN.
! ****************************************************************
      nz=1
      if(hydrogen_fraction.eq.0.0) then
         be7_electron_frac=0.
         c12_alpha_frac=0.
         goto 50
      endif
      nz=8
! **************************************************************
!  CALCULATE REACTION RATES FOR THE THREE PRINCIPAL REACTIONS OF
!   THE PP CHAIN: PP, HE3+HE3, HE3 +HE4, AND THE FOUR PROTON
!   BURNING REACTIONS ON C12: C13, N14, AND O16.
! **************************************************************
! R1 IS (T SUB 9)^(-3/2) TIMES (S SUB EFF)/(S SUB 0). THE CORRECT
!  EXPRESSION FOR (S SUB EFF)/(S SUB 0) IS GIVEN IN EQUATION 3.14 OF
!  NEUTRINO ASTROPHYSICS.  THE NUMERICAL FORM THAT IS USED IS EQUATION
!  52 OF FOWLER, CAUHLAN, AND ZIMMERMAN, VOL. 5, 1967.
! RATE(I) IS THE RATE OF THE DIFFERENT REACTIONS PER SECOND PER GRAM,
!  EXCEPT THAT THE MASS FRACTIONS ARE OMITTED.
! THE PREVIOUS YALE VERSION HAD THE RATE FOR THE O16 + P REACTION
!  MULTIPLIED BY T9**(-1/7). THIS FACTOR IS INCORRECT AND HAS BEEN
!  REMOVED; IT APPEARED BEFORE AS AN IF STATEMENT REFERRING ONLY TO
!  RATE(7).
      do i=1,7
!         R1=T9M23+Q1(I)*T9M13+Q2(I)+Q3(I)*T9P13+Q4(I)*T9P23+Q5(I)*T9
! MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
! RATE DERIVATIVES
         r1=t9m23+q1(i)*t9m13+qs0e_scale(i)*(q2(i)+q3(i)*t9p13)+ &
            qqs0ee_scale(i)*(q4(i)*t9p23+q5(i)*t9)
         rate(i)=density*r1*exp(q6(i)*t9m13+q7(i)+(q8(i)*t9)**2+screening_factor(i))
         rate(i) = rate(i)*cross_section_scale(i)
         if(rate(i).lt.1.E-30) rate(i)=0.0d0
   40 continue
      end do
! ***************************************************************
! END OF CALCULATION OF REACTION RATES FOR FIRST 7 REACTIONS.
! ***************************************************************
! ***********************************************
! BE7 BURNING
! ***********************************************
! CALCULATE THE BURNING OF BE7 BY PROTONS (WHICH PRODUCES THE MOST
!  EXPERIMENTALLY ACCESSIBLE SOLAR NEUTRINOS) AND THE BURNING OF BE7
!  BY ELECTRON CAPTURE (THE DOMINANT PROCESS). THE PREVIOUSLY USED
!  EXPRESSION CONTAINED A TERM THAT WAS PHYSICALLY INCORRECT
!  AND ALSO USED CROSS SECTION FACTORS THAT WERE OUT OF DATE.
! THE ELECTRON CAPTURE RATE IN SEC^(-1) IS GIVEN BY EQUATION (3.18) OF
!  NEUTRINO ASTROPHYSICS.  CAN OMIT THE FACTOR OF RHO (IN CGS UNITS)
!  WHICH ALSO APPEARS IN THE BE7PROTON RATE. THE BE7 PROTON CAPTURE
!  RATE IS GIVEN BY TABLE 3.2 AND EQUATION (3.12).  NOTE THAT 1 OVER MU
!  SUB E IS EQUAL TO EMU.
! HERE WE USE THE NOTATION BE7 + PROTON = BE7PROTON AND
!  BE7 + E = BE7ELECTRON.
! F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY ELECTRON CAPTURE.
! 1-F1 IS THE FRACTION OF THE BE7 THAT IS BURNED BY PROTON CAPTURE.
! F3 IS THE FRACTION OF THE N14 THAT IS BURNED BY P, ALPHA REACTION.
! 1-F3 IS THE FRACTION OF THE N15 THAT IS BURNED BY P, GAMMA REACTION.
! SEE TABLE 21 OF BAHCALL AND ULRICH (1988), REV. MOD. PHYS. 60.
! 10/13/97. I changed Temp3 (i.e., tau) by 5/10^5 as a result of
! using bare nuclear masses. Previously coefficient was -10.26202.
! I did not change Be7electron for the slightly different S0 for
! electron capture, since that is done in SStandard. Previous coefficient
! in Be7electron expression was (3.126571E+5). 10/14/97.
!
         be7_electron_rate = (1.752E-10)*t9m12*(1.0 + 0.004*(1000.*t9-16.))
         be7_electron_rate = be7_electron_rate*mu_e_inv*cross_section_scale(15)
         be7_temp_factor = (-10.2625*t9m13)
         be7_proton_rate = (3.128813E+5)*hydrogen_fraction*cross_section_scale(16)*exp(be7_temp_factor)

! INCLUDE FOR BE7PROTON THE T9M23 FACTOR AND ALL CORRECTIONS PROPORTIONAL TO
!  Q1,...,Q5 FROM EQUATION 3.14 OF NEUTRINO ASTROPHYSICS. THESE
!  CORRECTIONS ARE DEFINED EARLIER IN THIS SUBROUTINE.
!         QRBE7 = T9M23 + Q1(8)*T9M13 + Q2(8)+ Q3(8)*T9P13
!     $          + Q4(8)*T9P23 + Q5(8)*T9
! MHP 9/14 ADDED THE ABILITY TO ALTER DERIVATIVES INDEPENDENTLY
         be7_q_factor = t9m23 + q1(8)*t9m13 + qs0e_scale(8)*(q2(8)+ q3(8)*t9p13) &
                + qqs0ee_scale(8)*(q4(8)*t9p23 + q5(8)*t9)
         be7_proton_rate = be7_proton_rate*be7_q_factor
! CALCULATE THE SCREENING CORRECTION FOR BE7 + P REACTION.  USE WEAK AND
!  INTERMEDIATE SCREENING FORMULAE.
         be7p_charge_product = 4.0
         be7p_z86 = 5.7790
         weak_screening_u = lambda0_zcurl*be7p_charge_product
         if(weak_screening_u.le.weak_screening_threshold) then
            be7p_screening_u = weak_screening_u
         else
            intermediate_screening_u = 0.38*lambda0_86*xtr*be7p_z86/(mu_ion_inv*z_curl_58*z_bar_28)
            be7p_screening_u = intermediate_screening_u
          endif
         be7_proton_rate = be7_proton_rate*exp(be7p_screening_u)
! END OF CALCULATION OF SCREENING CORRECTION FOR BE7 + P REACTION.
! MULTIPLY RATES BY FACTOR OF RHO/[(ATOMIC MASS UNIT)*A(BE7)] TO GET
!  IN UNITS OF GM^{-1}.  CALL FACTOR CAMUBE7. Corrected for Be7 bare
! mass number. Previously used 8.582295E+22 in the expression below.
! These corrections were made on 10/14/97.
!
         be7_mass_factor = density*8.584981E+22
         be7_proton_rate = be7_mass_factor*be7_proton_rate
         be7_electron_rate = be7_mass_factor*be7_electron_rate
! END OF MULTIPLICATION INSERTED NOVEMBER 6, 1990.
         be7_electron_frac = be7_electron_rate/(be7_electron_rate + be7_proton_rate)
! *********************************************************************
! END OF CALCULATION OF CRUCIAL BE7 ELECTRON CAPTURE AND PROTON CAPTURE
!  RATES AND THEIR RATIO.
! *********************************************************************
! ***************************
! N15 + P BRANCHING.
! ***************************
! THE FOLLOWING STATEMENTS COMPUTE THE BRANCHING OF N15(P,ALPHA)C12
!  AND N15(P,GAMMA)O16. THESE STATEMENTS REPLACE OUTDATED STATEMENTS
!  IN THE YALE CODE. THE CNO CROSS-SECTION FACTORS ARE FROM TABLE 3.4
!  OF NEUTRINO ASTROPHYSICS.  THE RATIO OF THE REACTIONS DEPENDS ONLY
!  UPON THE EFFECTIVE ZERO ENERGY S-FACTOR, WHICH IS S0(ZERO ENERGY)
!  TIMES THE COMBINATION OF TEMPERATURE AND S-FACTOR DERIVATIVES/S0
!  THAT WAS USED PREVIOUSLY AS R1 IN THE RATE CALCULATIONS. THE
!  NUMERICAL COEFFICIENTS THAT APPEAR IN THE RATE WERE REPRESENTED
!  BY THE Q1(J),...Q5(J) FOR THE OTHER REACTIONS.
!  THE QVALUES FOR THE N15 REACTIONS HAVE BEEN COMPUTED SEPARATELY.
! DBG 8/94 APPLIED MHP UPDATE TO NUCLEAR REACTIONS
! To agree with Solar Fusion Workshop paper, the value of S0 in keV-b
! has been changed from 78000. to 675000. on 9/25/97. On 10/14/97,
! JNB changed cues.f so as to compute the Q-coefficients for the N15 + p
! reactions.  Also, checked that the coefficients are the same as the
! ones given earlier in energy.f when we use the older CNO data.
!
      o16_gamma_rate = t9m23 + 0.0273016*t9m13 + 0.14374 + 0.027490*t9p13 &
                + 6.14685*t9p23 + 2.98940*t9
! MULTIPLY BY THE VALUE OF S0 IN KEV-B.
!      O16GAMMA = O16GAMMA*64.
! MHP 9/14 ADDED THE OPTION TO MODIFY THE RELATIVE CROSS SECTIONS
! FOR N15+P -> C12+ALPHA AND O16+GAMMA
      o16_gamma_rate = o16_gamma_rate*64*o16_gamma_scale
!
      c12_alpha_n15p_rate = t9m23 + 0.0273016*t9m13 + 2.01186 + 0.384763*t9p13 &
                + 17.0579*t9p23 + 8.29580*t9
!      C12ALPHA = C12ALPHA*67500.
      c12_alpha_n15p_rate = c12_alpha_n15p_rate*67500*c12_alpha_scale
      c12_alpha_frac = c12_alpha_n15p_rate/(c12_alpha_n15p_rate + o16_gamma_rate)
      o16_gamma_frac = 1.0d0 - c12_alpha_frac
! END OF NEW ROUTINE FOR THE BRANCHING OF N15 + P .
   50 do i=nz,num_reactions
         rate(i)=0.
   60 continue
   end do
!***MHP 3/91 ALPHA CAPTURE REACTIONS UPDATED TO CAUGHLAN AND FOWLER(1988)
!   RATES.  THE RATES ARE EXPRESSED IN THE SAME TERMS USED BY CZ, WITH
!   THE CONVERSION FACTOR IN THE FRONT OBTAINED FROM VANDENBERG'S
!   NOTES ON THE REACTION RATES.
!  RATE(8)  HE4+C13
!  RATE(10) HE4+C12=>O16
!  RATE(11) HE4+N14=>O18
!  RATE(12) TRIPLE ALPHA
      if (.not. (log_temperature.lt.tcut(4))) then
! C13(ALPHA,N) O16
      r1=t9m23+0.0129d0*t9m13+2.04d0+0.184d0*t9p13
      a1 = 6.77d15*exp(-32.329d0*t9m13-(t9/1.284d0)**2)
      a2 = 3.82d5*exp(-9.373*t9m1)
      a3 = 1.41d6*exp(-11.873*t9m1)
      a4 = 2.00d9*exp(-20.409*t9m1)
      a5 = 2.92d9*exp(-29.283*t9m1)
      rate(8) = 1.157126d22*density*exp(screening_factor(8))*(a1*r1+t9m32* &
           (a2+a3+a4+a5))
! C12(ALPHA,GAMMA)O16
      r1 = 1.0d0/(1.0d0+0.0489d0*t9m23)
      r2 = 1.0d0/(1.0d0+0.2654d0*t9m23)
      a1 = t9m2*exp(-32.120*t9m13)
      a2 = 1.04d8*r1**2*exp(-(t9/3.496)**2)
      a3 = 1.76d8*r2**2
      a4 = 1.25d3*t9m32*exp(-27.499*t9m1)
      a5 = 1.43d-2*t9**5*exp(-15.541*t9m1)
      rate(10) = 1.25388d22*density*exp(screening_factor(10))*(a1*(a2+a3)+a4+a5)
! N14(ALPHA,GAMMA)F18 + F18=>O18+EPLUS+NU
      r1 = t9m23+0.012d0*t9m13+1.45d0+0.177d0*t9p13+1.97d0*t9p23 &
           +0.406d0*t9
      a1 = 7.78d9*exp(-36.031d0*t9m13-(t9/0.881d0)**2)
      a2 = t9m32*2.36d-10*exp(-2.798d0*t9m1)
      a3 = t9m32*2.03d0*exp(-5.054d0*t9m1)
      a4 = t9m23*1.15d4*exp(-12.310*t9m1)
      rate(11)= 1.07452d22*density*exp(screening_factor(11))*(a1*r1+a2+a3+a4)
! TRIPLE ALPHA
      rate(12) = 1.565315d21*density**2*t9m1*t9m2*2.79E-8* &
                 exp(-4.4027*t9m1+screening_factor(12))
      end if
  100 continue
      rate(9) = 0.0d0
      rate(13) = 0.0d0
! END OF XEROING OUT OF REACTIONS 9 AND 13.
      do i=1,num_reactions
         if(rate(i).le.1.E-5) rate(i) = 0.0
  130 continue
      end do
! ******************************************************
! RATES PER 10^9 YEARS PER ATOMIC MASS UNIT: HRK(IU)
! ******************************************************
! HR1, ..., HR13 ARE THE RATES OF THE INDIVIDUAL REACTIONS.
!  THE INTERPRETATION OF WHICH REACTION GOES WITH WHICH SYMBOL CAN BE
!  MADE EASILY BY LOOKING AT THE DEFINITIONS OF THE EG(I)'S.
!  THE ABUNDANCES ARE UPDATED IN SUBROUTINE KEMCOM USING THESE MATRICES.
! C21 IS THE PRODUCT OF (10^9 YEARS/1 SECOND)*(1 ATOMIC MASS UNIT/1
!  GRAM). I HAVE USED HERE SIDEREAL YEAR IN CONVERTING TO SECONDS.
  200 rate_pp(zone_idx)=rate(1)*c21
      rate_he3_he3(zone_idx)=rate(2)*c21
      rate_he3_he4(zone_idx)=rate(3)*c21
      rate_c12_p(zone_idx)=rate(4)*c21
      rate_c13_p(zone_idx)=rate(5)*c21
      rate_n14_p(zone_idx)=rate(6)*c21
      rate_o16_p(zone_idx)=rate(7)*c21
      rate_c13_alpha(zone_idx)=rate(8)*c21
      rate_zero9(zone_idx)=rate(9)*c21
      rate_c12_alpha(zone_idx)=rate(10)*c21
      rate_n14_alpha(zone_idx)=rate(11)*c21
      rate_triple_alpha(zone_idx)=rate(12)*c21
      rate_zero13(zone_idx)=rate(13)*c21
      frac_c12_alpha(zone_idx)=c12_alpha_frac
      frac_be7_electron(zone_idx)=be7_electron_frac
! ****************************************
! END OF COMPUTATION OF HRK(IU).
! ****************************************
      return
end subroutine rates

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! KC 2025-05-31 SAFEDIVEXP
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
subroutine safedivexp(numerator, exponent)
      double precision :: numerator, exponent

      if (exponent .lt. 709.7827d0) then
         numerator = numerator/exp(exponent)
      else
         numerator = 0.0d0
      endif

      return
end subroutine safedivexp

!----------------------------------------------------------------------
! ifermi12
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ifermi12.f; only variable names, source form, and comment
! style were updated.
!
! This routine applies a rational function expansion to get the
! inverse Fermi-Dirac integral of order 1/2 (i.e. the electron
! degeneracy parameter eta such that F_{1/2}(eta) = fermi_half_integral)
! when it is equal to fermi_half_integral.
! Maximum error is 4.19d-9. Reference: Antia, ApJS 84,101 1993.
double precision function ifermi12(fermi_half_integral)

      implicit none

! fermi_half_integral: value of the order-1/2 Fermi-Dirac integral
!                       whose inverse (the degeneracy parameter) is
!                       sought.
      double precision, intent(in) :: fermi_half_integral

!..declare
      integer :: term_idx, deg_num_small, deg_den_small, deg_num_large, &
           deg_den_large
      double precision :: fd_order, coef_num_small(12), coef_den_small(12), &
           coef_num_large(12), coef_den_large(12), numerator, denominator, &
           scaled_arg
!     unused in the original: z,drn
!..load the coefficients of the expansion
      data  fd_order,deg_num_small,deg_den_small,deg_num_large,deg_den_large &
           /0.5d0, 4, 3, 6, 5/
      data  (coef_num_small(term_idx),term_idx=1,5)/ 1.999266880833d4,   5.702479099336d3, &
           6.610132843877d2,   3.818838129486d1, &
           1.0d0/
      data  (coef_den_small(term_idx),term_idx=1,4)/ 1.771804140488d4,  -2.014785161019d3, &
           9.130355392717d1,  -1.670718177489d0/
      data  (coef_num_large(term_idx),term_idx=1,7)/-1.277060388085d-2,  7.187946804945d-2, &
           -4.262314235106d-1,  4.997559426872d-1, &
           -1.285579118012d0,  -3.930805454272d-1, &
           1.0d0/
      data  (coef_den_large(term_idx),term_idx=1,6)/-9.745794806288d-3,  5.485432756838d-2, &
           -3.299466243260d-1,  4.077841975923d-1, &
           -1.145531476975d0,  -6.067091689181d-2/

      if (fermi_half_integral .lt. 4.0d0) then
       numerator  = fermi_half_integral + coef_num_small(deg_num_small)
       do term_idx=deg_num_small-1,1,-1
        numerator  = numerator*fermi_half_integral + coef_num_small(term_idx)
       enddo
       denominator = coef_den_small(deg_den_small+1)
       do term_idx=deg_den_small,1,-1
        denominator = denominator*fermi_half_integral + coef_den_small(term_idx)
       enddo
       ifermi12 = log(fermi_half_integral * numerator/denominator)

      else
       scaled_arg = 1.0d0/fermi_half_integral**(1.0d0/(1.0d0 + fd_order))
       numerator = scaled_arg + coef_num_large(deg_num_large)
       do term_idx=deg_num_large-1,1,-1
        numerator = numerator*scaled_arg + coef_num_large(term_idx)
       enddo
       denominator = coef_den_large(deg_den_large+1)
       do term_idx=deg_den_large,1,-1
        denominator = denominator*scaled_arg + coef_den_large(term_idx)
       enddo
       ifermi12 = numerator/(denominator*scaled_arg)
      end if
      return
end function ifermi12

!----------------------------------------------------------------------
! zfermim12
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original zfermim12.f; only variable names, source form, and comment
! style were updated.
!
!..this routine applies a rational function expansion to get the
!..fermi-dirac integral of order -1/2 evaluated at degeneracy_parameter.
!..maximum error is 1.23d-12. reference: antia apjs 84,101 1993.
double precision function zfermim12(degeneracy_parameter)

      implicit none

! degeneracy_parameter: electron degeneracy parameter eta at which the
!                        order -1/2 Fermi-Dirac integral is evaluated.
      double precision, intent(in) :: degeneracy_parameter

!..declare
      integer :: term_idx, deg_num_small, deg_den_small, deg_num_large, &
           deg_den_large
!       double precision fd_order,coef_num_small(12),coef_den_small(12),coef_num_large(12),coef_den_large(12),numerator,denominator,xx  ! KC 2025-05-31
      double precision :: coef_num_small(12), coef_den_small(12), &
           coef_num_large(12), coef_den_large(12), numerator, denominator, &
           scaled_arg
!..load the coefficients of the expansion
!       data  fd_order,deg_num_small,deg_den_small,deg_num_large,deg_den_large /-0.5d0, 7, 7, 11, 11/  ! KC 2025-05-31
      data  deg_num_small,deg_den_small,deg_num_large,deg_den_large /7, 7, 11, 11/
      data  (coef_num_small(term_idx),term_idx=1,8)/ 1.71446374704454d7,    3.88148302324068d7, &
           3.16743385304962d7,    1.14587609192151d7, &
           1.83696370756153d6,    1.14980998186874d5, &
           1.98276889924768d3,    1.0d0/
      data  (coef_den_small(term_idx),term_idx=1,8)/ 9.67282587452899d6,    2.87386436731785d7, &
           3.26070130734158d7,    1.77657027846367d7, &
           4.81648022267831d6,    6.13709569333207d5, &
           3.13595854332114d4,    4.35061725080755d2/
      data (coef_num_large(term_idx),term_idx=1,12)/-4.46620341924942d-15, -1.58654991146236d-12, &
           -4.44467627042232d-10, -6.84738791621745d-8, &
           -6.64932238528105d-6,  -3.69976170193942d-4, &
           -1.12295393687006d-2,  -1.60926102124442d-1, &
           -8.52408612877447d-1,  -7.45519953763928d-1, &
           2.98435207466372d0,    1.0d0/
      data (coef_den_large(term_idx),term_idx=1,12)/-2.23310170962369d-15, -7.94193282071464d-13, &
           -2.22564376956228d-10, -3.43299431079845d-8, &
           -3.33919612678907d-6,  -1.86432212187088d-4, &
           -5.69764436880529d-3,  -8.34904593067194d-2, &
           -4.78770844009440d-1,  -4.99759250374148d-1, &
           1.86795964993052d0,    4.16485970495288d-1/

      if (degeneracy_parameter .lt. 2.0d0) then
       scaled_arg = exp(degeneracy_parameter)
       numerator = scaled_arg + coef_num_small(deg_num_small)
       do term_idx=deg_num_small-1,1,-1
        numerator = numerator*scaled_arg + coef_num_small(term_idx)
       enddo
       denominator = coef_den_small(deg_den_small+1)
       do term_idx=deg_den_small,1,-1
        denominator = denominator*scaled_arg + coef_den_small(term_idx)
       enddo
       zfermim12 = scaled_arg * numerator/denominator
!..
      else
       scaled_arg = 1.0d0/(degeneracy_parameter*degeneracy_parameter)
       numerator = scaled_arg + coef_num_large(deg_num_large)
       do term_idx=deg_num_large-1,1,-1
        numerator = numerator*scaled_arg + coef_num_large(term_idx)
       enddo
       denominator = coef_den_large(deg_den_large+1)
       do term_idx=deg_den_large,1,-1
        denominator = denominator*scaled_arg + coef_den_large(term_idx)
       enddo
       zfermim12 = sqrt(degeneracy_parameter)*numerator/denominator
      end if
      return
end function zfermim12

end module net_lib
