!----------------------------------------------------------------------
! eqrelv
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqrelv.f; only variable names, source form, and comment
! style were updated.
!
! Fully-ionized-gas equation of state ("relativistic" electron gas via
! the tabulated Fermi-Dirac integrals in common/ccr/, following
! Prather's thesis). Given log10(T), log10(P), and the mean molecular
! weights, iterates on density until the electron+ion+radiation
! pressure sum matches the given P, then returns density and (if
! in_atmosphere is false, via eqstat2.f90's call) thermodynamic
! derivatives. Called from eqstat2.f90 both for the T >= saha cutoff
! branch (fully ionized gas assumed) and to interpolate against the
! eqsaha/eqscve result near the ionization cutoff.
!
!          QDTT,QDTP = NAT-LOG DERIVATIVES OF QDT=D(LOG D )/D(LOG T )
!          QAT,QAP = NAT-LOG DERIVATIVES OF THE ADIABATIC DERIVATIVE
! DBG 7/92 MODIFIED TO INCLUDE DEBYE-HUCKEL CORRECTION.
!     CORRECTION IS APPLIED TO PRESSURE AND ALL THERMODYNAMIC VARIABLES.
!     TREATMENT IS ALSO VALID FOR ELECTRON DEGENERACY.  A RAMP FUNCTION
!     BETWEEN THE NONDEGENERATE AND DEGERNATE ELECTRON GAS CASES OF THE
!     D.H. CORRECTION IS USED (THIS IS NOT CORRECT BUT, HOPEFULLY IS
!     CLOSE ENOUGH.  CONSTANT TERM CDH IS CALCULATED IN SETUPS.
subroutine eqrelv(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, beta, ion_mean_weight_inverse, &
     electron_mean_weight_inverse, electron_degeneracy_parameter, &
     dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, &
     dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt, &
     adiabatic_gradient_dp, specific_heat_cp_dt, specific_heat_cp_dp)

      use const_lib
      use luout_lib
      implicit none

      double precision, intent(inout) :: log10_temperature
      double precision, intent(in) :: temperature
      double precision, intent(inout) :: log10_pressure
      double precision, intent(in) :: pressure
      double precision, intent(inout) :: log10_density
      double precision, intent(out) :: density
      double precision, intent(in) :: beta, ion_mean_weight_inverse, &
           electron_mean_weight_inverse
      double precision, intent(out) :: electron_degeneracy_parameter, &
           dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, &
           dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt, &
           adiabatic_gradient_dp, specific_heat_cp_dt, specific_heat_cp_dp

! common/ccr/: this batch's own block (no reuse precedent found
! elsewhere in the already-converted sources). fermi_table_x_grid is
! the tabulated X(=log D - log mu_e - 3/2 log T) grid; fermi_table_eta
! is the tabulated electron degeneracy parameter (eta) at those grid
! points; fermi_table_data(quantity,x_index,t_index) is the main
! table of Fermi-Dirac-integral-derived quantities (electron pressure,
! internal energy, and their T-derivatives; see the kk=1..5 loop
! below); fermi_table_x_lookup is a coarse index used to jump straight
! to the right 3-point interpolation window in fermi_table_x_grid for
! a given X.
      double precision :: fermi_table_x_grid(43), fermi_table_eta(43), &
           fermi_table_data(5,43,20)
      integer :: fermi_table_x_lookup(261)
      common/ccr/ fermi_table_x_grid, fermi_table_eta, fermi_table_data, &
           fermi_table_x_lookup



! common/const3/: not used in this file.
      double precision :: cdelrl, cmixl, cmixl2, cmixl3, clndp, &
           seconds_per_year
      common/const3/ cdelrl, cmixl, cmixl2, cmixl3, clndp, seconds_per_year

! DBG 7/92 CDH IS CONSTANT TERM DEFINED IN SETUPS
!     RAMP FUNCTION BETWEEN NO ELECTRON DEGENERACY ETA .LT. ETADH0
!     AND FULL ELECTRON DEGENERACY ETA .GT. ETADH1 VERSION OF D.H.
!     CORRECTION.
!     IF LDH = .TRUE. THEN APPLY D.H. CORRECTION.
!     ZDH IS ARRAY OF RELATIVE MASS FRACTIONS OF LAOL METAL MIXTURE
!     SUMMED TO 1.0.
! DBG 7/92 COMMON BLOCK ADDED TO COMPUTE DEBYE-HUCKEL CORRECTION.
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh

      save

! --- locals ---
! Scratch algebra kept close to the original short names (rather than
! invented descriptive names) given the density of unlabeled physics
! here; see the inline comments (carried over from the original) for
! what each quantity represents. qpd (d ln P/d ln rho, a byproduct of
! the Newton iteration) is intentionally a separate local from the
! dlnrho_dlnp output argument (d ln rho/d ln P, its reciprocal) --
! they are related but distinct quantities in the original code (QPD
! vs QDP) and must not be aliased to the same name.
      double precision :: pinv, beta1, tl8, pl8, pr, ramu, ramut, emul, &
           emul8, ttud, ttcu, ttdcu, dml, almix, cmfdh0, cmfdh1, cmfdh, &
           xx, pr6, pr9, pr8, pr7, cden1, cden2, cden3, dpel2, dx1, dx2, &
           dx3, cl1, cl2, cl3, dcl1, dcl2, dcl3, t0, t1, slope, pel, &
           dpel1, dl8, d8, pa, pe, pdh, pt, pdhp, ptl, factor, pr8l, &
           pr9l, factl, fact, corr, peq, qped, qqpedd, pr1, pr2, pr3, &
           pr4, pr5, qpeqd, qpeqt, qpt, qpd, qqptt, qqpdt, qqpdd, qdpp, &
           ue, ueq, queqt, pd, pdt, ua, ur, udh, u, udhu, ui, pdu, qutd, &
           qudt, qutp, qqutt, qqudt, acp, bcp, qcpi
      integer :: jt1, jt2, jt3, nden, id1p, ii, kk, id1, id2, id3, nn, id
      double precision :: df1(5,3), df2(5,3), ff(5,3), ffe(5)

!  THE ELECTRON PRESSURE (TABULATED IN FERMI TABLE) FOR A GIVEN (RHO,T)
!  IS FOUND HERE; THE DENSITY IS ADJUSTED SO THAT THE RESULTING TOTAL
!  PRESSURE (PE + PA + PR) EQUALS THE INITIAL P.  FOR DETAILS
!  CONSULT PRATHER'S THESIS.
!  SUFFIX E=ELECTRON;R=RADIATION;A=NUCLEON;I.E.PA=PRESSURE DUE TO NUCLEI
      pinv = 1.0d0/pressure
      beta1 = 1.0d0 - beta
      tl8 = log10_temperature
      pl8 = log10_pressure
      pr = beta1*pressure
      ramu = gas_constant*ion_mean_weight_inverse
      ramut = ramu*temperature
      emul = dlog10(electron_mean_weight_inverse)
      emul8 = emul
!  PE IS TABULATED AS A FUNCTION OF X(=LOG D - LOG MUE - 3/2LOG T) AND
!  Y(= LOG T). FIND INDEX(JT1,JT2,JT3) FOR 3 PT INTERPOLATION IN Y.
! MHP 10/02 CORRECTED MIXED TYPE
      jt1 = int(5.0d0*(log10_temperature - 5.0d0))
      jt1 = max0(1,min0(jt1,18))  ! Force jt1 to be in the table
      jt2 = jt1 + 1
      jt3 = jt1 + 2
      ttud = 5.0d0*(tl8 - (5.0d0 + 0.20d0*dfloat(jt1)))
      ttcu = 0.50d0*ttud*(ttud - 1.0d0)
      ttdcu = ttud - 0.50d0
      dml = log10_density + emul
      nden = 1
      id1p = -100
! DBG 7/92 CALCULATE D.H. TERMS: ALPHA-MIXTURE, MASS-FRAC0, MASS-FRAC1
      if (ldh) then
         almix = 0.0d0
         do ii=1,18
            almix = almix + dhnue(ii)*zdh(ii)
         end do
         cmfdh0 = 2.0d0*xxdy+1.5d0*yydh+0.5d0*(zzdh+almix)
         cmfdh1 = xxdy+yydh+0.5d0*almix
      end if
   10 continue
! BEGIN ITERATION LOOP FOR CORRECT DENSITY
!  FIND INDEX (ID1,ID2,ID3) FOR 3-PT INTERPOLATION IN X
      xx = dml - 1.50d0*tl8
      pr6= 20.0d0*(xx - fermi_table_x_grid(1)) + 1.0d0
! MHP 10/02 CORRECTED MIXED TYPE
      kk = int(pr6)
!      KK = PR6
      kk = min0(261,max0(1,kk))
      id1 = fermi_table_x_lookup(kk)
!  IF INDEX UNCHANGED FROM PREVIOUS LOOP,SKIP THIS SECTION
      if (id1.eq.id1p)  go to 30
      id1p = id1
      id2 = id1 + 1
      id3 = id1 + 2
      pr9 = 1.0d0/( fermi_table_x_grid(id3) - fermi_table_x_grid(id2) )
      pr8 = 1.0d0/( fermi_table_x_grid(id2) - fermi_table_x_grid(id1) )
      pr7 = 1.0d0/( fermi_table_x_grid(id1) - fermi_table_x_grid(id3) )
      cden1 = -pr7*pr8
      cden2 = -pr8*pr9
      cden3 = -pr9*pr7
!  INTERPOLATION IN Y IS NEWTONIAN(EQUAL SPACING IN Y)
      do 20 nn = 1,3
       id = id1 + nn - 1
       df1(1,nn) = fermi_table_data(1,id,jt2) - fermi_table_data(1,id,jt1)
       df2(1,nn) = fermi_table_data(1,id,jt3)-fermi_table_data(1,id,jt2)-df1(1,nn)
       ff(1,nn) = fermi_table_data(1,id,jt1) + ttud*df1(1,nn) + ttcu*df2(1,nn)
   20 continue
      dpel2 = 2.0d0*(ff(1,1)*cden1 + ff(1,2)*cden2 + ff(1,3)*cden3)
   30 dx1 = xx - fermi_table_x_grid(id1)
      dx2 = xx - fermi_table_x_grid(id2)
      dx3 = xx - fermi_table_x_grid(id3)
!  INTERPOLATION IN X IS LAGRANGIAN(UNEQUAL TABLE SPACING IN X)
      cl1 = dx2*dx3*cden1
      cl2 = dx3*dx1*cden2
      cl3 = dx1*dx2*cden3
      dcl1 = (dx2 + dx3)*cden1
      dcl2 = (dx3 + dx1)*cden2
      dcl3 = (dx1 + dx2)*cden3
      if (ldh) then
         electron_degeneracy_parameter = cl1*fermi_table_eta(id1) + &
              cl2*fermi_table_eta(id2) + cl3*fermi_table_eta(id3)
         if (electron_degeneracy_parameter .lt. etadh0) then
            cmfdh = cmfdh0*sqrt(cmfdh0)
         else if (electron_degeneracy_parameter .gt. etadh1) then
            cmfdh = cmfdh0*sqrt(cmfdh1)
         else
            t0 = cmfdh0*sqrt(cmfdh0)
            t1 = cmfdh0*sqrt(cmfdh1)
            slope = (t1-t0)/(etadh1-etadh0)
            cmfdh = t0+(electron_degeneracy_parameter-etadh0)*slope
         end if
      end if
!  NOW FIND PE AND ITS 1ST AND 2ND DERIVS.W/R/T DENSITY
      pel =  cl1*ff(1,1) +  cl2*ff(1,2) +  cl3*ff(1,3)
      dpel1 = dcl1*ff(1,1) + dcl2*ff(1,2) + dcl3*ff(1,3)
      dl8 = dml - emul8
      d8 = dexp(ln10*dl8)
      pa = ramut*d8
      pe = dexp(ln10*pel)
      if (ldh) then
         pdh = cdh*cmfdh*sqrt(d8/temperature)*d8
      end if
      pt = pa + pe + pr
      if (ldh) then
         pt = pt + pdh
         pdhp = pdh/pt
      end if
      ptl = dlog10(pt)
      pr7 = (pa + pe*dpel1)/pt
      if (ldh) then
         pr7 = pr7 + 1.5d0*pdhp
      end if
      pr8 = (ptl - pl8)/pr7
      pr9 = ln10*( (pa + pe*(dpel1**2+dpel2/ln10))/(pt*pr7) - pr7 )
      if (ldh) then
         pr9 = pr9 + 2.25d0*pdhp/pr7
      end if
! FACTOR IS THE DIFFERENCE BETWEEN THE P PREDICTED FROM THE GIVEN RHO
! AND THE TRUE P; IF TOO LARGE INCREMENT RHO AND TRY AGAIN.
      factor=dabs(pr8)
      if(factor .lt. 1d-15) factor = 1d-15
      pr8l=dlog10(factor)
      pr9l=dlog10(dabs(pr9))
      factl=dlog10(0.5d0)+pr9l+2.d0*pr8l
      fact=dexp(ln10*factl)
      fact=dsign(fact,pr9)
      corr=-pr8-fact
!  limit range of changes in density to a factor of 10**4   llp  1/31/07
      if(corr .gt.  4d0) corr =  4d0
      if(corr .lt. -4d0) corr = -4d0
      if(dabs(corr).ge.1.0d-08) then
       dml = dml + corr
       nden = nden + 1
       if(nden.le.20) go to 10
       write(short_file_unit,40) log10_temperature,log10_pressure,ptl,dml,corr
   40    format('EQRELV: Did not Converge: T,P,Pcalc,Dcalc,CORR', &
                4F10.6,F20.12)
!       PAUSE
       return
!         STOP 'ERRELV failed'
      end if
      log10_density = dl8
      density = d8
      do 60 kk = 2,5
       do 50 nn = 1,3
          id = id1 + nn - 1
          df1(kk,nn) = fermi_table_data(kk,id,jt2) - fermi_table_data(kk,id,jt1)
          df2(kk,nn)=(fermi_table_data(kk,id,jt3)-fermi_table_data(kk,id,jt2))-df1(kk,nn)
          ff(kk,nn)=fermi_table_data(kk,id,jt1)+ttud*df1(kk,nn)+ttcu*df2(kk,nn)
   50    continue
       ffe(kk) = cl1*ff(kk,1) + cl2*ff(kk,2) + cl3*ff(kk,3)
   60 continue
! DERIVATIVES OF P
      peq = dexp(ln10*ffe(3))
      qped = ffe(4)
      qqpedd =(dcl1*ff(4,1) + dcl2*ff(4,2) + dcl3*ff(4,3))*clni
      pr1 = df1(3,1)+ttdcu*df2(3,1)
      pr2 = df1(3,2)+ttdcu*df2(3,2)
      pr3 = df1(3,3)+ttdcu*df2(3,3)
      pr4 = dcl1*ff(3,1) + dcl2*ff(3,2) + dcl3*ff(3,3)
      pr5 =   5.0d0*(cl1*pr1 + cl2*pr2 + cl3*pr3)
      qpeqd = peq*pr4
      qpeqt = peq*(pr5 - 1.50d0*pr4)
      qpt = (pa + peq)*pinv + 4.0d0*beta1
      if (ldh) then
         qpt = qpt - 0.5d0*pdhp
      end if
      qpd = (pa + pe*qped)*pinv
      if (ldh) then
         qpd = qpd + 1.5d0*pdhp
      end if
      qqptt = (pa + qpeqt)*pinv + 16.0d0*beta1 - qpt**2
      if (ldh) then
         qqptt = qqptt+(0.5d0*qpt+0.25d0)*pdhp
      end if
      qqpdt = (pa + qpeqd)*pinv - qpt*qpd
      if (ldh) then
         qqpdt = qqpdt+(0.5d0*qpt-0.75d0)*pdhp
      end if
      qqpdd = (pa + pe*(qped**2 + qqpedd))*pinv - qpd**2
      if (ldh) then
        qqpdd = qqpdd + (-1.5d0*qpd+2.25d0)*pdhp
      end if
! DERIVATIVES OF DENSITY
      dlnrho_dlnp = 1.0d0/qpd
      dlnrho_dlnt = -qpt*dlnrho_dlnp
      pr1 = dlnrho_dlnp**2
      qdpp = -qqpdd*dlnrho_dlnp*pr1
      dlnrho_dlnp_dt = -qpt*qdpp - qqpdt*pr1
      dlnrho_dlnt_dt = -qpt*dlnrho_dlnp_dt - dlnrho_dlnp*(qqptt + dlnrho_dlnt*qqpdt)
! INTERNAL ENERGY TERM(U) AND ITS DERIVATIVES
      ue  = dexp(ln10*(ffe(2) + emul))
      ueq = dexp(ln10*(ffe(5) + emul))
      pr1 = df1(5,1)+ttdcu*df2(5,1)
      pr2 = df1(5,2)+ttdcu*df2(5,2)
      pr3 = df1(5,3)+ttdcu*df2(5,3)
      pr4 = dcl1*ff(5,1) + dcl2*ff(5,2) + dcl3*ff(5,3)
      pr5 =   5.0d0*(cl1*pr1 + cl2*pr2 + cl3*pr3)
      queqt = ueq*(pr5 - 1.50d0*pr4)
      pd = pressure/density
      pdt = pd/temperature
      ua = 1.50d0*ramut
      ur = 3.0d0*beta1*pd
      if (ldh) then
         udh = 3.0d0*pdh/d8
      end if
      u = ua + ue + ur
      if (ldh) then
         u = u + udh
         udhu = udh/u
      end if
      ui = 1.0d0/u
      pdu = pd*ui
      qutd = (ua + ueq + 4.0d0*ur)*ui
      if(ldh) then
         qutd = qutd -0.5d0*udhu
      end if
      qudt = pdu*(1.0d0 - qpt)
      if (ldh) then
         qudt = qudt + 0.5d0*udhu
      end if
      qutp = qutd + dlnrho_dlnt*qudt
      qqutt = (ua + queqt + 16.0d0*ur)*ui - qutd**2
      if (ldh) then
         qqutt = qqutt + (0.5d0*qutd+0.25d0)*udhu
      end if
      qqudt = qudt*(qpt - qutd) - pdu*qqptt
      if (ldh) then
         qqudt = qqudt + (0.5d0*qudt-0.25d0)*udhu
      end if
! SPECIFIC HEAT(QCP) AND ITS DERIVATIVES
      pr1 =-pdt*dlnrho_dlnt
      pr2 = u/temperature
      pr3 = -(1.0d0 + dlnrho_dlnt) +  dlnrho_dlnt_dt/dlnrho_dlnt
      pr4 = (1.0d0 - dlnrho_dlnp) +  dlnrho_dlnp_dt/dlnrho_dlnt
      acp = pr1*qpt
      bcp = pr2*qutd
      specific_heat_cp = acp + bcp
      qcpi = 1.0d0/specific_heat_cp
      specific_heat_cp_dt = acp*pr3 + pr1*(qqptt + dlnrho_dlnt*qqpdt) + &
              bcp*(qutp - 1.0d0) + pr2*(qqutt + dlnrho_dlnt*qqudt)
      specific_heat_cp_dp = acp*pr4 + dlnrho_dlnp*(pr1*qqpdt + bcp*qudt + pr2*qqudt)
      specific_heat_cp_dt = specific_heat_cp_dt*qcpi
      specific_heat_cp_dp = specific_heat_cp_dp*qcpi
! ADIABATIC GRADIENT(DELA) AND DERIVATIVES ( ALL DERIVS ARE LN(DELA)/LN(
      adiabatic_gradient = pr1*qcpi
      adiabatic_gradient_dt = pr3 - specific_heat_cp_dt
      adiabatic_gradient_dp = pr4 - specific_heat_cp_dp
      electron_degeneracy_parameter = cl1*fermi_table_eta(id1) + &
           cl2*fermi_table_eta(id2) + cl3*fermi_table_eta(id3)

      return
end subroutine eqrelv
