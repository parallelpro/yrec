!----------------------------------------------------------------------
! condopacp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original condopacp.f; only variable names, source form, and comment
! style were updated. This file holds three subroutines: condopacp
! (the Potekhin conductive-opacity table lookup, called from
! condopacpint.f90), and its two internal helpers cinterp3 (cubic/
! quadratic 1-D interpolation with derivatives) and hunt (Numerical
! Recipes' bisection/hunt table search).
!
!  This subroutine interpolates the electron thermal conductivity
!               from the data file "condall.d"
!  ---------------------------------------------------  Version 23.05.99
! Input: ion_charge - ion charge, log10_temperature - lg(T[K]),
!        log10_density - lg(rho[g/cc])
! Output: log10_conductivity - Log_{10} thermal conductivity (kappa) [CGS units]
!         dlnkappa_dlnrho - d log kappa / d log rho
!         dlnkappa_dlnt - d log kappa / d log T
!!!      (it is also possible to obtain all second derivatives)      !!!
!
!     This Subroutine interpolates thermal conductivity from
!     a precalculated table condall.d available at
!     http://www.ioffe.rssi.ru/astro/conduct
!     For theoretical base of this calculation, see
!      A.Y.Potekhin, D.A.Baiko, P.Haensel, D.G.Yakovlev, 1999
!                  Astron. Astrophys. 346, 345.
!     Extension from strongly- to weakly-degenerate regime
!     has been done using the thermal averaging - for example
!         A.Y.Potekhin, 1999, Astron. Astrophys. 351, 787.
!     Please quote these publication when using this program
!     SOURCE DATA FILE 'condall.d' MUST BE IN YOUR
!      CURRENT DIRECTORY
!     Please address any questions/comments to Alexander Potekhin:
!     e-mail: palex@astro.ioffe.rssi.ru
subroutine condopacp(ion_charge, log10_temperature, log10_density, &
     log10_conductivity, dlnkappa_dlnrho, dlnkappa_dlnt, ierr)
      use star_info_lib, only: star

      use conductive_table_lib, only: cond_table, cond_n_temp, cond_n_rho, cond_n_z, cond_table_loaded_marker
      use math_lib
      implicit none

      double precision, intent(in) :: ion_charge, log10_temperature, &
           log10_density
      double precision, intent(out) :: log10_conductivity, &
           dlnkappa_dlnrho, dlnkappa_dlnt
      integer, intent(out) :: ierr
! (the grid sizes cond_n_* and the table itself, cond_table, live in
! kap/conductive_table_lib.f90 and must be consistent with "condall.d")
      integer :: file_unit, z_index, t_index, r_index
      double precision :: z_grid_value, log10_ion_charge
      integer :: t_index_lo, t_index_hi, r_index_lo, r_index_hi
      double precision :: cktmz0, drktmz0, dr2ktmz0, xr
      double precision :: ckt0z0, drkt0z0, dr2kt0z0
      double precision :: ckt1z0, drkt1z0, dr2kt1z0
      double precision :: cktpz0, drktpz0, dr2ktpz0
      double precision :: cktmz1, drktmz1, dr2ktmz1
      double precision :: ckt0z1, drkt0z1, dr2kt0z1
      double precision :: ckt1z1, drkt1z1, dr2kt1z1
      double precision :: cktpz1, drktpz1, dr2ktpz1
      double precision :: xz1, xz0
      double precision :: cktm, drktm, dr2ktm
      double precision :: ckt0, drkt0, dr2kt0
      double precision :: ckt1, drkt1, dr2kt1
      double precision :: cktp, drktp, dr2ktp
      double precision :: dt2k, drtk, drt2k, xt

      if (cond_table%table_loaded_flag.ne.cond_table_loaded_marker) then   ! Reading
         file_unit = star%ctrl%icondopacp
! MHP 8/25 file opening moved to parmin
!         open(IP,file=FcondOpacP,status='OLD')
!         print*,'Reading thermal conductivity data...'
         read(file_unit,'(A)') ! skip the first line
        do z_index=1,cond_n_z
           read(file_unit,*) z_grid_value,(cond_table%temp_grid(t_index),t_index=1,cond_n_temp)
           cond_table%z_grid(z_index)=log10(z_grid_value)
          do r_index=1,cond_n_rho
             read(file_unit,*) cond_table%rho_grid(r_index), &
                  (cond_table%log10_kappa_table(t_index,r_index,z_index),t_index=1,cond_n_temp)
          enddo
      enddo
         close(file_unit)
         cond_table%table_loaded_flag=cond_table_loaded_marker
! KC 2025-05-30 fixed -Winteger-division
         z_index=int(cond_n_z/2.)+1
         t_index=int(cond_n_temp/2.)+1
         r_index=int(cond_n_rho/2.)+1
!         print*,'Potekhin Conductivity File read in.'
      endif
      ierr = 0
      log10_ion_charge=log10(ion_charge)
      call hunt(cond_table%z_grid,cond_n_z,log10_ion_charge,z_index)
      if (z_index.eq.0.or.z_index.eq.cond_n_z) then
! 2026 audit (section 7): stop -> ierr (runtime table-range guard;
! see kap_lib's kap_get).
         write(*,*) 'CONINTER: Z out of range'
         ierr = 1
         return
      end if
      call hunt(cond_table%temp_grid,cond_n_temp,log10_temperature,t_index)

!      if (IT.eq.0.or.IT.eq.MAXT) stop 'CONINTER: T out of range'
      if (t_index.eq.0.or.t_index.eq.cond_n_temp) then
                  print*, cond_table%z_grid
                  print*, cond_n_z, log10_temperature, z_index
! 2026 (ROADMAP.md stage 3): stop 'CONINTER: T out of range'
! converted to ierr (see kap_lib's kap_get); message preserved.
                  print*, 'CONINTER: T out of range'
                  ierr = 1
                  return
      endif

      call hunt(cond_table%rho_grid,cond_n_rho,log10_density,r_index)
      if (r_index.eq.0.or.r_index.eq.cond_n_rho) then
! 2026 audit (section 7): stop -> ierr (runtime table-range guard;
! see kap_lib's kap_get).
         write(*,*) 'CONINTER: rho out of range'
         ierr = 1
         return
      end if
      t_index_lo=max0(1,t_index-1)
      t_index_hi=min0(cond_n_temp,t_index+2)
      r_index_lo=max0(1,r_index-1)
      r_index_hi=min0(cond_n_rho,r_index+2)
! Cubic interpolation in RLG:
! Z0:
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_lo,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index_lo,r_index,z_index), &
           cond_table%log10_kappa_table(t_index_lo,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index_lo,r_index_hi,z_index), &
           cktmz0,drktmz0,dr2ktmz0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index,r_index,z_index), &
           cond_table%log10_kappa_table(t_index,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index,r_index_hi,z_index), &
           ckt0z0,drkt0z0,dr2kt0z0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index+1,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index+1,r_index,z_index), &
           cond_table%log10_kappa_table(t_index+1,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index+1,r_index_hi,z_index), &
           ckt1z0,drkt1z0,dr2kt1z0,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_hi,r_index_lo,z_index), &
           cond_table%log10_kappa_table(t_index_hi,r_index,z_index), &
           cond_table%log10_kappa_table(t_index_hi,r_index+1,z_index), &
           cond_table%log10_kappa_table(t_index_hi,r_index_hi,z_index), &
           cktpz0,drktpz0,dr2ktpz0,xr,ierr)
      if (ierr /= 0) return
! Z1:
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_lo,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index_lo,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index_lo,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index_lo,r_index_hi,z_index+1), &
           cktmz1,drktmz1,dr2ktmz1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index,r_index_hi,z_index+1), &
           ckt0z1,drkt0z1,dr2kt0z1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index+1,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index+1,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index+1,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index+1,r_index_hi,z_index+1), &
           ckt1z1,drkt1z1,dr2kt1z1,xr,ierr)
      if (ierr /= 0) return
      call cinterp3(cond_table%rho_grid(r_index_lo),cond_table%rho_grid(r_index),cond_table%rho_grid(r_index+1), &
           cond_table%rho_grid(r_index_hi),log10_density,r_index,cond_n_rho, &
           cond_table%log10_kappa_table(t_index_hi,r_index_lo,z_index+1), &
           cond_table%log10_kappa_table(t_index_hi,r_index,z_index+1), &
           cond_table%log10_kappa_table(t_index_hi,r_index+1,z_index+1), &
           cond_table%log10_kappa_table(t_index_hi,r_index_hi,z_index+1), &
           cktpz1,drktpz1,dr2ktpz1,xr,ierr)
      if (ierr /= 0) return
! Linear interpolation in ZLG:
      xz1=(log10_ion_charge-cond_table%z_grid(z_index))/(cond_table%z_grid(z_index+1)-cond_table%z_grid(z_index))
      xz0=1d0-xz1
      cktm=xz0*cktmz0+xz1*cktmz1
      drktm=xz0*drktmz0+xz1*drktmz1
      dr2ktm=xz0*dr2ktmz0+xz1*dr2ktmz1
      ckt0=xz0*ckt0z0+xz1*ckt0z1
      drkt0=xz0*drkt0z0+xz1*drkt0z1
      dr2kt0=xz0*dr2kt0z0+xz1*dr2kt0z1
      ckt1=xz0*ckt1z0+xz1*ckt1z1
      drkt1=xz0*drkt1z0+xz1*drkt1z1
      dr2kt1=xz0*dr2kt1z0+xz1*dr2kt1z1
      cktp=xz0*cktpz0+xz1*cktpz1
      drktp=xz0*drktpz0+xz1*drktpz1
      dr2ktp=xz0*dr2ktpz0+xz1*dr2ktpz1
! Cubic interpolation in TLG:
      call cinterp3(cond_table%temp_grid(t_index_lo),cond_table%temp_grid(t_index),cond_table%temp_grid(t_index+1), &
           cond_table%temp_grid(t_index_hi),log10_temperature,t_index,cond_n_temp, &
           cktm,ckt0,ckt1,cktp, & ! input: values of lg kappa
           log10_conductivity,dlnkappa_dlnt,dt2k,xt,ierr) ! lg kappa, d lg k / d lg T, d2 lg k / d2 lg T
      if (ierr /= 0) return
      call cinterp3(cond_table%temp_grid(t_index_lo),cond_table%temp_grid(t_index),cond_table%temp_grid(t_index+1), &
           cond_table%temp_grid(t_index_hi),log10_temperature,t_index,cond_n_temp, &
           drktm,drkt0,drkt1,drktp, & ! input: values of d lg k / d lg rho
           dlnkappa_dlnrho,drtk,drt2k,xt,ierr) ! d lg k / d lg rho, d2 lgk/(d lgT d lg rho)
      if (ierr /= 0) return
      return
end subroutine condopacp

! Given 4 values of Z and 4 values of V, find VF corresponding to 5th Z
!                                                       Version 23.05.99
!   Output: interp_value - interpolated value of function
!           interp_deriv - interpolated derivative
!           interp_deriv2 - interpolated second derivative
!           interp_frac - fraction of the path from grid_index to
!                         grid_index+1
subroutine cinterp3(grid_lo2,grid_lo1,grid_hi1,grid_hi2,grid_target, &
     grid_index,grid_size,val_lo2,val_lo1,val_hi1,val_hi2, &
     interp_value,interp_deriv,interp_deriv2,interp_frac, ierr)

      implicit none

      double precision, intent(in) :: grid_lo2, grid_lo1, grid_hi1, &
           grid_hi2, grid_target
      integer, intent(in) :: grid_index, grid_size
      double precision, intent(in) :: val_lo2, val_lo1, val_hi1, val_hi2
      double precision, intent(out) :: interp_value, interp_deriv, &
           interp_deriv2, interp_frac
      integer, intent(out) :: ierr

      double precision :: x, h, hm, v01, hp, v11, c2, c3

      ierr = 0
      if (grid_index.le.0.or.grid_index.ge.grid_size) then
! 2026 audit (section 7): stop -> ierr (runtime table-range guard;
! see kap_lib's kap_get).
         write(*,*) 'CINTERP: N0 out of range'
         ierr = 1
         return
      end if
      x=grid_target-grid_lo1
      h=grid_hi1-grid_lo1   ! basic interval
      interp_frac=x/h
      if (grid_index.gt.1) then
         hm=grid_lo1-grid_lo2  ! left adjoint interval
         v01=((val_hi1-val_lo1)/h**2+(val_lo1-val_lo2)/hm**2)/ &
         (1d0/h+1d0/hm) ! left derivative
      endif
      if (grid_index.lt.grid_size-1) then
         hp=grid_hi2-grid_hi1 ! right adjoint interval
         v11=((val_hi1-val_lo1)/h**2+(val_hi2-val_hi1)/hp**2)/ &
         (1d0/h+1d0/hp) ! right derivative
      endif
      if (grid_index.gt.1.and.grid_index.lt.grid_size-1) then   ! Cubic interpolation
         c2=3d0*(val_hi1-val_lo1)-h*(v11+2d0*v01)
         c3=h*(v01+v11)-2d0*(val_hi1-val_lo1)
         interp_value=val_lo1+v01*x+c2*interp_frac**2+c3*interp_frac**3
         interp_deriv=v01+(2d0*c2*interp_frac+3d0*c3*interp_frac**2)/h
         interp_deriv2=(2d0*c2+6d0*c3*interp_frac)/h**2
         return
      endif
      if (grid_index.eq.1) then   ! Quadratic interpolation
         c2=val_lo1-val_hi1+v11*h
         interp_value=val_hi1-v11*(h-x)+c2*(1d0-interp_frac)**2
         interp_deriv=v11-2d0*c2*(1d0-interp_frac)/h
         interp_deriv2=2d0*c2/h**2
      else  ! N0=MXNV-1
         c2=val_hi1-val_lo1-v01*h
         interp_value=val_lo1+v01*x+c2*interp_frac**2
         interp_deriv=v01+2d0*c2*interp_frac/h
         interp_deriv2=2d0*c2/h**2
      endif
      return
end subroutine cinterp3

!   W.H.Press, B.P.Flannery, S.A.Teukolsky, W.T.Vetterling
!   Numerical Receipes(Cambridge Univ., 1986)
!     Given an array XX of length N, and given a value X,
!     returns a value JLO such that X is between XX(JLO) and XX(JLO+1).
!     XX must be monotonic, either increasing or decreasing.
!     JLO=0 or JLO=N is returned to indicate that X is out of range.
!     JLO on input is taken as the initial guess for JLO on output.
subroutine hunt(table,n,target_value,index)

      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: table(*)
      double precision, intent(in) :: target_value
      integer, intent(inout) :: index

      logical :: ascending
      integer :: index_hi, increment, index_mid

      ascending=table(n).gt.table(1) ! true if ascending order, false otherwise
      if (index.le.0.or.index.gt.n) then ! Input guess not useful.
         index=0
         index_hi=n+1  ! go immediately to bisection
      else
      increment=1 ! set the hunting increment
      if (target_value.ge.table(index).eqv.ascending) then ! Hunt up:
        do
          index_hi=index+increment
        if (index_hi.gt.n) then ! Done hunting, since off end of table
           index_hi=n+1
           exit
        elseif (target_value.ge.table(index_hi).eqv.ascending) then ! Not done hunting
           index=index_hi
           increment=increment+increment
        else
           exit
        endif
        end do
      else ! Hunt down:
         index_hi=index
         do
          index=index_hi-increment
        if (index.lt.1) then ! Done hunting, since off end of table
           index=0
           exit
        elseif (target_value.lt.table(index).eqv.ascending) then ! Not done hunting
           index_hi=index
           increment=increment+increment ! so double the increment
        else
           exit
        endif ! Done hunting, value bracketed
         end do
      endif
      endif
!   Hunt is done, so begin the final bisection phase:
      do
      if (index_hi-index.eq.1) return
      index_mid=(index_hi+index)/2
      if (target_value.ge.table(index_mid).eqv.ascending) then
         index=index_mid
      else
         index_hi=index_mid
      endif
      end do
end subroutine hunt
