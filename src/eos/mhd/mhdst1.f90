!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MHDST1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------------------------------------------------------------------
! mhdst1
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdst1.f; only variable names, source form, and comment
! style were updated.
!
! Reads one MHD equation-of-state table unit. table_kind=0 reads a
! single ZAMS-type table pair (lower/upper T region). table_kind=1
! reads a centre-type table plus the down/up-composition tables
! needed to build numerical X-derivatives, and stores the extended
! variable set (indices ivar1+1..ivar2) used in T-rho regions with
! inhomogeneous composition.
subroutine mhdst1(table_unit,table_kind,nt1m,nr1m,ivar1,nt2m,nr2m,ivar2,nchem0, &
                  num_t1,num_r1,num_t2,num_r2,log10t1,log10t2,table_vars1,table_vars2, &
                  drho1,drho2,num_chem,atomic_weight,number_abundance,mass_fraction,mean_molecular_weight, &
                  log10t_down,log10t_up,table_vars_centroid,table_vars_down,table_vars_up, &
                  atomic_weight_down,atomic_weight_up, &
                  number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
      use const_lib
      implicit none
      integer, intent(in) :: table_unit, table_kind, nt1m, nr1m, ivar1, &
           nt2m, nr2m, ivar2, nchem0
      integer, intent(inout) :: num_t1, num_r1, num_t2, num_r2, num_chem
      double precision, intent(inout) :: log10t1(nt1m), &
           table_vars1(nt1m,nr1m,ivar1), &
           log10t2(nt2m), table_vars2(nt2m,nr2m,ivar2)
      double precision, intent(inout) :: drho1, drho2
      double precision, intent(inout) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0), mean_molecular_weight

      double precision :: log10t_down(nt2m), log10t_up(nt2m)
      double precision :: table_vars_centroid(nt2m,nr2m,ivar1)
      double precision :: table_vars_down(nt2m,nr2m,ivar1)
      double precision :: table_vars_up(nt2m,nr2m,ivar1)
      double precision :: atomic_weight_down(nchem0), number_abundance_down(nchem0), &
           mass_fraction_down(nchem0)
      double precision :: atomic_weight_up(nchem0), number_abundance_up(nchem0), &
           mass_fraction_up(nchem0)
      integer :: num_composition_reads, composition_pass, num_vars_read, &
           table_kind_read, composition_flag_read, composition_flag_expected
      double precision :: delta_x, composition_tolerance
      double precision :: unused_mean_molecular_weight_down, &
           unused_mean_molecular_weight_up
      integer :: species_index, temp_check_index, temp_deriv_index, &
           density_deriv_index, var_index

      integer, intent(out) :: ierr

      ierr = 0

      if (table_kind.eq.0) then
          num_composition_reads = 1
      else
          num_composition_reads = 3
      end if
      do composition_pass=1,num_composition_reads
!     READ(IR,98,END=1000) IVARR,IDXR,IRESCR,DDX
      read(table_unit,   end=1000) num_vars_read,table_kind_read,composition_flag_read,delta_x
      if (ivar1.lt.num_vars_read) then
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (table_kind.ne.table_kind_read) then
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (composition_pass.eq.1) then
          composition_flag_expected= 0
      else if (composition_pass.eq.2) then
          composition_flag_expected=-1
      else if (composition_pass.eq.3) then
          composition_flag_expected= 1
      end if
      if (composition_flag_expected.ne.composition_flag_read) then
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      if (composition_pass.eq.1) call rabu(table_unit,nchem0,num_chem,atomic_weight,number_abundance,mass_fraction,mean_molecular_weight,ierr)
      if (composition_pass.eq.2) call rabu(table_unit,nchem0,num_chem,atomic_weight_down,number_abundance_down,mass_fraction_down,unused_mean_molecular_weight_down,ierr)
      if (composition_pass.eq.3) call rabu(table_unit,nchem0,num_chem,atomic_weight_up,number_abundance_up,mass_fraction_up,unused_mean_molecular_weight_up,ierr)
      if (ierr /= 0) return
!     READ(IR,1001) NT1,NT2,DRH1,DRH2
      read(table_unit     ) num_t1,num_t2,drho1,drho2
      if (table_kind.eq.1 .and. num_t1.ne.0) then
          ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
          ! facades stop when their caller passes no ierr.
          ierr = 1
          return
      end if
      if (num_t1.gt.0) then
         call rtab(table_unit,nt1m,nr1m,ivar1,num_t1,num_r1,log10t1,table_vars1, ierr)
         if (ierr /= 0) return
      end if
      if (table_kind.eq.1) then
       if (composition_pass.eq.1) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t2,table_vars_centroid, ierr)
       if (ierr /= 0) return
       if (composition_pass.eq.2) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t_down,table_vars_down, ierr)
       if (ierr /= 0) return
       if (composition_pass.eq.3) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t_up,table_vars_up, ierr)
       if (ierr /= 0) return
      else if (table_kind.eq.0) then
       call rtab(table_unit,nt2m,nr2m,ivar2,num_t2,num_r2,log10t2,table_vars2, ierr)
       if (ierr /= 0) return
      end if
 400  continue
      end do
      if (.not. (table_kind.eq.0)) then
!     IF IDX=1: CHECK TABLES FOR CORRECT COMPOSITION CONSTRUCTION
!     AND PERFORM NUMERICAL DERIVATIVES W.R.T. X
      composition_tolerance = 0.05d0*abs(delta_x)
      if ( abs(mass_fraction(1)-mass_fraction_down(1)-delta_x).gt.composition_tolerance  .or. abs(mass_fraction(1)-mass_fraction_up(1)+delta_x).gt.composition_tolerance  .or. abs(mass_fraction(2)-mass_fraction_down(2)+delta_x).gt.composition_tolerance  .or. abs(mass_fraction(2)-mass_fraction_up(2)-delta_x).gt.composition_tolerance ) then
         continue
         
         
         ierr = 1
         return
      end if
      do species_index=3,num_chem
      if ( abs(mass_fraction(species_index)-mass_fraction_up(species_index)).gt.composition_tolerance ) then
         continue
         
         
         ierr = 1
         return
      end if
      if ( abs(mass_fraction(species_index)-mass_fraction_down(species_index)).gt.composition_tolerance ) then
         continue
         
         
         ierr = 1
         return
      end if
 420  continue
      end do
      do temp_check_index=1,num_t2
      if (log10t2(temp_check_index).ne.log10t_down(temp_check_index)) then
         continue
         
         
         ierr = 1
         return
      end if
      if (log10t2(temp_check_index).ne.log10t_up(temp_check_index)) then
         continue
         
         
         ierr = 1
         return
      end if
 430  continue
      end do
!     NUMERICAL DERIVATIVES W.R.T. X
      do temp_deriv_index =1,num_t2
! KC 2025-05-30 fixed "Shared DO termination label"
!       DO 440 M =1,NR2
      do density_deriv_index =1,num_r2
      do var_index=1,ivar1
! KC 2025-05-30 fixed "DO termination statement which is not END DO or CONTINUE"
! 435   TDVAR2(N,M,IV)=TDDIF0(N,M,IV)
        table_vars2(temp_deriv_index,density_deriv_index,var_index)=table_vars_centroid(temp_deriv_index,density_deriv_index,var_index)
435   continue
      end do
!
!     EXTENDED SET OF VARIABLES (TDVAR2(N,M,IVAR1+1...IVAR2))
!     FOR T-RHO REGIONS WITH INHOMOGENEOUS COMPOSITION.
!     IN THE COMMENTS,R AND T DENOTE LOG10(RHO) AND LOG10(T).
!     DLOG10(P)/DX,DLOG10(U)/DX,DDELAD/DX,DLOG10(CP)/DX
      table_vars2(temp_deriv_index,density_deriv_index,21)=(table_vars_up(temp_deriv_index,density_deriv_index, 2)-table_vars_down(temp_deriv_index,density_deriv_index, 2))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,22)=(table_vars_up(temp_deriv_index,density_deriv_index, 3)-table_vars_down(temp_deriv_index,density_deriv_index, 3))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,23)=(table_vars_up(temp_deriv_index,density_deriv_index, 8)-table_vars_down(temp_deriv_index,density_deriv_index, 8))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,24)=(table_vars_up(temp_deriv_index,density_deriv_index, 9)-table_vars_down(temp_deriv_index,density_deriv_index, 9))/(2.d0*delta_x)
!     SPACE-HOLDER VARIABLE (LIKE VAR(20))
      table_vars2(temp_deriv_index,density_deriv_index,25)=8888844444.d0
  441 continue
      end do
  440 continue
      end do
!     NORMAL EXIT
      end if
 450  continue
      return
!     ERROR EXIT AND ERROR MESSAGES
  500 continue
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return
 600  continue
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return
 1000 continue
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return
!  98   FORMAT(1X,3I5,F13.5)
!  99   FORMAT(1X,I5,(/1X,3E15.7))
! 1001  FORMAT(2I5,2F10.6)
! 8001  FORMAT(' CORRECT TABLE CONSTRUCTION FOR X-DERIVATIVES.',
!      1       ' CENTROID COMPOSITION IS:'//)
! 8002  FORMAT('      AT. WEIGHT     NUMBER ',
!      1 'ABUNDANCE  MASS FRACTION',(/1X,1P3G16.7))
! 8003  FORMAT(/' MEAN MOLECULAR WEIGHT = ',F12.7//)
! 9006  FORMAT(' ERROR IN MHDST1. IVARR READ FROM TABLE IS',
!      1 ' BIGGER THAN THE VALUE USED IN THE COMMONS.',
!      2 ' IVAR,IVARR= ',/1X,2I8)
! 9007  FORMAT(' ERROR IN MHDST1. IDXR READ FROM TABLE IS INCORRECT',
!      1 ' IDX,IDXR= ',/1X,2I8)
! 9008  FORMAT(' ERROR IN MHDST1. IRESCR READ FROM TABLE IS INCORRECT',
!      1 ' IRESCO,IRESCR= ',/1X,2I8)
! 9010  FORMAT(' ERROR IN MHDST1. NT1 AND IDX ARE INCONSISTENT',
!      1 ' IDX,NT1,NT2 ',/1X,3I8)
! 9800  FORMAT(' ERROR IN TABLE CONSTRUCTION FOR X-DERIVATIVES',
!      1       ' CENTRAL, LOWER, UPPER TABLE: N(ELEMENT),ABFRCS(N)'//)
! 9810  FORMAT(1X,I5,F15.9)
! 9820  FORMAT(/)
! 9850  FORMAT(' ERROR IN TABLE CONSTRUCTION FOR X-DERIVATIVES:',
!      1       ' TEMPERATURES WRONG: J,TLOW(J),TCENT(J),TUPP(J)'//)
! 9860  FORMAT(1X,I5,3F15.9)
! 9900  FORMAT(' EOF REACHED IN INPUT FILE. ERROR STOP. IR,IDX = ',2I5)
end subroutine mhdst1
