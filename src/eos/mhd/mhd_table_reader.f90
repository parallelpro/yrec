!----------------------------------------------------------------------
! mhd_table_reader
!----------------------------------------------------------------------
! The MHD tables' reader, read_mhd_tables (eos_init calls it, through
! MHD's table-source entry), and Daeppen's reading routines under their
! original names (mhdst, mhdst1, rabu, rtab): the files' I/O, kept apart
! from the MHD EOS that interpolates them (mhd_eos prints nothing). The
! tables live in mhd_eos_tables (mhd).
module mhd_table_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use mhd_eos_tables, only: i_mhd_dgrad_ad_dx, i_mhd_dlog10cp_dx, i_mhd_dlog10p_dx, &
           i_mhd_dlog10u_dx, i_mhd_grad_ad, i_mhd_log10_cp, i_mhd_log10_p, i_mhd_log10_u, &
           i_mhd_placeholder, mhd, mhd_ivarc, mhd_ivarx, mhd_n_centre, mhd_n_zams, &
           mhd_nchem0, mhd_nr1m, mhd_nr2m, mhd_nrxm, mhd_nt1m, mhd_nt2m, mhd_ntxm
      use eos_messages, only: eos_report, message_length, unwritten_line
      implicit none
      private
      public :: read_mhd_tables

contains

!----------------------------------------------------------------------
! read_mhd_tables (formerly mhdtbl)
!----------------------------------------------------------------------
! Opens the 8 MHD equation-of-state table files (the 3 ZAMS-type tables
! for the compositions A, B and C, and the 5 centre-type tables), each
! on a new unit, reads them with mhdst and closes them. A file that
! cannot be opened, or a table that fails mhdst's checks, returns
! ierr = 1 with a message to the terminal and the run log.
subroutine read_mhd_tables(zams_paths, centre_paths, ierr)
      character(len=*), intent(in) :: zams_paths(mhd_n_zams), centre_paths(mhd_n_centre)
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      integer :: unit_zams(mhd_n_zams), unit_centre(mhd_n_centre), ic

      ierr = 0
      unit_zams = 0
      unit_centre = 0
      do ic = 1, mhd_n_zams
         if (ierr == 0) call open_mhd_table(trim(zams_paths(ic)), unit_zams(ic), ierr)
      end do
      do ic = 1, mhd_n_centre
         if (ierr == 0) call open_mhd_table(trim(centre_paths(ic)), unit_centre(ic), ierr)
      end do
      if (ierr == 0) then
         call mhdst(unit_zams, unit_centre, ierr)
         if (ierr /= 0) then
            eos_lines = unwritten_line
            write(eos_lines,*) 'eos_init: an MHD table is corrupt or does not match its layout'
            call eos_report(eos_lines)
         end if
      end if
      do ic = 1, mhd_n_zams
         if (unit_zams(ic) /= 0) close(unit_zams(ic))
      end do
      do ic = 1, mhd_n_centre
         if (unit_centre(ic) /= 0) close(unit_centre(ic))
      end do
end subroutine read_mhd_tables

! open one MHD table (unformatted) on a new unit; ierr = 1, with a
! message to the terminal and the run log, when it cannot be opened
subroutine open_mhd_table(path, table_unit, ierr)
      character(len=*), intent(in) :: path
      integer, intent(out) :: table_unit, ierr
      character(len=message_length) :: eos_lines(4)
      open(newunit=table_unit, file=path, status='OLD', form='UNFORMATTED', &
           iostat=ierr)
      if (ierr /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: cannot open the MHD table ', trim(path)
         call eos_report(eos_lines)
         table_unit = 0
         ierr = 1
      end if
end subroutine open_mhd_table

!----------------------------------------------------------------------
! mhdst
!----------------------------------------------------------------------
! Reads the 3 ZAMS-type tables (compositions A, B, C; unit_zams) and the
! 5 centre-type tables (unit_centre) into mhd_eos_tables's mhd and sets
! the temperature limits mhdpx1 uses: the lower ZAMS tables up to their
! last row, the upper ZAMS tables up to the first centre row, the
! centre tables above.
subroutine mhdst(unit_zams, unit_centre, ierr)

      integer, intent(in) :: unit_zams(mhd_n_zams), unit_centre(mhd_n_centre)

      integer, intent(out) :: ierr

      integer, parameter :: ivarc = mhd_ivarc
      integer, parameter :: ivarx = mhd_ivarx
      integer, parameter :: nchem0 = mhd_nchem0
      integer, parameter :: nt1m = mhd_nt1m
      integer, parameter :: nt2m = mhd_nt2m
      integer, parameter :: ntxm = mhd_ntxm
      integer, parameter :: nr1m = mhd_nr1m
      integer, parameter :: nr2m = mhd_nr2m
      integer, parameter :: nrxm = mhd_nrxm

      real(dp) :: unused_zams_table(nt1m,nr1m,ivarc), unused_zams_log10t(nt1m)
!     Working storage for numerical X-derivatives. Common to all
!     X-TABLES
      real(dp) :: log10t_down(nt2m), log10t_up(nt2m)
      real(dp) :: table_centre_vars(nt2m,nr2m,ivarc)
      real(dp) :: table_down_vars(nt2m,nr2m,ivarc)
      real(dp) :: table_up_vars(nt2m,nr2m,ivarc)
      real(dp) :: atomic_weight_down(nchem0), number_abundance_down(nchem0), &
           mass_fraction_down(nchem0)
      real(dp) :: atomic_weight_up(nchem0), number_abundance_up(nchem0), &
           mass_fraction_up(nchem0)
!     Placeholder actuals for the mhdst1 arguments that the centre-table
!     calls do not need (the ZAMS lower-table slots).
      integer :: num_chem_species, unused_num_t, unused_num_r, table_index
      real(dp) :: unused_drho
      integer :: ic

      ierr = 0

      num_chem_species = 0
      unused_drho = 0.d0
      unused_num_t = 0
      unused_num_r = 0
!     Read ZAMS tables
      do ic = 1, mhd_n_zams
         table_index = 0
         call mhdst1(unit_zams(ic),table_index,nt1m,nr1m,ivarc,nt2m,nr2m,ivarc,nchem0, &
                     mhd%zams_lower_num_t,mhd%zams_lower_num_r,mhd%zams_upper_num_t,mhd%zams_upper_num_r, &
                     mhd%zams_lower_log10t,mhd%zams_upper_log10t, &
                     mhd%zams_lower_table(:,:,:,ic),mhd%zams_upper_table(:,:,:,ic), &
                     mhd%zams_lower_drho,mhd%zams_upper_drho,num_chem_species, &
                     mhd%zams_atomic_weight(:,ic),mhd%zams_number_abundance(:,ic), &
                     mhd%zams_mass_fraction(:,ic),mhd%zams_mean_molecular_weight(ic), &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
         if (ierr /= 0) return
      end do
!     Read centre tables
      do ic = 1, mhd_n_centre
         table_index = 1
         call mhdst1(unit_centre(ic),table_index,nt1m,nr1m,ivarc,ntxm,nrxm,ivarx,nchem0, &
                     unused_num_t,unused_num_r,mhd%centre_num_t,mhd%centre_num_r,unused_zams_log10t, &
                     mhd%centre_log10t,unused_zams_table,mhd%centre_table(:,:,:,ic), &
                     unused_drho,mhd%centre_drho,num_chem_species, &
                     mhd%centre_atomic_weight(:,ic),mhd%centre_number_abundance(:,ic), &
                     mhd%centre_mass_fraction(:,ic),mhd%centre_mean_molecular_weight(ic), &
                     log10t_down,log10t_up,table_centre_vars,table_down_vars,table_up_vars, &
                     atomic_weight_down,atomic_weight_up, &
                     number_abundance_down,number_abundance_up,mass_fraction_down,mass_fraction_up, ierr)
         if (ierr /= 0) return
      end do
! every family needs 4 rows and 4 density points for its 4-point
! interpolation
      if (min(mhd%zams_lower_num_t, mhd%zams_upper_num_t, mhd%centre_num_t, &
           mhd%zams_lower_num_r, mhd%zams_upper_num_r, mhd%centre_num_r) < 4) then
         ierr = 1
         return
      end if
!     Temperature limits
      mhd%table_log10t_min = mhd%zams_lower_log10t(1)
      mhd%zams_lower_upper_boundary_log10t = mhd%zams_lower_log10t(mhd%zams_lower_num_t)
      mhd%zams_centre_boundary_log10t = mhd%centre_log10t(1)
      mhd%table_log10t_max = mhd%centre_log10t(mhd%centre_num_t)
end subroutine mhdst

!----------------------------------------------------------------------
! mhdst1
!----------------------------------------------------------------------
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
      integer, intent(in) :: table_unit, table_kind, nt1m, nr1m, ivar1, &
           nt2m, nr2m, ivar2, nchem0
      integer, intent(inout) :: num_t1, num_r1, num_t2, num_r2, num_chem
      real(dp), intent(inout) :: log10t1(nt1m), &
           table_vars1(nt1m,nr1m,ivar1), &
           log10t2(nt2m), table_vars2(nt2m,nr2m,ivar2)
      real(dp), intent(inout) :: drho1, drho2
      real(dp), intent(inout) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0), mean_molecular_weight

      integer, intent(out) :: ierr

      real(dp) :: log10t_down(nt2m), log10t_up(nt2m)
      real(dp) :: table_vars_centroid(nt2m,nr2m,ivar1)
      real(dp) :: table_vars_down(nt2m,nr2m,ivar1)
      real(dp) :: table_vars_up(nt2m,nr2m,ivar1)
      real(dp) :: atomic_weight_down(nchem0), number_abundance_down(nchem0), &
           mass_fraction_down(nchem0)
      real(dp) :: atomic_weight_up(nchem0), number_abundance_up(nchem0), &
           mass_fraction_up(nchem0)
      integer :: num_composition_reads, composition_pass, num_vars_read, &
           table_kind_read, composition_flag_read, composition_flag_expected
      real(dp) :: delta_x, composition_tolerance
      real(dp) :: unused_mean_molecular_weight_down, &
           unused_mean_molecular_weight_up
      integer :: species_index, temp_check_index, temp_deriv_index, &
           density_deriv_index, var_index

      ierr = 0

      if (table_kind == 0) then
          num_composition_reads = 1
      else
          num_composition_reads = 3
      end if
      do composition_pass=1,num_composition_reads
      read(table_unit, iostat=ierr) num_vars_read,table_kind_read,composition_flag_read,delta_x
      if (ierr /= 0) then
         ierr = 1
         return
      end if
      if (ivar1 < num_vars_read) then
         ierr = 1
         return
      end if
      if (table_kind /= table_kind_read) then
         ierr = 1
         return
      end if
      if (composition_pass == 1) then
          composition_flag_expected= 0
      else if (composition_pass == 2) then
          composition_flag_expected=-1
      else if (composition_pass == 3) then
          composition_flag_expected= 1
      end if
      if (composition_flag_expected /= composition_flag_read) then
         ierr = 1
         return
      end if
      if (composition_pass == 1) call rabu(table_unit,nchem0,num_chem,atomic_weight,number_abundance,mass_fraction,mean_molecular_weight,ierr)
      if (composition_pass == 2) call rabu(table_unit,nchem0,num_chem,atomic_weight_down,number_abundance_down,mass_fraction_down,unused_mean_molecular_weight_down,ierr)
      if (composition_pass == 3) call rabu(table_unit,nchem0,num_chem,atomic_weight_up,number_abundance_up,mass_fraction_up,unused_mean_molecular_weight_up,ierr)
      if (ierr /= 0) return
      read(table_unit, iostat=ierr) num_t1,num_t2,drho1,drho2
      if (ierr /= 0) then
         ierr = 1
         return
      end if
      if (table_kind == 1 .and. num_t1 /= 0) then
          ierr = 1
          return
      end if
      if (num_t1 > 0) then
         call rtab(table_unit,nt1m,nr1m,ivar1,num_t1,num_r1,log10t1,table_vars1, ierr)
         if (ierr /= 0) return
      end if
      if (table_kind == 1) then
       if (composition_pass == 1) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t2,table_vars_centroid, ierr)
       if (ierr /= 0) return
       if (composition_pass == 2) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t_down,table_vars_down, ierr)
       if (ierr /= 0) return
       if (composition_pass == 3) call rtab(table_unit,nt2m,nr2m,ivar1,num_t2,num_r2,log10t_up,table_vars_up, ierr)
       if (ierr /= 0) return
      else if (table_kind == 0) then
       call rtab(table_unit,nt2m,nr2m,ivar2,num_t2,num_r2,log10t2,table_vars2, ierr)
       if (ierr /= 0) return
      end if
      end do
      if (table_kind /= 0) then
!     If idx=1: Check tables for correct composition construction
!     and perform numerical derivatives w.r.t. X
      composition_tolerance = 0.05d0*abs(delta_x)
      if ( abs(mass_fraction(1)-mass_fraction_down(1)-delta_x) > composition_tolerance  .or. abs(mass_fraction(1)-mass_fraction_up(1)+delta_x) > composition_tolerance  .or. abs(mass_fraction(2)-mass_fraction_down(2)+delta_x) > composition_tolerance  .or. abs(mass_fraction(2)-mass_fraction_up(2)-delta_x) > composition_tolerance ) then
         ierr = 1
         return
      end if
      do species_index=3,num_chem
      if ( abs(mass_fraction(species_index)-mass_fraction_up(species_index)) > composition_tolerance ) then
         ierr = 1
         return
      end if
      if ( abs(mass_fraction(species_index)-mass_fraction_down(species_index)) > composition_tolerance ) then
         ierr = 1
         return
      end if
      end do
      do temp_check_index=1,num_t2
      if (log10t2(temp_check_index) /= log10t_down(temp_check_index)) then
         ierr = 1
         return
      end if
      if (log10t2(temp_check_index) /= log10t_up(temp_check_index)) then
         ierr = 1
         return
      end if
      end do
!     Numerical derivatives w.r.t. X
      do temp_deriv_index =1,num_t2
      do density_deriv_index =1,num_r2
      do var_index=1,ivar1
        table_vars2(temp_deriv_index,density_deriv_index,var_index)=table_vars_centroid(temp_deriv_index,density_deriv_index,var_index)
      end do
!
!     Extended set of variables (TDVAR2(n,m,IVAR1+1...IVAR2))
!     for T-rho regions with inhomogeneous composition.
!     In the comments,R and T denote log10(rho) and log10(T).
!     log10(P)/DX,log10(U)/DX,DDELAD/DX,log10(CP)/DX
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dlog10p_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_log10_p)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_log10_p))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dlog10u_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_log10_u)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_log10_u))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dgrad_ad_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_grad_ad)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_grad_ad))/(2.d0*delta_x)
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_dlog10cp_dx)=(table_vars_up(temp_deriv_index,density_deriv_index,i_mhd_log10_cp)-table_vars_down(temp_deriv_index,density_deriv_index,i_mhd_log10_cp))/(2.d0*delta_x)
!     Space-holder variable (like var(20))
      table_vars2(temp_deriv_index,density_deriv_index,i_mhd_placeholder)=8888844444.d0
      end do
      end do
      end if
end subroutine mhdst1

!----------------------------------------------------------------------
! rabu
!----------------------------------------------------------------------
! Reads one composition-table record (element count, atomic weights,
! number abundances, mass fractions, and mean molecular weight) from
! an unformatted file unit; a short record or more than nchem0
! elements is ierr = 1.
subroutine rabu(table_unit, nchem0, num_chem, atomic_weight, &
     number_abundance, mass_fraction, mean_molecular_weight, ierr)


      integer, intent(in) :: table_unit, nchem0
      integer, intent(out) :: num_chem
      real(dp), intent(out) :: atomic_weight(nchem0), &
           number_abundance(nchem0), mass_fraction(nchem0)
      real(dp), intent(out) :: mean_molecular_weight
      integer, intent(out) :: ierr

      integer :: ic

      ierr = 0

      read(table_unit, iostat=ierr) num_chem,(atomic_weight(ic),number_abundance(ic), &
           mass_fraction(ic), ic=1,num_chem),mean_molecular_weight
      if (ierr /= 0 .or. nchem0 < num_chem) then
         ierr = 1
         return
      end if
end subroutine rabu

!----------------------------------------------------------------------
! rtab
!----------------------------------------------------------------------
!
! Reads an unformatted opacity/EOS-style table of shape
! (num_t_points, num_rho_points, num_vars) from file_unit: for each of
! num_t_points temperature records, reads the record's density-point
! count (asserted equal for every record) and log10(T), then reads
! num_vars variables at each of its density points. A short record, a
! table larger than its arrays or a ragged one is ierr = 1.
subroutine rtab(file_unit,max_t_points,max_rho_points,num_vars, &
     num_t_points,num_rho_points,log_t,table_data, ierr)

      integer, intent(in) :: file_unit, max_t_points, max_rho_points, &
           num_vars, num_t_points
      integer, intent(out) :: num_rho_points

!     NT is input; NR,TL,tdvar are output
      real(dp), intent(out) :: log_t(max_t_points), &
           table_data(max_t_points,max_rho_points,num_vars)

      integer, intent(out) :: ierr
! --- locals ---
      integer :: t_idx, rho_idx, var_idx, rho_count_read

      ierr = 0

      do t_idx = 1, num_t_points
      read(file_unit, iostat=ierr) rho_count_read,log_t(t_idx)
      if (ierr /= 0) then
         ierr = 1
         return
      end if
      if (t_idx == 1) num_rho_points=rho_count_read
!     Check limits and number of density points of table
      if ( num_t_points > max_t_points .or.  num_rho_points > max_rho_points ) then
         ierr = 1
         return
      end if
      if ( t_idx > 1 .and. rho_count_read /= num_rho_points ) then
         ierr = 1
         return
      end if
      do rho_idx = 1, num_rho_points
      read(file_unit, iostat=ierr) (table_data(t_idx,rho_idx,var_idx),var_idx=1,num_vars)
      if (ierr /= 0) then
         ierr = 1
         return
      end if
      end do
      end do
end subroutine rtab

end module mhd_table_reader
