!----------------------------------------------------------------------
! allard_atmosphere
!----------------------------------------------------------------------
! The surface boundary condition from the Allard model atmospheres
! (NextGen, AMES-Cond, AMES-Dusty): log P, and with allard_use_tau100
! log P and log T at tau = 100, at (log Teff, log g) for one [Fe/H] and
! [alpha/Fe].
!     read_allard_tables          read a table file (the old NextGen 1
!                                 layout or the newer one) at atm_init,
!                                 then find_allard_table_bounds
!     interpolate_allard_surface  log P (and log T) at (log Teff, log g);
!                                 lookup_failed outside the tables
! The tables live in atm_tables (atm_table%allard_*).
module allard_atmosphere
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: atm_settings, atm_table, atm_table_nga, &
           atm_table_nta
      use io_units, only: run_log_unit
      use numerics_lib, only: locate, polint
      use math_lib
      implicit none
      private
      public :: interpolate_allard_surface, read_allard_tables

contains

!----------------------------------------------------------------------
! interpolate_allard_surface (formerly alsurfp)
!----------------------------------------------------------------------
! PURPOSE
! To obtain the desired properties of stellar and sub-stellar atmospheres from atmospheric
! tables developed by France Allard, et al. Given Log(Teff) and GL, it returns Log(P) at
! this Teff and GL, as well as Log(P) and Log(T) at Tau = 100 for this Teff and GL.
!
! Currently supported input formats include the old Allard, et al, circa 1999 NextGen format,
! to which we are adding the newer 2007 NextGen format. As other formats become available, we
! anticipate being able to handle those also.
!
! Results are obtained by doing 4-point Lagrange interpolation in two dimensions, in TeffL and
! GL. This means we will be interpolating (or extrapolating) over sixteen points in four rows
! of four. We require that the desired Teff and GL are within the limits of the extended range
! of the table. The extended range for GL is from one column's width past the lowest GL on the
! left to one column's width past the highest GL on the right for the current row. The extended
! range for TeffL is defined in the same manner as for the columns. Note that the limits are on
! a row by row and column by column basis. A TeffL, GL point outside the table sets
! lookup_failed (or ierr, below the table in TeffL), with a diagnostic.
!
! Arguments of table lookup subroutine
!
! log_teff      - Input. Value of TL at which the results are desired.
! log_g         - Input. Value of GL at which the results are desired.
! print_to_files - Input. Selects which of the two run-log banners (formats
!              70/71 or 73/74) accompanies the interpolated Log(P), Log(T).
! lookup_failed - Output. If .TRUE. the specified TeffL/GL is above the table in
!              TeffL or outside it in GL and no result was stored; the caller
!              (atm_lib's atm_get_surface_pt) then reports the failure.
! ierr          - Output. 1 when TeffL is below the table minimum or a polint
!              call fails; 0 otherwise.
!
! RESULTS (stored in atm_table, atm_tables.f90)
!
! atm_log10_pressure    - Log(Pressure) for the input TeffL and GL. If
!              atm_settings%allard_use_tau100 is .FALSE. it is PL at TEFFL,GL;
!              if .TRUE. it is PL at TAU=100 for the starting TEFFL,GL.
! atm_log10_temperature - Log(Temperature): TEFFL itself, or TL at TAU=100 when
!              atm_settings%allard_use_tau100 is .TRUE.
!
! Internal Arrays
! TEFFLs      1-D Array of T=TeffL table elements - one for each TeffL. (Columns)
! GLs            1-D Array of G = GL table elements - one for each GL. (Rows)
! PLs            2-D array of assoiated PL's for the matching T in TEFFLs and the matching G in GLA.
!            - One for each TeffL by One for each GL
! P100Ls      2-D array of assoiated P100L's for the matching T in TEFFLs and  the matching G in GLA.
!             - One for each TeffL by One for each GL
! T100Ls      2-D array of assoiated T100L's for the matching T in TEFFLs and  the matching G in GLA.
!             - One for each TeffL by One for each GL
!
! In addition, we define the following auxiliarry arrays
!
! For each row (or GL):
! GLmin      Minimum acceptable GL (including extrapolation) for this row.   - One for each TeffL
! GLmax      Maximum acceptable GL (including extrapolation) for this row.   - One for each TeffL
! IGLmin      Row index of first acceptable GL in table (excludes extrapolation) - One for each TeffL
! IGLmax      Row index of first acceptable GL in table (excludes extrapolation) - One for each TeffL
!
! For teh Table
! TEFFLmin            Minimum acceptable TeffL (Including extrapolation) for the table
! TEFFLmax            Maximum acceptable TeffL (including extrapolation) for the table.
!
! OVERALL LOGIC
!
! I      Determine if TeffL and GL are within overall limits of the table. If TEFFL is greater than the
!     max of all table TL's or less than the minimum, or GL is of the same nature, the program has
!     failed.  To fail, write diagnostic message and go to error exit.
!
! II      Locate TL and GL in table. This gives us a column index, iGL, and a row index, iTL. The subtable
!     used for Lagrange interpolation will        have four rows of four elements. Find indices such that
!     each row is in the table.  Start with TeffL. If it is out of the table, adjust the TeffL index
!     such that all four rows are in the table. Then for each row asjust the GL index such that the
!     row is in the table. If this is not possible, the program has failed. To fail, write diagnostic
!     message and go to error exit.
!
! III      Next we do four-point-Lagrange-interpolation in GL in the first row, then the second, then the
!     third, then the fourth rows. We obtain, for each row, the associated Log(P), Log(P100) and
!     Log(T100) for each row's specific TL.
!
!  IV      To get our desired results, we interpolate the results from III. PL is the result of 4-point
!     Lagrange interpolation in TL using the four Log(P)'s obtained in III above. Likewist P100L and
!     T100L.
!
!       If LPRT  is .TRUE. print Log(P) at the associated Log(Teff) to ISHORT and IMODPT.
subroutine interpolate_allard_surface(log_teff, log_g, print_to_files, lookup_failed, ierr)


      real(dp), intent(in) :: log_teff, log_g
      logical, intent(in) :: print_to_files
      logical, intent(out) :: lookup_failed

      integer :: gl_index(4)
! sized 20 (only the first 4 elements are ever used) to match polint's
! dummy arguments xa(20)/ya(20) -- see numerics/polint.f90's header
! comment and rotation/shape/rotation_shape_factors.f90's identical convention.
      real(dp) :: pressure_row(20), pressure_tau100_row(20), &
           temp_tau100_row(20), pressure_col(20), pressure_tau100_col(20), &
           temp_tau100_col(20)
      logical :: bad_point
      integer :: teffl_index, i, i1, j, j1, j2, j3, n
      real(dp) :: unused_dy

!       Determine if TeffL and GL are within overall limits of the table. If TEFFL is greater than the
!    max of all table TL's or GL is greater than the max or less than the min. Under these
!    circumstances we want to switch to another kind of atmosphere by changing KTTAU. Most likely
!    the calling routine should set KTTAU to 0, to specify a gray atmosphere.

      integer, intent(out) :: ierr

      ierr = 0

      bad_point = .false.       ! presume we have started at a good point

      if (log_teff >= atm_table%allard_teffl_max) then
         bad_point = .true.
         write(run_log_unit,*)
         write(run_log_unit,*)'ALSURFP: TEFFL greater than TeffLmax: '
         write(run_log_unit,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               atm_table%allard_teffl_min,log_teff,atm_table%allard_teffl_max
         write(*,*)
         write(*,*)'ALSURFP: TEFFL greater than TeffLmax: '
         write(*,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               atm_table%allard_teffl_min,log_teff,atm_table%allard_teffl_max
      endif


      if ((log_g >= atm_table%allard_gl_max) .or. (log_g <= atm_table%allard_gl_min)) then
         bad_point = .true.
         write(run_log_unit,*)
         write(*,*)
         write(run_log_unit,*)'ALSURFP: GL out of max range: '
         write(run_log_unit,*)'    GLXmin,GL,GLXmax; ', &
               atm_table%allard_gl_min,log_g,atm_table%allard_gl_max
         write(*,*)'ALSURFP: GL out of max range: '
         write(*,*)'    GLXmin,GL,GLXmax; ', &
               atm_table%allard_gl_min,log_g,atm_table%allard_gl_max
      endif

      if (bad_point) then
         lookup_failed = .true.  ! We are above the table in Teffl, or outside the table in GL
                           ! Set LAlFail so calling routine can take appropiate action
           return
        else
         lookup_failed = .false.  ! We are not above the table in Teffl, and inside the table in GL
      endif

!    If TeffL is less than TeffLmin, the program has failed.  We write a diagnostic message and go
!    to the error exit.

      if (log_teff <= atm_table%allard_teffl_min) then
         write(run_log_unit,*)
         write(run_log_unit,*)'ALSURFP: TEFFL less than TEFFLmin: '
         write(run_log_unit,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               atm_table%allard_teffl_min,log_teff,atm_table%allard_teffl_max
         write(*,*)
         write(*,*)'ALSURFP: TEFFL less than TEFFLmin: '
        write(*,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               atm_table%allard_teffl_min,log_teff,atm_table%allard_teffl_max
!       End program at Error Exit (was go to 9999)
        write(*,*)
        write(*,*)'******** ALSURFP: Program Terminated ********'
        write(*,*)
        write(run_log_unit,*)
        write(run_log_unit,*)'******** ALSURF: Program Terminated ********'
        write(run_log_unit,*)
        ierr = 1
        return
      endif

! II      Locate TEFFL and GL in table. This gives us a column index, iGL, and a row index, iTEFFL. The
!     subtable used for Lagrange interpolation will have four rows of four elements. Find indices such
!     that each row is in the table.  Start with TeffL. If it is out of the table, adjust the TeffL
!     index such that all four rows are in the table. Then for each row asjust the GL index such that the
!     row is in the table. If this is not possible, the program has failed. To fail, write diagnostic
!     message and go to error exit.

!     First locate TEFFL
      call locate(atm_table%allard_teffl_grid,atm_table%allard_num_teff,log_teff,teffl_index) ! returns row 2 in Teffs for 4x4 Lagrange
      teffl_index = teffl_index-1  ! Change to row 1 of 4x4 square in TEFFLs
      if (teffl_index < 1) teffl_index = 1                ! Be sure it is safe to do 4-point Lagrange interpolation
      if (teffl_index > (atm_table%allard_num_teff - 3)) teffl_index = atm_table%allard_num_teff - 3

!     Now find the right index for each row.
      do i = 1, 4      ! Start with iTEFFL and go up to iTEFFL=3
         i1 =  teffl_index+i-1      ! current row pointer - index in TEFFLs aray (0 to 3)
         j1 = atm_table%allard_gl_index_min(i1)         ! lowest valid entry in this row
         j2 = atm_table%allard_gl_index_max(i1)         ! higest valid entry in this row
         n = j2-j1+1          ! number of OK columns
         call locate(atm_table%allard_gl_grid(j1),n,log_g,j3)  ! Returns offset in GLs such that GL is >= GLs(jl1+j3-1)
                        ! and GL,GLs(JL1+j3)
         gl_index(i) = j1+j3-1      ! Save the offset
         gl_index(i) = gl_index(i)-1  ! Change from col 2 to column 1 in GLs direction of 4x4 square
!        FIrst make sure it is safe to do 4-point Lagrange interpolaltion
         if (gl_index(i) < j1) gl_index(i) = j1
         if (gl_index(i) > (j2 - 3)) gl_index(i) = j2 - 3

!        Next, make sue that GL is in range in the current row.  If not, set out of table flag
         if ((log_g < atm_table%allard_gl_row_min(i1)) .or. (log_g > atm_table%allard_gl_row_max(i1))) then      ! and exit
            write(run_log_unit,*)
            write(run_log_unit,*)'ALSURFP: GL out of extended range: '
            write(run_log_unit,*)'    GLmin,GL,GLmax; ', &
               atm_table%allard_gl_row_min(i1),log_g,atm_table%allard_gl_row_max(i1)
            write(*,*)
            write(*,*)'ALSURFP: GL out of extended range: '
           write(*,*)'    GLmin,GL,GLmax; ', &
               atm_table%allard_gl_row_min(i1),log_g,atm_table%allard_gl_row_max(i1)
           lookup_failed = .true.
           return
        endif
      enddo


! III      Next we do four-point-Lagrange-interpolation in GL in the first row, then the second, then the
!     third, then the fourth rows. We obtain, for each row, the associated Log(P), Log(P100) and
!     Log(T100) for each row's specific TL.

!     Do the Lagrange interpolation in GL
      do i = 1, 4      ! Start with iTEFFL and go up to iTEFFL=3, a row at a time
         i1 = teffl_index + i -1   ! Index in TEFFL array (0 to 3)
         do j= 1, 4
            j1 = gl_index(i) + j -1  ! Index in GL array - from iGL(i) to iGL(i)+3
            pressure_row(j) = atm_table%allard_log10_pressure(i1,j1)        ! Get this row's four working PL's
            pressure_tau100_row(j) = atm_table%allard_log10_pressure_tau100(i1,j1)  ! Get this row's four working P100L's
            temp_tau100_row(j) = atm_table%allard_log10_temp_tau100(i1,j1)  ! Get this row's four working T100L's
         enddo
         j1 = gl_index(i)
         call polint(atm_table%allard_gl_grid(j1),pressure_row,4,log_g,pressure_col(i),unused_dy, ierr)
         if (ierr /= 0) return
         call polint(atm_table%allard_gl_grid(j1),pressure_tau100_row,4,log_g,pressure_tau100_col(i),unused_dy, ierr)
         if (ierr /= 0) return
         call polint(atm_table%allard_gl_grid(j1),temp_tau100_row,4,log_g,temp_tau100_col(i),unused_dy, ierr)
         if (ierr /= 0) return
      enddo
!     Now do final 4-point Lagrange inerpolations in TEFFL
      if (.not. atm_settings%allard_use_tau100) then
!        Current standard alternative. PL for TEFFL,GL, TL=TEFFL
         call polint(atm_table%allard_teffl_grid(teffl_index),pressure_col,4,log_teff,atm_table%atm_log10_pressure,unused_dy, ierr)
         if (ierr /= 0) return
         atm_table%atm_log10_temperature = log_teff
      else
         call polint(atm_table%allard_teffl_grid(teffl_index),pressure_tau100_col,4,log_teff,atm_table%atm_log10_pressure,unused_dy, ierr)
         if (ierr /= 0) return
         call polint(atm_table%allard_teffl_grid(teffl_index),temp_tau100_col,4,log_teff,atm_table%atm_log10_temperature,unused_dy, ierr)
         if (ierr /= 0) return
      endif

!     We now have obtained the needed temparatures and pressures

!     If requested, WRITE OUT INFORMATION TO THE MODEL FILE.
      if (print_to_files) then
         if (atm_settings%allard_use_tau100) then
            write(run_log_unit,70)
   70       format('********PL,TL at Tau=100 INTERPOLATED' &
                 , ' FROM ALLARD TABULATED VALUES********')
         else
            write(run_log_unit,71)
   71       format('********PL,TL at T=TEFF,GL INTERPOLATED ' &
                , 'FROM ALLARD TABULATED VALUES********')
         endif
         write(run_log_unit,72) log_teff,log_g,atm_table%atm_log10_temperature,atm_table%atm_log10_pressure
   72    format(' ',20x,'LOG(Teffl) =',f10.5,' ,LOG(G) =' ,f10.5, &
              ', LOG(T) =',f10.5,', Log(P) =', f10.5)
      endif

!     If not requested via LPRT, WRITE OUT INFORMATION TO THE SHORT FILE.
      if (.not. print_to_files) then
        if (atm_settings%allard_use_tau100) then
           write(run_log_unit,73)
   73      format('ALSURFP:  PL,TL at Tau=100 INTERPOLATED' &
                , ' FROM ALLARD TABULATED VALUES:')
        else
           write(run_log_unit,74)
   74      format('ALSURFP: PL,TL at T=TEFF,GL INTERPOLATED ' &
                , 'FROM ALLARD TABULATED VALUES:')
        endif
        write(run_log_unit,75) log_teff,log_g,atm_table%atm_log10_temperature,atm_table%atm_log10_pressure
   75   format(' ',10x,'LOG(Teffl) =',f10.5,', LOG(G) =',f10.5, &
             ', LOG(T) =',f10.5,', Log(P) =',f10.5)
      endif


end subroutine interpolate_allard_surface

!----------------------------------------------------------------------
! read_allard_tables (formerly alfilein)
!----------------------------------------------------------------------
! PURPOSE
! To read Allard atmosphere information files into the appropriate input arrays
! and calculate associated auxilliary information. These files contain the
! desired properties of stellar and substellar atmospheres from atmospheric
! tables developed by France Allard, et al. These properties include Log(P)
! at Teff and GL, as well as Log(P) and Log(T) at Tau = 100 for this Teff and
! GL. Currently supported input formats include the old Allard, et al, circa
! 1999 NextGen format, to which we are adding the newer 2007 NextGen format.
! As other formats become available, weanticipate being able to handle those also.
!
! Input files may be in Old NextGen file format or in the new Allard Atmosphere
! file format. This format, is of a conceptually different style from the old.
! The data part of the new file format has one record for each Teff, GL, FE/h,
! Alpha. The rest of that record contains the associated PL @ T=Teff,
! PL @ Tau=100 and TL @ Tau100. In other words the new format has one data record
! for each item in the table. The old format had one record for each GL in the
! table. This new format is more flexible in supporting future needs.
!
! The contents of these files are put into the arrays and elements used by the
! Allard Atmosphere table lookup routine, interpolate_allard_surface. In addition, we fill any
! vacant entries in the table with -999.0. A value of -999. in the table
! identifies the entry as "invalid." In addition, key auxililiary arrays, as
! described below, are created.
subroutine read_allard_tables(allard_table_path, ierr)
!
! Parameters NTA and NGA are respectively the maximum expected numbers
! of Teff's and GL's we expect to encounter, even in future tables
! The maximum nuber of Teff's in the current table is nTeff, and the
! associated max number of GL's is nGL.


      character(len=*), intent(in) :: allard_table_path
      character(len=35) :: first_record
      character(len=256) :: header_line
      real(dp) :: local_teffs(atm_table_nta)

      integer :: i, j, nhdr, irecno, i1, j1
      real(dp) :: teff_value, gl_value, feh_value, alpha_value
      real(dp) :: pressure_value, pressure_tau100_value, temp_tau100_value
! latmtptau100 is never set, so the check below that an old NextGen 1
! file cannot give the tau = 100 point never fires: in the original it
! was a misspelling of the allard_use_tau100 flag (LALTPTau100), kept
! as it was (an open item in docs/ROADMAP.md).
      logical :: latmtptau100

      integer, intent(out) :: ierr
      integer :: ios, allard_unit
      logical :: found_first

      ierr = 0

!     set output arrays to invalid values
      do i = 1,atm_table_nta
         atm_table%allard_teffl_grid(i) = -999.d0
         do j = 1,atm_table_nga
            atm_table%allard_gl_grid(j) = -999.d0
            atm_table%allard_log10_pressure(i,j) = -999.d0
            atm_table%allard_log10_pressure_tau100(i,j) = -999.d0
            atm_table%allard_log10_temp_tau100(i,j) = -999.d0
         enddo
      enddo
!     the input file name is in FALLARD and its unit number is IOATMA
      open(newunit=allard_unit,file=allard_table_path,status='OLD',iostat=ios)
      if (ios /= 0) then
       write(*,*)
       write(*,*) 'ALFILEIN: Allard File OPEN Failure'
       write(*,*) '     File name: ', trim(allard_table_path)
       call fail_9999
       ierr = 1
       return
      end if

!     Decide what kind of input file
      read(allard_unit,901) first_record
  901  format(a)
      rewind(allard_unit)
      if (first_record(4:7) == '3.50') then

!     Process old-stype nextgen input file.
       atm_table%allard_is_old_nextgen = .true.
       write(run_log_unit,*) 'ALFileIn: File Description: 1999 NEXTGEN', &
           ' (Old NEXTGEN)'

!     Ensure that we are not requesting PL,TL at Tau-100. These are not
!     present in the old Nextgen 1 files. If requested, fail
      if (latmtptau100) then
         write(*,*)
         write(*,*) 'ALFileIN: Invalid old Allard input file ', &
               'for requested PT,TL at Tau=100'
         write(run_log_unit,*)
         write(run_log_unit,*) 'ALFileIN: Invalid old Allard input file ', &
               'for requested PT,TL at Tau=100'
         call fail_9999
         ierr = 1
         return
      endif
!     Read range of gravities
      atm_table%allard_num_teff = 54
      atm_table%allard_num_gl = 5
      read(allard_unit,905) (atm_table%allard_gl_grid(i),i=1,atm_table%allard_num_gl)
  905  format(5f7.2/)
!       read(ioatma,*)

!     Read range of temperatures
      read(allard_unit,907) (atm_table%allard_teffl_grid(i),i=1,atm_table%allard_num_teff)
  907  format(10(1p5e16.8,/),1p4e16.8,/)
!       read(ioatma,*)
      read(allard_unit,909) ((atm_table%allard_log10_pressure(j,i),i=1,atm_table%allard_num_gl),j=1,atm_table%allard_num_teff)
  909  format(1p5e16.8)
      close(allard_unit)
!     Old Nextgen input file has been read an put into the proper arrays

      else if (first_record(4:9) == 'ALLARD') then

!     Process new-style Allard Atmosphere file
       atm_table%allard_is_old_nextgen = .false.
      read(allard_unit,911) nhdr, header_line  ! nhdr is number of header recoreds to skip over
  911      format(i2,x,a)
       write(run_log_unit,912) 'ALFilein: New Allard Atm: File Description: ',header_line(1:47)
  912  format(2a)
!     Skip over rest of header, if any
      if (nhdr > 1) then
        do i = 2, nhdr
          read(allard_unit,*)
        enddo
      endif
      atm_table%allard_num_gl = 0
      atm_table%allard_num_teff = 0
      atm_table%allard_num_feh = 0
      atm_table%allard_num_alpha = 0
      irecno = 0

!     First we read the file to count the # of Teff,GL,Fe/H and Alpha's

      found_first = .false.
      do
      read(allard_unit,915,iostat=ios) teff_value,gl_value,feh_value,alpha_value
  915      format(f6.0,3f6.2)
      if (ios > 0) then
!     File Read error exit
       write(*,*)
       write(*,*)'AtTabInit Pass 1 File Error '
       write(*,*)
       call fail_9999
       ierr = 1
       return
      end if
      if (ios < 0) exit
      irecno=irecno+1

!     If the record does not have the correct FeH and Alpha, we skip it
      if ((abs(feh_value-atm_settings%allard_target_feh) > 1d-5) .or. &
             (abs(alpha_value-atm_settings%allard_target_alpha) > 0d0)) then
         cycle  ! On to next record. Stay in this part until we get an acceptable record
      endif
      if (.not. found_first) then
      found_first = .true.
      atm_table%allard_num_teff = 1    ! Initialize counters and variables after reading first record
      atm_table%allard_num_gl = 1
      atm_table%allard_num_feh = 1
      atm_table%allard_num_alpha = 1
      local_teffs(atm_table%allard_num_teff) = teff_value
      atm_table%allard_gl_grid(atm_table%allard_num_gl) = gl_value
      atm_table%allard_feh_grid(atm_table%allard_num_feh) = feh_value
      atm_table%allard_alpha_grid(atm_table%allard_num_alpha) = alpha_value
      cycle
      end if

!      Now check Teffs, increase counter if needed
       do i = 1,atm_table%allard_num_teff  ! Skip out if any old Teff is a match
        if (abs(teff_value-local_teffs(i)) < 1d-6) exit
      enddo
      if (i > atm_table%allard_num_teff) then
        atm_table%allard_num_teff = atm_table%allard_num_teff+1      ! count the nextTeff's
        local_teffs(atm_table%allard_num_teff) = teff_value  ! save the new Teff in array Teffs
      end if

!     Now check GLs, increase counter if needed
      do i = 1,atm_table%allard_num_gl  ! Skip out if any old GL is a match
        if (abs(gl_value-atm_table%allard_gl_grid(i)) < 1d-6) exit
      enddo
      if (i > atm_table%allard_num_gl) then
        atm_table%allard_num_gl = atm_table%allard_num_gl+1      ! count the next GL
        atm_table%allard_gl_grid(atm_table%allard_num_gl) = gl_value  ! save the new GL in array GLs
      end if

!     Now check FeHs, increase counter if needed
!     (NOTE: these two scans are bounded by allard_num_gl in the
!     original; preserved as-is.)
      do i = 1,atm_table%allard_num_gl  ! Skip out if any old FeH is a match
        if (abs(feh_value-atm_table%allard_feh_grid(i)) < 1d-6) exit
      enddo
      if (i > atm_table%allard_num_gl) then
        atm_table%allard_num_feh = atm_table%allard_num_feh+1      ! count the next FeH
        atm_table%allard_feh_grid(atm_table%allard_num_feh) = feh_value  ! save the new FeH in array FeHs
      end if

!     Now check ALPHAs, increase counter if needed
      do i = 1,atm_table%allard_num_gl  ! Skip out if any old ALPHA is a match
        if (abs(alpha_value-atm_table%allard_alpha_grid(i)) < 1d-6) exit
      enddo
      if (i > atm_table%allard_num_gl) then
        atm_table%allard_num_alpha = atm_table%allard_num_alpha+1      ! count the next ALPHA
        atm_table%allard_alpha_grid(atm_table%allard_num_alpha) = alpha_value  ! save the new ALPHA in array ALPHAs
      end if

      end do ! go on to next record

!       write(*,*) 'nTeff,nGL,nFeH,nALPHA: ',nTeff,nGL,nFeH,nALPHA

!      Now we have unsorted Teff,GL,Fe/H and Alpha's.
!     The next step is to sort them
      call sort_shell(atm_table%allard_num_teff,local_teffs)
      call sort_shell(atm_table%allard_num_gl,atm_table%allard_gl_grid)
      call sort_shell(atm_table%allard_num_feh,atm_table%allard_feh_grid)
      call sort_shell(atm_table%allard_num_alpha,atm_table%allard_alpha_grid)

!     Now we can rewind the input file and start over for real
      rewind(allard_unit)  ! go back to the beginning of the file

!      Skip over header
      do i = 1,nhdr
         read(allard_unit,*)
      enddo

      do
      read(allard_unit,920,iostat=ios) teff_value,gl_value,feh_value,alpha_value, &
              pressure_value,pressure_tau100_value,temp_tau100_value
  920      format(f6.0,3f6.2,1p4d16.8)
      if (ios > 0) then
       write(*,*)
       write(*,*)'AtTabInit File Read error'
       write(*,*)
       call fail_9999
       ierr = 1
       return
      end if
      if (ios < 0) exit  ! We heve finished with the input file and entered all inputs

!     If the record does not have the correct FeH and Alpha, we skip it
      if ((abs(feh_value-atm_settings%allard_target_feh) > 1d-5) .or. &
             (abs(alpha_value-atm_settings%allard_target_alpha) > 0d0)) then
         cycle  ! On to next record
      endif

      do i = 1,atm_table%allard_num_teff  ! Skip out when we find the matching Teff
        if (abs(teff_value-local_teffs(i)) < 1d-6) then
           i1 = i
           exit
        endif
      enddo
      if (i > atm_table%allard_num_teff) then
      write(*,*) 'ALFilein: Impossible failure #1'
      write(*,*) 'Teff of ', teff_value, ' not in Teff table.'
      call fail_9999
      ierr = 1
      return
      end if

      do j = 1,atm_table%allard_num_gl  ! Skip out when we find the a matching GL
        if (abs(gl_value-atm_table%allard_gl_grid(j)) < 1d-6) then
           j1 = j
           exit
        endif
      enddo
      if (j > atm_table%allard_num_gl) then
      write(*,*) 'ALFilein: Impossible failure #2'
      write(*,*) 'GL of ', gl_value, ' not in GL table.'
      call fail_9999
      ierr = 1
      return
      end if

!     we now verify that we have the correct FeH and Alpha

      if (abs(feh_value -atm_settings%allard_target_feh) < 1d-6) then
      if (abs(alpha_value -atm_settings%allard_target_alpha) < 1d-6) then


!     We now have the correct indices for our tables, i1 for the Teff-direction
!     and j1 for the GL direction

       atm_table%allard_teffl_grid(i1) = log10(teff_value)  ! We need log(teff), but up to now we've been using Teff
       atm_table%allard_gl_grid(j1) = gl_value
       atm_table%allard_log10_pressure(i1,j1) = pressure_value
       atm_table%allard_log10_pressure_tau100(i1,j1) = pressure_tau100_value
       atm_table%allard_log10_temp_tau100(i1,j1) = temp_tau100_value

      end if
      end if
      end do   ! Back to process the next record
      close(allard_unit)

      else
!     If we get here, the input file is invalid
      write(run_log_unit,*)'*** Invalid Allard Atmosphere input file ***'
      write(run_log_unit,*) "First record was '",first_record,"'"
      write(*,*)'*** Invalid Allard Atmosphere input file ***'
      write(*,*) "First record was '",first_record,"'"
      call fail_9999
      ierr = 1
      return
      end if

      call find_allard_table_bounds(ierr)   ! Initialize Allard tables
      if (ierr /= 0) return

       return

contains

! THE ERROR EXIT (was label 9999): announce termination; every caller
! then sets ierr = 1 and returns.
subroutine fail_9999
       write(*,*)
       write(*,*)
       write(*,*) '**************** PROGRAM ALFilein TERMINATED ', &
             '@ 9999 **************'
       write(*,*)
       write(*,*)
       write(*,*)
       write(run_log_unit,*)
       write(run_log_unit,*) '**************** PROGRAM ALFilein TERMINATED ', &
            '@ 9999 *************'
       write(run_log_unit,*)
       write(run_log_unit,*)
end subroutine fail_9999

end subroutine read_allard_tables      ! END OF ALINITTAB

! ************* End of read_allard_tables *************************


! Numerical Recipes Shell sort

pure subroutine sort_shell(num_elements, values)

      integer, intent(in) :: num_elements
      real(dp), intent(inout) :: values(:)
! Sorts an array A(1:nn) into ascending numerical order by Shell's method (diminishing
! incremental sort) nn is input. A is replaced on output by its sorted arrangement
      integer :: n, inc, i, j
      real(dp) :: v

      if (num_elements == 1) return                  ! Exit if only one element
      n=num_elements
      inc=1
      do
         inc=3*inc+1                         ! Determin starting increment
      if (inc > n) exit
      end do
      do
         inc=inc/3                        ! Loop over the partial sorts
         do i=inc+1,n
            v=values(i)                        ! Outer loop of straight insertion
            j=i
    do while (values(j-inc) > v)
               values(j)=values(j-inc)
               j=j-inc
               if (j <= inc) exit
    end do
            values(j)=v
         enddo
      if (inc <= 1) exit
      end do
end subroutine sort_shell

! ************* End of sort_shell *************************

!----------------------------------------------------------------------
! find_allard_table_bounds (formerly altabinit)
!----------------------------------------------------------------------
! PURPOSE
! To analyze the internal Allard-format tables provided by read_allard_tables,
! verifying them and creating needed additional tables and elements.
!
! Several steps are taken:
! 1. Find the minimum (TEFFLmin) and maximum (TEFFLmax) permissable values
!    of TEFFL. These are one row's width below the bottom and above the
!    top of the table. Because the first level of interpolation is in GL,
!    only a single minimum and maximum value of TEFFL are needed.
! 2. For each row in GL, find the index of the smallest element (iGLmin)
!    and the index of the largest element (iGLmax). There is one pair of
!    these for each TEFFL.
! 3. For each row in GL, find the the minimum (GLmin) and maximum (GLmax)
!    permissable value of GL. These are one column's width less than the
!    row minimum in GL and one column's width greater than the column
!    maximum. There is one pair of these for each TEFFL.
! 4. Validate the table. (a) Ensure that there are no invalid elements
!    inside the table, i.e., no invalid element between the row's iGLmin
!    and iGLmax. This check is made for every TEFFL. (b) Ensure that
!    every row has at least 4 valid entries and that there are at least
!    4 columns.
subroutine find_allard_table_bounds(ierr)


      integer, intent(out) :: ierr
      integer :: i, j, j1, j2
      logical :: table_is_bad

!     1. Find the minimum (TEFFLmin) and maximum (TEFFLmax) permissable values of TEFFL. These are
!        one row's width below the bottom and above the top of the table. Because the first level
!        of interpolation is in GL, only a single minimum and maximum value of TEFFL are needed.

      ierr = 0

      atm_table%allard_al_teffl_min = atm_table%allard_teffl_grid(1)
      atm_table%allard_al_teffl_max = atm_table%allard_teffl_grid(atm_table%allard_num_teff)
      atm_table%allard_teffl_min = atm_table%allard_teffl_grid(1) - (atm_table%allard_teffl_grid(2)-atm_table%allard_teffl_grid(1))
      atm_table%allard_teffl_max = atm_table%allard_teffl_grid(atm_table%allard_num_teff) + &
           (atm_table%allard_teffl_grid(atm_table%allard_num_teff)-atm_table%allard_teffl_grid(atm_table%allard_num_teff-1))

!     2. For each row in GL, find the index of the smallest element (iGLmin) and the index of the
!        largest element (iGLmax). There is one pair of these for each TEFFL.
      do i = 1, atm_table%allard_num_teff
         do j = 1, atm_table%allard_num_gl
            if (atm_table%allard_log10_pressure(i,j) > -998d0) then
               atm_table%allard_gl_index_min(i) = j
               exit
            endif
         enddo
         do j = atm_table%allard_num_gl, 1, -1
            if (atm_table%allard_log10_pressure(i,j) > -998d0) then
               atm_table%allard_gl_index_max(i) = j
               exit
            endif
         enddo

      enddo

!     3. For each row in GL, find the the minimum (GLmin) and maximum (GLmax) permissable value of
!        GL. These are one column's width less than the row minimum in GL and one column's width
!        greater than the column maximum. There is one pair of these for each TEFFL.

      atm_table%allard_gl_min = 999d0
      atm_table%allard_gl_max = -999d0
      do i = 1, atm_table%allard_num_teff
         j1 = atm_table%allard_gl_index_min(i)
         j2 = atm_table%allard_gl_index_max(i)
         atm_table%allard_gl_row_min(i) = atm_table%allard_gl_grid(j1) - 4d0*(atm_table%allard_gl_grid(j1+1) - atm_table%allard_gl_grid(j1))
         if (atm_table%allard_gl_row_min(i) < atm_table%allard_gl_min) atm_table%allard_gl_min = atm_table%allard_gl_row_min(i)
         atm_table%allard_gl_row_max(i) = atm_table%allard_gl_grid(j2) + (atm_table%allard_gl_grid(j2) - atm_table%allard_gl_grid(j2-1))
         if (atm_table%allard_gl_row_max(i) > atm_table%allard_gl_max) atm_table%allard_gl_max = atm_table%allard_gl_row_max(i)
      enddo

!     4. Validate the table. (a) Ensure that there are no invalid elements inside the table, i.e.,
!        no invalid element between the row's iGLmin and iGLmax. This check is made for every TEFFL.

      table_is_bad = .false.
      do i = 1, atm_table%allard_num_teff
         j1 = atm_table%allard_gl_index_min(i)
         j2 = atm_table%allard_gl_index_max(i)
         do j = j1, j2
            if (atm_table%allard_log10_pressure(i,j) < -998d0) then
              table_is_bad = .true.
              write(run_log_unit,900) 'ALTABINIT: Bad input Allard Table: ', &
                 'TEFF, GL: ', pow(10d0, atm_table%allard_teffl_grid(i)), atm_table%allard_gl_grid(j)
  900              format(2a,f5.0,f7.2)
            endif
         enddo
      enddo

!        (b) Ensure that every row has at least 4 valid entries and that there are at least 4 columns.
      if (atm_table%allard_num_teff < 4) then
         table_is_bad = .true.
         write(run_log_unit,910) 'ALTABINIT: Bad input Allard Table: ', &
            'Less than 4 rows: nTeff = ',atm_table%allard_num_teff
  910         format(a,i3)
      endif
      do i = 1, atm_table%allard_num_teff
         j1 = atm_table%allard_gl_index_min(i)
         j2 = atm_table%allard_gl_index_max(i)
         if ((j2 - j1 + 1) < 4) then
            table_is_bad = .true.
            write(run_log_unit,920) 'ALTABINIT: Bad input Allard Table: ', &
              'Row with less that 4 elements: i,#,Teff,GLMin,GLmax: ', &
               i,j2-j1+1,pow(10d0, atm_table%allard_teffl_grid(i)),atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i)
  920            format(2a,2i4,2x,3f7.2)
         endif
      enddo

      if (.not. table_is_bad) then

      return                  ! If good table, return

      end if
       write(*,*)
       write(*,*)'******** ALTABINIT: Program Terminated ********'
       write(*,*)
       write(run_log_unit,*)
       write(run_log_unit,*)'******** ALTABINIT: Program Terminated ********'
       write(run_log_unit,*)
       call print_allard_tables
       ierr = 1

end subroutine find_allard_table_bounds

!----------------------------------------------------------------------
! print_allard_tables
!----------------------------------------------------------------------
! Prints the internal Allard-format atmosphere tables and auxiliary
! information provided by read_allard_tables and verified by find_allard_table_bounds, to the
! run log (run_log_unit).
subroutine print_allard_tables

      integer :: i, j

      atm_table%allard_teffl_min = atm_table%allard_teffl_grid(1) - &
           (atm_table%allard_teffl_grid(2)-atm_table%allard_teffl_grid(1))
      atm_table%allard_teffl_max = atm_table%allard_teffl_grid(atm_table%allard_num_teff) + &
           (atm_table%allard_teffl_grid(atm_table%allard_num_teff)-atm_table%allard_teffl_grid(atm_table%allard_num_teff-1))

      write(run_log_unit,*)
      write(run_log_unit,*)' **************** Allard Table  ***************'
      write(run_log_unit,*)
      write(run_log_unit,*)'nTeff= ',atm_table%allard_num_teff,' ,nGL=',atm_table%allard_num_gl
      write(run_log_unit,*)'nFeH= ',atm_table%allard_num_feh,' ,nAlpha=',atm_table%allard_num_alpha
      write(run_log_unit,*)
      write(run_log_unit,*)'TEFFLmin= ',atm_table%allard_teffl_min,' ,TEFFLmax= ',atm_table%allard_teffl_max
      write(run_log_unit,*)'GLXmin=',atm_table%allard_gl_min,' GLXmax=',atm_table%allard_gl_max
      write(run_log_unit,900)' PL @ Teff','GL ->',(atm_table%allard_gl_grid(j),j=1,atm_table%allard_num_gl)
  900      format(A10,16X,A,20(5X,F6.2,5X))
      write(run_log_unit,*)
      write(run_log_unit,*) 'PL at Teff,GL'
      write(run_log_unit,*)
      do i = 1, atm_table%allard_num_teff
         write(run_log_unit,910) i,pow(10D0, atm_table%allard_teffl_grid(i)), &
              atm_table%allard_gl_index_min(i),atm_table%allard_gl_index_max(i), &
              atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i), &
              (atm_table%allard_log10_pressure(i,j),j=1,atm_table%allard_num_gl)
  910     format(I5,F6.0,2I4,2F6.2,1P20E16.8)
      enddo
      write(run_log_unit,*)
      write(run_log_unit,*) 'PL at Tau=100'
      write(run_log_unit,*)
      do i = 1, atm_table%allard_num_teff
         write(run_log_unit,910) i,pow(10D0, atm_table%allard_teffl_grid(i)), &
              atm_table%allard_gl_index_min(i),atm_table%allard_gl_index_max(i), &
              atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i), &
              (atm_table%allard_log10_pressure_tau100(i,j),j=1,atm_table%allard_num_gl)
      enddo
      write(run_log_unit,*)
      write(run_log_unit,*) 'TL at Tau=100'
      write(run_log_unit,*)
      do i = 1, atm_table%allard_num_teff
         write(run_log_unit,910) i,pow(10D0, atm_table%allard_teffl_grid(i)), &
              atm_table%allard_gl_index_min(i),atm_table%allard_gl_index_max(i), &
              atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i), &
              (atm_table%allard_log10_temp_tau100(i,j),j=1,atm_table%allard_num_gl)
      enddo
      write(run_log_unit,*)

end subroutine print_allard_tables

end module allard_atmosphere
