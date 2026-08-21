!----------------------------------------------------------------------
! alsurfp
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alsurfp.f; only variable names, source form, and comment
! style were updated.
!
! PURPOSE
! To obtain the desired properties of stellar and sub-stellar atmospheres from atmospheric
! tables developed by France Allard, et al. Given Log(Teff) and GL, ASURFP returns Log(P) at
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
! a row by row and column by column basis. If the request TeffL, GL point is outside the table,
! the progam is terminated with a definitive error message.
!
! ARGUMENTS OF TABLE LOOKUP SUBROUTINE
!
! TeffL      - Input argument. Value of TL at which the results are desired.
! GL            - Input argument. Value of GL at which the results are desired.
! AP            - Output argument to COMMON /ATMPRT/. Value of Log(Pressure) for the
!              associated input TeffL and GL If LATMTPTau100 is .FALSE. it is PL
!              at TEFFL,GL. If .TRUE. it is PL at TAU=100 for the starting TEFFL,GL.
! AT            - Output argument to COMMON /ATMPRT/. Value of Log(Temperature) for the
!              associated input TeffL and GL If LATMTPTau100 is .FALSE. it is TEFFL.
!              If .TRUE. it is TL at TAU=100 for the starting TEFFL,GL.
! LPRT        Input argument: If .TRUE. print Log(P) at the associated Log(Teff)  to ISHORT and IMODPT.
! LAlFail   Output argument. If .TRUE. the specified TeffL/GL is outside the Allard tables. The
!           calling routine should set KTTAU to 0 to change to a gray atmosphere, or KTTAU to 3 for
!           aKurucz atmosphere.
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
subroutine alsurfp(log_teff, log_g, print_to_files, lookup_failed)

      use locate_mod
      use polint_mod
      implicit none
      integer, parameter :: nta = 250
      integer, parameter :: nga = 25

      double precision, intent(in) :: log_teff, log_g
      logical, intent(in) :: print_to_files
      logical, intent(out) :: lookup_failed

      integer :: ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

! common/atmprt/: only atm_log10_pressure and atm_log10_temperature
! (AP, AT) are used/set here; the remaining members are unused
! placeholders preserving the shared storage layout.
      double precision :: atm_tau, atm_log10_pressure, atm_log10_temperature, &
           atm_log10_density, atm_opacity, atm_ion_fraction(3)
      common/atmprt/atm_tau, atm_log10_pressure, atm_log10_temperature, &
           atm_log10_density, atm_opacity, atm_ion_fraction

! Shared: ALFILEIN, ALTABINIT and ALSURFP
      double precision :: allard_teffl_grid(nta), allard_gl_grid(nga), &
           allard_feh_grid(nga), allard_alpha_grid(nga), &
           allard_log10_pressure(nta,nga), allard_log10_pressure_tau100(nta,nga), &
           allard_log10_temp_tau100(nta,nga)
      logical :: allard_is_old_nextgen
      integer :: allard_num_teff, allard_num_gl, allard_num_feh, allard_num_alpha
      common /alatm01/ allard_teffl_grid, allard_gl_grid, allard_feh_grid, &
           allard_alpha_grid, allard_log10_pressure, allard_log10_pressure_tau100, &
           allard_log10_temp_tau100, allard_is_old_nextgen, allard_num_teff, &
           allard_num_gl, allard_num_feh, allard_num_alpha
! Shared: ALTABINIT and ALSURFP
      double precision :: allard_gl_row_min(nta), allard_gl_row_max(nta)
      integer :: allard_gl_index_min(nta), allard_gl_index_max(nta)
      double precision :: allard_teffl_min, allard_teffl_max, allard_gl_min, &
           allard_gl_max
      common /alatm02/ allard_gl_row_min, allard_gl_row_max, allard_gl_index_min, &
           allard_gl_index_max, allard_teffl_min, allard_teffl_max, &
           allard_gl_min, allard_gl_max
! Shared: ALFILEIN, ALSURFP and PARMIN
      double precision :: allard_target_feh, allard_target_alpha
      logical :: allard_use_tau100
      integer :: allard_table_unit
      common /alatm03/ allard_target_feh, allard_target_alpha, allard_use_tau100, &
           allard_table_unit
! MHP 8/25 Removed character file names from common block
      double precision :: alatm04_placeholder1, alatm04_placeholder2, &
           alatm04_placeholder3, alatm04_placeholder4
      common /alatm04/ alatm04_placeholder1, alatm04_placeholder2, &
           alatm04_placeholder3, alatm04_placeholder4

      integer :: gl_index(4)
! sized 20 (only the first 4 elements are ever used) to match polint's
! dummy arguments xa(20)/ya(20) -- see numerics/polint.f90's header
! comment and rotation/fpft.f90's identical convention.
      double precision :: pressure_row(20), pressure_tau100_row(20), &
           temp_tau100_row(20), pressure_col(20), pressure_tau100_col(20), &
           temp_tau100_col(20)

      save

      logical :: bad_point
      integer :: teffl_index, i, i1, j, j1, j2, j3, n
      double precision :: unused_dy

!       Determine if TeffL and GL are within overall limits of the table. If TEFFL is greater than the
!    max of all table TL's or GL is greater than the max or less than the min. Under these
!    circumstances we want to switch to another kind of atmosphere by changing KTTAU. Most likely
!    the calling routine should set KTTAU to 0, to specify a gray atmosphere.

      bad_point = .false.       ! presume we have started at a good point

      if (log_teff .ge. allard_teffl_max) then
         bad_point = .true.
         write(short_file_unit,*)
         write(short_file_unit,*)'ALSURFP: TEFFL greater than TeffLmax: '
         write(short_file_unit,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               allard_teffl_min,log_teff,allard_teffl_max
         write(*,*)
         write(*,*)'ALSURFP: TEFFL greater than TeffLmax: '
         write(*,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               allard_teffl_min,log_teff,allard_teffl_max
      endif


      if ((log_g .ge. allard_gl_max) .or. (log_g .le. allard_gl_min)) then
         bad_point = .true.
         write(short_file_unit,*)
         write(*,*)
         write(short_file_unit,*)'ALSURFP: GL out of max range: '
         write(short_file_unit,*)'    GLXmin,GL,GLXmax; ', &
               allard_gl_min,log_g,allard_gl_max
         write(*,*)'ALSURFP: GL out of max range: '
         write(*,*)'    GLXmin,GL,GLXmax; ', &
               allard_gl_min,log_g,allard_gl_max
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

      if (log_teff .le. allard_teffl_min) then
         write(short_file_unit,*)
         write(short_file_unit,*)'ALSURFP: TEFFL less than TEFFLmin: '
         write(short_file_unit,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               allard_teffl_min,log_teff,allard_teffl_max
         write(*,*)
         write(*,*)'ALSURFP: TEFFL less than TEFFLmin: '
        write(*,*)'    TEFFLmin,TEFFL,TEFFLmax; ', &
               allard_teffl_min,log_teff,allard_teffl_max
        go to 9999        ! End program at Error Exit
      endif

! II      Locate TEFFL and GL in table. This gives us a column index, iGL, and a row index, iTEFFL. The
!     subtable used for Lagrange interpolation will have four rows of four elements. Find indices such
!     that each row is in the table.  Start with TeffL. If it is out of the table, adjust the TeffL
!     index such that all four rows are in the table. Then for each row asjust the GL index such that the
!     row is in the table. If this is not possible, the program has failed. To fail, write diagnostic
!     message and go to error exit.

!     First locate TEFFL
      call locate(allard_teffl_grid,allard_num_teff,log_teff,teffl_index) ! returns row 2 in Teffs for 4x4 Lagrange
      teffl_index = teffl_index-1  ! Change to row 1 of 4x4 square in TEFFLs
      if (teffl_index .lt. 1) teffl_index = 1                ! Be sure it is safe to do 4-point Lagrange interpolation
      if (teffl_index .gt. (allard_num_teff - 3)) teffl_index = allard_num_teff - 3

!     Now find the right index for each row.
      do i = 1, 4      ! Start with iTEFFL and go up to iTEFFL=3
         i1 =  teffl_index+i-1      ! current row pointer - index in TEFFLs aray (0 to 3)
         j1 = allard_gl_index_min(i1)         ! lowest valid entry in this row
         j2 = allard_gl_index_max(i1)         ! higest valid entry in this row
         n = j2-j1+1          ! number of OK columns
         call locate(allard_gl_grid(j1),n,log_g,j3)  ! Returns offset in GLs such that GL is >= GLs(jl1+j3-1)
                        ! and GL,GLs(JL1+j3)
         gl_index(i) = j1+j3-1      ! Save the offset
         gl_index(i) = gl_index(i)-1  ! Change from col 2 to column 1 in GLs direction of 4x4 square
!        FIrst make sure it is safe to do 4-point Lagrange interpolaltion
         if (gl_index(i) .lt. j1) gl_index(i) = j1
         if (gl_index(i) .gt. (j2 - 3)) gl_index(i) = j2 - 3

!        Next, make sue that GL is in range in the current row.  If not, set out of table flag
         if ((log_g .lt. allard_gl_row_min(i1)) .or. (log_g .gt. allard_gl_row_max(i1))) then      ! and exit
            bad_point = .true.
            write(short_file_unit,*)
            write(short_file_unit,*)'ALSURFP: GL out of extended range: '
            write(short_file_unit,*)'    GLmin,GL,GLmax; ', &
               allard_gl_row_min(i1),log_g,allard_gl_row_max(i1)
            write(*,*)
            write(*,*)'ALSURFP: GL out of extended range: '
           write(*,*)'    GLmin,GL,GLmax; ', &
               allard_gl_row_min(i1),log_g,allard_gl_row_max(i1)
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
            pressure_row(j) = allard_log10_pressure(i1,j1)        ! Get this row's four working PL's
            pressure_tau100_row(j) = allard_log10_pressure_tau100(i1,j1)  ! Get this row's four working P100L's
            temp_tau100_row(j) = allard_log10_temp_tau100(i1,j1)  ! Get this row's four working T100L's
         enddo
         j1 = gl_index(i)
         call polint(allard_gl_grid(j1),pressure_row,4,log_g,pressure_col(i),unused_dy)
         call polint(allard_gl_grid(j1),pressure_tau100_row,4,log_g,pressure_tau100_col(i),unused_dy)
         call polint(allard_gl_grid(j1),temp_tau100_row,4,log_g,temp_tau100_col(i),unused_dy)
      enddo
!     Now do final 4-point Lagrange inerpolations in TEFFL
      if (.not. allard_use_tau100) then
!        Current standard alternative. PL for TEFFL,GL, TL=TEFFL
         call polint(allard_teffl_grid(teffl_index),pressure_col,4,log_teff,atm_log10_pressure,unused_dy)
         atm_log10_temperature = log_teff
      else
         call polint(allard_teffl_grid(teffl_index),pressure_tau100_col,4,log_teff,atm_log10_pressure,unused_dy)
         call polint(allard_teffl_grid(teffl_index),temp_tau100_col,4,log_teff,atm_log10_temperature,unused_dy)
      endif

!     We now have obtained the needed temparatures and pressures

!     If requested, WRITE OUT INFORMATION TO THE MODEL FILE.
      if (print_to_files) then
         if (allard_use_tau100) then
            write(short_file_unit,70)
            write(istor,70)
   70       format('********PL,TL at Tau=100 INTERPOLATED' &
                 , ' FROM ALLARD TABULATED VALUES********')
         else
            write(short_file_unit,71)
            write(istor,71)
   71       format('********PL,TL at T=TEFF,GL INTERPOLATED ' &
                , 'FROM ALLARD TABULATED VALUES********')
         endif
         write(short_file_unit,72) log_teff,log_g,atm_log10_temperature,atm_log10_pressure
         write(istor,72) log_teff,log_g,atm_log10_temperature,atm_log10_pressure
   72    format(' ',20x,'LOG(Teffl) =',f10.5,' ,LOG(G) =' ,f10.5, &
              ', LOG(T) =',f10.5,', Log(P) =', f10.5)
      endif

!     If not requested via LPRT, WRITE OUT INFORMATION TO THE SHORT FILE.
      if (.not. print_to_files) then
        if (allard_use_tau100) then
           write(short_file_unit,73)
   73      format('ALSURFP:  PL,TL at Tau=100 INTERPOLATED' &
                , ' FROM ALLARD TABULATED VALUES:')
        else
           write(short_file_unit,74)
   74      format('ALSURFP: PL,TL at T=TEFF,GL INTERPOLATED ' &
                , 'FROM ALLARD TABULATED VALUES:')
        endif
        write(short_file_unit,75) log_teff,log_g,atm_log10_temperature,atm_log10_pressure
   75   format(' ',10x,'LOG(Teffl) =',f10.5,', LOG(G) =',f10.5, &
             ', LOG(T) =',f10.5,', Log(P) =',f10.5)
      endif




      return

 9999      continue
       write(*,*)
       write(*,*)'******** ALSURFP: Program Terminated ********'
       write(*,*)
       write(short_file_unit,*)
       write(short_file_unit,*)'******** ALSURF: Program Terminated ********'
       write(short_file_unit,*)
       stop     ! eliminate for testing

end subroutine alsurfp

! End if inputs are valid
