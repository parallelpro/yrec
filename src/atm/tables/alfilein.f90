!----------------------------------------------------------------------
! alfilein
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alfilein.f; only variable names, source form, and comment
! style were updated. FORMAT strings and the READ/OPEN statements
! that must byte-for-byte match the on-disk table file are preserved
! verbatim.
!
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
! Allard Atmosphere table lookup routine, alsurfp.f90. In addition, we fill any
! vacant entries in the table with -999.0. A value of -999. in the table
! identifies the entry as "invalid." In addition, key auxililiary arrays, as
! described below, are created.
subroutine alfilein(allard_table_path, ierr)
!
! Parameters NTA and NGA are respectively the maximum expected numbers
! of Teff's and GL's we expect to encounter, even in future tables
! The maximum nuber of Teff's in the current table is nTeff, and the
! associated max number of GL's is nGL.

      use atm_table_lib
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: nta = 250
      integer, parameter :: nga = 25

      character(len=256), intent(in) :: allard_table_path
      character(len=35) :: first_record
      character(len=256) :: header_line
      double precision :: local_teffs(nta)



      external sort_shell
      integer :: i, j, nhdr, irecno, i1, j1
      double precision :: teff_value, gl_value, feh_value, alpha_value
      double precision :: pressure_value, pressure_tau100_value, temp_tau100_value
! NOTE: the original alfilein.f has a stray blank line followed by a
! badly-indented "     EXTERNAL sort_shell" whose column-6 character
! lands in the fixed-form continuation column, silently merging it
! into the end of the preceding COMMON /ALATM04/ statement (turning
! "...,DUMMY4" into a single garbled identifier
! "DUMMY4XTERNALSORT_SHELL" -- still one double-precision scalar, so
! the /ALATM04/ storage layout is unaffected). The net effect is that
! the original file never actually declares sort_shell EXTERNAL; that
! has no behavioral consequence since sort_shell is only CALLed
! directly, never passed as an actual argument. This free-form
! version declares it properly above (a no-op semantically) rather
! than reproducing the column-6 accident, since there is no
! free-form equivalent of that fixed-form parsing quirk to preserve.
!
! NOTE: line ~105 below tests a variable spelled LATMTPTau100, which
! is NOT the same as the common-block flag LALTPTau100 (allard_use_tau100,
! common/alatm03/) despite the similar name -- this looks like a typo
! in the original code that silently declares a second, distinct,
! never-assigned local logical (implicitly typed, SAVE'd) instead of
! referencing the intended common flag. Preserved exactly as a
! separate, always-default-valued local below.
      logical :: latmtptau100

!     set output arrays to invalid values
      integer, intent(out) :: ierr

      ierr = 0

      do i = 1,nta
         atm_table%allard_teffl_grid(i) = -999.d0
         do j = 1,nga
            atm_table%allard_gl_grid(j) = -999.d0
            atm_table%allard_log10_pressure(i,j) = -999.d0
            atm_table%allard_log10_pressure_tau100(i,j) = -999.d0
            atm_table%allard_log10_temp_tau100(i,j) = -999.d0
         enddo
      enddo
!     the input file name is in FALLARD and its unit number is IOATMA
      open(allard_table_unit,file=allard_table_path,status='OLD',err= 899)
      goto 900

  899      continue
       write(*,*)
       write(*,*) 'ALFILEIN: Allard File OPEN Failure'
       write(*,*) '     Unit-Number, FIle-Name are: ', &
           allard_table_unit,allard_table_path
       goto 9999      ! Error exit

  900      continue
!     Decide what kind of input file
      read(allard_table_unit,901) first_record
  901  format(a)
      rewind(allard_table_unit)
      if (first_record(4:7) .eq. '3.50') goto 100  ! Old NextGen file type
      if (first_record(4:9) .eq. 'ALLARD') goto 200  ! Allard file type
!     If we get here, the input file is invalid
      write(short_file_unit,*)'*** Invalid Allard Atmosphere input file ***'
      write(short_file_unit,*) "First record was '",first_record,"'"
      write(*,*)'*** Invalid Allard Atmosphere input file ***'
      write(*,*) "First record was '",first_record,"'"
      goto 9999   ! The error exit

!     Process old-stype nextgen input file.
  100      continue
       atm_table%allard_is_old_nextgen = .true.
       write(short_file_unit,*) 'ALFileIn: File Description: 1999 NEXTGEN', &
           ' (Old NEXTGEN)'

!     Ensure that we are not requesting PL,TL at Tau-100. These are not
!     present in the old Nextgen 1 files. If requested, fail
      if (latmtptau100) then
         write(*,*)
         write(*,*) 'ALFileIN: Invalid old Allard input file ', &
               'for requested PT,TL at Tau=100'
         write(short_file_unit,*)
         write(short_file_unit,*) 'ALFileIN: Invalid old Allard input file ', &
               'for requested PT,TL at Tau=100'
         goto 9999  ! The ERROR EXIT
      endif
!     READ RANGE OF GRAVITIES
      atm_table%allard_num_teff = 54
      atm_table%allard_num_gl = 5
      read(allard_table_unit,905) (atm_table%allard_gl_grid(i),i=1,atm_table%allard_num_gl)
  905  format(5f7.2/)
!       read(ioatma,*)

!     READ RANGE OF TEMPERATURES
      read(allard_table_unit,907) (atm_table%allard_teffl_grid(i),i=1,atm_table%allard_num_teff)
  907  format(10(1p5e16.8,/),1p4e16.8,/)
!       read(ioatma,*)
      read(allard_table_unit,909) ((atm_table%allard_log10_pressure(j,i),i=1,atm_table%allard_num_gl),j=1,atm_table%allard_num_teff)
  909  format(1p5e16.8)
      close(allard_table_unit)
      goto 500      ! NORMAL EXIT
!     Old Nextgen input file has been read an put into the proper arrays

!     Process new-style Allard Atmosphere file
  200      continue
       atm_table%allard_is_old_nextgen = .false.
      read(allard_table_unit,911) nhdr, header_line  ! nhdr is number of header recoreds to skip over
  911      format(i2,x,a)
       write(short_file_unit,912) 'ALFilein: New Allard Atm: File Description: ',header_line(1:47)
  912  format(2a)
!     Skip over rest of header, if any
      if (nhdr .gt. 1) then
        do i = 2, nhdr
          read(allard_table_unit,*)
        enddo
      endif
      atm_table%allard_num_gl = 0
      atm_table%allard_num_teff = 0
      atm_table%allard_num_feh = 0
      atm_table%allard_num_alpha = 0
      irecno = 0

!     First we read the file to count the # of Teff,GL,Fe/H and Alpha's

  299      read(allard_table_unit,915,end=400,err=399) teff_value,gl_value,feh_value,alpha_value ! Process first record(s)
  915      format(f6.0,3f6.2)
      irecno=irecno+1

!     If the record does not have the correct FeH and Alpha, we skip it
      if ((dabs(feh_value-allard_target_feh).gt.1d-5) .or. &
             (dabs(alpha_value-allard_target_alpha).gt.0d0)) then
         goto 299  ! On to next record. Stay in this part until we get an acceptable record
      endif
      atm_table%allard_num_teff = 1    ! Initialize counters and variables after reading first record
      atm_table%allard_num_gl = 1
      atm_table%allard_num_feh = 1
      atm_table%allard_num_alpha = 1
      local_teffs(atm_table%allard_num_teff) = teff_value
      atm_table%allard_gl_grid(atm_table%allard_num_gl) = gl_value
      atm_table%allard_feh_grid(atm_table%allard_num_feh) = feh_value
      atm_table%allard_alpha_grid(atm_table%allard_num_alpha) = alpha_value

  300      read(allard_table_unit,915,end=400,err=399) teff_value,gl_value,feh_value,alpha_value
      irecno=irecno+1

!     If the record does not have the correct FeH and Alpha, we skip it
      if ((dabs(feh_value-allard_target_feh).gt.1d-5) .or. &
             (dabs(alpha_value-allard_target_alpha).gt.0d0)) then
         goto 316  ! On to next record
      endif

!      Now check Teffs, increase counter if needed
       do i = 1,atm_table%allard_num_teff  ! Skip out if any old Teff is a match
        if (dabs(teff_value-local_teffs(i)) .lt.1d-6) goto 310
      enddo
        atm_table%allard_num_teff = atm_table%allard_num_teff+1      ! count the nextTeff's
        local_teffs(atm_table%allard_num_teff) = teff_value  ! save the new Teff in array Teffs
  310  continue

!     Now check GLs, increase counter if needed
      do i = 1,atm_table%allard_num_gl  ! Skip out if any old GL is a match
        if (dabs(gl_value-atm_table%allard_gl_grid(i)) .lt.1d-6) goto 312
      enddo
        atm_table%allard_num_gl = atm_table%allard_num_gl+1      ! count the next GL
        atm_table%allard_gl_grid(atm_table%allard_num_gl) = gl_value  ! save the new GL in array GLs
  312  continue


!     Now check FeHs, increase counter if needed
      do i = 1,atm_table%allard_num_gl  ! Skip out if any old FeH is a match
        if (dabs(feh_value-atm_table%allard_feh_grid(i)) .lt.1d-6) goto 314
      enddo
        atm_table%allard_num_feh = atm_table%allard_num_feh+1      ! count the next FeH
        atm_table%allard_feh_grid(atm_table%allard_num_feh) = feh_value  ! save the new FeH in array FeHs
  314  continue

!     Now check ALPHAs, increase counter if needed
      do i = 1,atm_table%allard_num_gl  ! Skip out if any old ALPHA is a match
        if (dabs(alpha_value-atm_table%allard_alpha_grid(i)) .lt.1d-6) goto 316
      enddo
        atm_table%allard_num_alpha = atm_table%allard_num_alpha+1      ! count the next ALPHA
        atm_table%allard_alpha_grid(atm_table%allard_num_alpha) = alpha_value  ! save the new ALPHA in array ALPHAs
  316  continue

      goto 300 ! go on to next record

!     File Read error exit
  399      continue
       write(*,*)
       write(*,*)'AtTabInit Pass 1 File Error '
       write(*,*)
       goto 9999

  400      continue
!       write(*,*) 'nTeff,nGL,nFeH,nALPHA: ',nTeff,nGL,nFeH,nALPHA

!      Now we have unsorted Teff,GL,Fe/H and Alpha's.
!     The next step is to sort them
      call sort_shell(atm_table%allard_num_teff,local_teffs)
      call sort_shell(atm_table%allard_num_gl,atm_table%allard_gl_grid)
      call sort_shell(atm_table%allard_num_feh,atm_table%allard_feh_grid)
      call sort_shell(atm_table%allard_num_alpha,atm_table%allard_alpha_grid)

!     Now we can rewind the input file and start over for real
      rewind(allard_table_unit)  ! go back to the beginning of the file

!      Skip over header
      do i = 1,nhdr
         read(allard_table_unit,*)
      enddo

  410       read(allard_table_unit,920,end=500,err=499) teff_value,gl_value,feh_value,alpha_value, &
              pressure_value,pressure_tau100_value,temp_tau100_value
  920      format(f6.0,3f6.2,1p4d16.8)

!     If the record does not have the correct FeH and Alpha, we skip it
      if ((dabs(feh_value-allard_target_feh).gt.1d-5) .or. &
             (dabs(alpha_value-allard_target_alpha).gt.0d0)) then
         goto 440  ! On to next record
      endif

      do i = 1,atm_table%allard_num_teff  ! Skip out when we find the matching Teff
        if (dabs(teff_value-local_teffs(i)) .lt.1d-6) then
           i1 = i
           goto 420
        endif
      enddo
      write(*,*) 'ALFilein: Impossible failure #1'
      write(*,*) 'Teff of ', teff_value, ' not in Teff table.'
      goto 9999  ! The error exit

  420  continue
      do j = 1,atm_table%allard_num_gl  ! Skip out when we find the a matching GL
        if (dabs(gl_value-atm_table%allard_gl_grid(j)) .lt.1d-6) then
           j1 = j
           goto 430
        endif
      enddo
      write(*,*) 'ALFilein: Impossible failure #2'
      write(*,*) 'GL of ', gl_value, ' not in GL table.'
      goto 9999  ! The error exit

  430      continue

!     we now verify that we have the correct FeH and Alpha

      if (dabs(feh_value -allard_target_feh) .ge. 1d-6)  goto 440     ! Wrong FeH, bypass item
      if (dabs(alpha_value -allard_target_alpha) .ge. 1d-6)  goto 440 ! Wrong Alpha, bypass item


!     We now have the correct indices for our tables, i1 for the Teff-direction
!     and j1 for the GL direction

       atm_table%allard_teffl_grid(i1) = log10(teff_value)  ! We need log(teff), but up to now we've been using Teff
       atm_table%allard_gl_grid(j1) = gl_value
       atm_table%allard_log10_pressure(i1,j1) = pressure_value
       atm_table%allard_log10_pressure_tau100(i1,j1) = pressure_tau100_value
       atm_table%allard_log10_temp_tau100(i1,j1) = temp_tau100_value

  440  continue
       goto 410   ! Back to process the next record

  499      continue
       write(*,*)
       write(*,*)'AtTabInit File Read error'
       write(*,*)
       goto 9999

  500      continue  ! We heve finished with the input file and entered all inputs

      call altabinit(ierr)   ! Initialize Allard tables
      if (ierr /= 0) return

       return

 9999      continue    ! THE ERROR EXIT
       write(*,*)
       write(*,*)
       write(*,*) '**************** PROGRAM ALFilein TERMINATED ', &
             '@ 9999 **************'
       write(*,*)
       write(*,*)
       write(*,*)
       write(short_file_unit,*)
       write(short_file_unit,*) '**************** PROGRAM ALFilein TERMINATED ', &
            '@ 9999 *************'
       write(short_file_unit,*)
       write(short_file_unit,*)
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the atm_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return

end subroutine alfilein      ! END OF ALINITTAB

! ************* End of alfilein *************************


! Numerical Recipes Shell sort

subroutine sort_shell(num_elements, values)

      implicit none
      integer, intent(in) :: num_elements
      double precision, intent(inout) :: values(num_elements)
! Sorts an array A(1:nn) into ascending numerical order by Shell's method (diminishing
! incremental sort) nn is input. A is replaced on output by its sorted arrangement
      integer :: n, inc, i, j
      double precision :: v

      if (num_elements .eq. 1) return                  ! Exit if only one element
      n=num_elements
      inc=1
    1 inc=3*inc+1                         ! Determin starting increment
      if (inc .le. n) goto 1
    2      continue
         inc=inc/3                        ! Loop over the partial sorts
         do i=inc+1,n
            v=values(i)                        ! Outer loop of straight insertion
            j=i
    3        if (values(j-inc) .gt. v) then      ! Inner loop of straight insertion
               values(j)=values(j-inc)
               j=j-inc
               if (j .le. inc) goto 4
            goto 3
            endif
    4       values(j)=v
         enddo
      if (inc .gt. 1) goto 2
      return
end subroutine sort_shell

! ************* End of sort_shell *************************
