!----------------------------------------------------------------------
! pindex
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original pindex.f; only variable names, source form, and comment
! style were updated.
!
! Builds id(), the array of shell indices to be printed out (e.g. if
! shells 1,3,5,... are to be printed, id(1)=1, id(2)=3, ...). Prints
! every zone within the H-burning shell (jxbeg..jxend) if lshell is
! set, plus every print_point_interval-th point elsewhere.
subroutine pindex(jxbeg, jxend, lshell, m, id, idm)

      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: jxbeg, jxend
      logical, intent(in) :: lshell
      integer, intent(in) :: m
      integer, intent(out) :: id(json)
      integer, intent(out) :: idm

! common/ccout/: not used in this file's logic; layout placeholder.
! Naming matches ccoeft.f90/wrtout.f90.
      logical :: lstore, lstatm, lstenv, lstmod, lstphys, lstrot, lscrib, &
           lstch, lphhd
      common/ccout/ lstore, lstatm, lstenv, lstmod, lstphys, lstrot, &
           lscrib, lstch, lphhd

! common/ccout1/: only print_point_interval (originally NPRTPT) is
! used here. Naming matches wrtmil.f90/hpoint.f90.
      integer :: npenv, nprtmod, print_point_interval, npoint
      common/ccout1/ npenv, nprtmod, print_point_interval, npoint

      save

      integer :: ibeg, iend, ixbeg, j

!  ID IS THE ARRAY OF THE SHELLS TO BE PRINTED OUT;
!  I.E. IF SHELLS 1,3,5,..TO BE PRINTED ID(1)=1,ID(2)=3,..
      id(1) = 1
      idm = 2
      if (lshell) then
!  PRINT OUT EVERY ZONE IN H-BURNING SHELL
       if (print_point_interval.lt.jxbeg) then
          ibeg = max(2,print_point_interval)
          iend = int(jxbeg/print_point_interval)*print_point_interval
          do 10 j = ibeg,iend,print_point_interval
             id(idm) = j
             idm = idm + 1
   10       continue
       end if
       if (iend .eq. jxbeg) then
           ixbeg = jxbeg + 1
       else
           ixbeg = jxbeg
       end if
       do 20 j = ixbeg,jxend
          id(idm) = j
          idm = idm + 1
   20    continue
       if (print_point_interval.lt.m) then
          ibeg = int(jxend/print_point_interval+1)*print_point_interval
          iend = int(m/print_point_interval)*print_point_interval
          do 30 j = ibeg,iend,print_point_interval
             id(idm) = j
             idm = idm + 1
   30       continue
       end if
      else if (print_point_interval.lt.m) then
!  GENERAL CASE; PRINT OUT EVERY NPRTPT POINTS.
       ibeg = max(2,print_point_interval)
       iend = int(m/print_point_interval)*print_point_interval
       if (iend .eq. m) iend = iend - print_point_interval
       do 40 j = ibeg,iend,print_point_interval
          id(idm) = j
          idm = idm + 1
   40    continue
      end if
      id(idm) = m
      return
end subroutine pindex
