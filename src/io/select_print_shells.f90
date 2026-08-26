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
subroutine select_print_shells(jxbeg, jxend, lshell, m, id, idm)
      use star_info_lib, only: star
      use star_info_lib, only: star, json
      implicit none

      integer, intent(in) :: jxbeg, jxend
      logical, intent(in) :: lshell
      integer, intent(in) :: m
      integer, intent(out) :: id(json)
      integer, intent(out) :: idm
      integer :: ibeg, iend, ixbeg, j

!  ID IS THE ARRAY OF THE SHELLS TO BE PRINTED OUT;
!  I.E. IF SHELLS 1,3,5,..TO BE PRINTED ID(1)=1,ID(2)=3,..
      id(1) = 1
      idm = 2
      if (lshell) then
!  PRINT OUT EVERY ZONE IN H-BURNING SHELL
       if (star%ctrl%print_point_interval.lt.jxbeg) then
          ibeg = max(2,star%ctrl%print_point_interval)
          iend = int(jxbeg/star%ctrl%print_point_interval)*star%ctrl%print_point_interval
          do j = ibeg,iend,star%ctrl%print_point_interval
             id(idm) = j
             idm = idm + 1
          end do
       end if
       if (iend .eq. jxbeg) then
           ixbeg = jxbeg + 1
       else
           ixbeg = jxbeg
       end if
       do j = ixbeg,jxend
          id(idm) = j
          idm = idm + 1
       end do
       if (star%ctrl%print_point_interval.lt.m) then
          ibeg = int(jxend/star%ctrl%print_point_interval+1)*star%ctrl%print_point_interval
          iend = int(m/star%ctrl%print_point_interval)*star%ctrl%print_point_interval
          do j = ibeg,iend,star%ctrl%print_point_interval
             id(idm) = j
             idm = idm + 1
          end do
       end if
      else if (star%ctrl%print_point_interval.lt.m) then
!  GENERAL CASE; PRINT OUT EVERY NPRTPT POINTS.
       ibeg = max(2,star%ctrl%print_point_interval)
       iend = int(m/star%ctrl%print_point_interval)*star%ctrl%print_point_interval
       if (iend .eq. m) iend = iend - star%ctrl%print_point_interval
       do j = ibeg,iend,star%ctrl%print_point_interval
          id(idm) = j
          idm = idm + 1
       end do
      end if
      id(idm) = m
      return
end subroutine select_print_shells
