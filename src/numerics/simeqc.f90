!----------------------------------------------------------------------
! simeqc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original simeqc.f; only variable names, source form, and comment
! style were updated.
!
! Gauss-Jordan elimination with partial pivoting, operating on a
! system matrix stored as a flat array: num_unknowns equations (rows),
! num_cols columns (num_cols > num_unknowns for one or more augmented
! right-hand-side columns), stored column-major in system_matrix(56).
subroutine simeqc(system_matrix, num_cols, num_unknowns, ierr)

      implicit none

      double precision, intent(inout) :: system_matrix(56)
      integer, intent(in) :: num_cols, num_unknowns
      integer :: jj, j, jy, it, i, ij, imax, ia, ib, ic, id, ix, jx, ny, &
           n1, ig, ih
      double precision :: biga, swap_val

      integer, intent(out) :: ierr

      ierr = 0

      jj=-num_unknowns
      do j=1,num_unknowns
      jy=j+1
      jj=jj+num_unknowns+1
      biga=0.0d0
      it=jj-j
      do i=j,num_unknowns
      ij=it+i
      if(dabs(biga).ge.dabs(system_matrix(ij))) cycle
      biga=system_matrix(ij)
      imax=i
   30 continue
      end do
      if(dabs(biga).eq.0.0d0) goto 1010
      goto 1012
 1010 write (5,1011)
 1011 format (1x,'STOPPED AT 1010')
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
      ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
      ! preserve the historical stop on a nonzero return.
      ierr = 1
      return
 1012 ia=j+num_unknowns*(j-2)
      it=imax-j
      do i=j,num_cols
      ia=ia+num_unknowns
      ib=ia+it
      swap_val=system_matrix(ia)
      system_matrix(ia)=system_matrix(ib)
      system_matrix(ib)=swap_val
      system_matrix(ia)=system_matrix(ia)/biga
   50 continue
      end do
      if(j.eq.num_unknowns) exit
      ia=num_unknowns*(j-1)
      do ix=jy,num_unknowns
      ib=ia+ix
      it=j-ix
      do jx=jy,num_cols
      ic=num_unknowns*(jx-1)+ix
      id=ic+it
      system_matrix(ic)=system_matrix(ic)-system_matrix(ib)*system_matrix(id)
   60 continue
      end do
   64 continue
      end do
   65 continue
      end do
   70 ny=num_unknowns-1
      it=num_unknowns*num_unknowns
      do j=1,ny
      ia=it-j
      ic=num_unknowns*num_cols
      ib=ic-j
      do i=1,j
      system_matrix(ib)=system_matrix(ib)-system_matrix(ia)*system_matrix(ic)
      n1=num_cols-1
      ig=ib
      ih=ic
   75 if(n1.le.num_unknowns) go to 78
      ig=ig-num_unknowns
      ih=ih-num_unknowns
      system_matrix(ig)=system_matrix(ig)-system_matrix(ia)*system_matrix(ih)
      n1=n1-1
      go to 75
   78 ia=ia-num_unknowns
      ic=ic-1
   80 continue
      end do
   85 continue
      end do
      return
end subroutine simeqc
