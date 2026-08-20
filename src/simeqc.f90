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
subroutine simeqc(system_matrix, num_cols, num_unknowns)

      implicit none

      double precision, intent(inout) :: system_matrix(56)
      integer, intent(in) :: num_cols, num_unknowns

      save

      integer :: jj, j, jy, it, i, ij, imax, ia, ib, ic, id, ix, jx, ny, &
           n1, ig, ih
      double precision :: biga, swap_val

      jj=-num_unknowns
      do 65 j=1,num_unknowns
      jy=j+1
      jj=jj+num_unknowns+1
      biga=0.0d0
      it=jj-j
      do 30 i=j,num_unknowns
      ij=it+i
      if(dabs(biga).ge.dabs(system_matrix(ij))) go to 30
      biga=system_matrix(ij)
      imax=i
   30 continue
      if(dabs(biga).eq.0.0d0) goto 1010
      goto 1012
 1010 write (5,1011)
 1011 format (1x,'STOPPED AT 1010')
      stop 29
 1012 ia=j+num_unknowns*(j-2)
      it=imax-j
      do 50 i=j,num_cols
      ia=ia+num_unknowns
      ib=ia+it
      swap_val=system_matrix(ia)
      system_matrix(ia)=system_matrix(ib)
      system_matrix(ib)=swap_val
      system_matrix(ia)=system_matrix(ia)/biga
   50 continue
      if(j.eq.num_unknowns) go to 70
      ia=num_unknowns*(j-1)
      do 64 ix=jy,num_unknowns
      ib=ia+ix
      it=j-ix
      do 60 jx=jy,num_cols
      ic=num_unknowns*(jx-1)+ix
      id=ic+it
      system_matrix(ic)=system_matrix(ic)-system_matrix(ib)*system_matrix(id)
   60 continue
   64 continue
   65 continue
   70 ny=num_unknowns-1
      it=num_unknowns*num_unknowns
      do 85 j=1,ny
      ia=it-j
      ic=num_unknowns*num_cols
      ib=ic-j
      do 80 i=1,j
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
   85 continue
      return
end subroutine simeqc
