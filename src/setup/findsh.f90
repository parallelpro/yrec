!----------------------------------------------------------------------
! findsh
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original findsh.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! SR findsh locates the outer edge of a central convection zone and
! the inner edge of a surface c.z. It also locates the beginning,
! center, and end of a hydrogen-burning shell if applicable.
! ** to be added **
! locating the exact edges of central and surface c.z.'s and adding
! points at the edges.
!
! input variables:
! declared:
!   composition : run of mass fraction of the various chemical species.
!   luminosity : run of luminosity as a function of mass (solar units).
!   is_convective : flag set t if a shell is convective.
!   num_points : number of points in the model.
! in common blocks :
!   atime(1) : user parameter; if x < atime(1) then the program considers
!        the star to have a hydrogen-burning shell rather than a hydrogen-
!        burning core.
! local variables :
!   luminosity_change_tol : the outer edge of the h-burning shell is defined
!        as the point where the change in luminosity from one shell to the
!        next is less than endshl.
!   hydrogen_surface_tol : the outer edge of the h-burning shell is reached
!        when x differs from the surface value by less than
!        hydrogen_surface_tol regardless of the luminosity test.
!
! output variables :
!   core_edge : outermost point in the convective core (1=no convective core)
!   envelope_edge : innermost point in the surface c.z. (m=no surface c.z.)
!   has_h_shell : flag set t if the program considers the model to have a
!        hydrogen-burning shell.
!   if has_h_shell=t then the following are computed :
!   shell_begin : first shell with x > atime(1).
!   shell_mid : first shell where x exceeds 1/2 the surface value.
!   shell_end : last shell where l(i)-l(i-1) > endshl or x is within 1.0e-5
!        of the surface value.
subroutine findsh(composition, luminosity, is_convective, num_points, &
     core_edge, envelope_edge, shell_begin, shell_end, shell_mid, &
     has_h_shell)
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: composition(15,json), luminosity(json)
      logical, intent(in) :: is_convective(json)
      integer, intent(in) :: num_points
      integer, intent(out) :: core_edge, envelope_edge, shell_begin, &
           shell_end, shell_mid
      logical, intent(out) :: has_h_shell


      double precision :: luminosity_change_tol, hydrogen_surface_tol
      double precision :: half_surface_x, luminosity_end_threshold
      integer :: i
      data luminosity_change_tol,hydrogen_surface_tol/1.0d-5,1.0d-5/
      save

!ccc h-shell values
      shell_begin = 1
      shell_mid = 1
      shell_end = 1
      has_h_shell = .false.
!  if central x below threshold then calculate h shell values
      if(composition(1,1).le.atime(1)) then
       has_h_shell = .true.
       half_surface_x = 0.50d0*composition(1,num_points)
       luminosity_end_threshold = luminosity_change_tol*luminosity(num_points)
!  find beginning(shell_begin), middle(shell_mid) and end(shell_end) of h shell
       do i = 1,num_points
!          IF(HCOMP(1,I).LE.1.0D-10) THEN  ! Changed after discussion with Marc
          if(composition(1,i).le.atime(1)) then ! to force consistency with above LLP 9/24/08
             shell_begin = shell_begin+1
             shell_mid = shell_mid+1
          else if(composition(1,i).le.half_surface_x) then
             shell_mid = shell_mid+1
          else if(luminosity(i) - luminosity(i-1).lt.luminosity_end_threshold) then
!               write(*,*)'luminosity criteria'
             goto 20
          else if(composition(1,num_points) - composition(1,i).lt.hydrogen_surface_tol) then
!               write(*,*)'composition criteria'
             goto 20
          endif
   10    continue
       end do
       i = num_points
   20    shell_end = i
      endif
!ccc find boundary of central convection zone.
      do i = 1,num_points
       if(.not.is_convective(i)) goto 40
   30 continue
      end do
   40 if(i.gt.1) then
       core_edge = i-1
      else
       core_edge = 1
      endif
!ccc find boundary of surface c.z.
      do i = num_points,1,-1
       if(.not.is_convective(i)) goto 60
   50 continue
      end do
   60 if(i.lt.num_points) then
       envelope_edge = i+1
      else
       envelope_edge = num_points
      endif
!  for a fully convective star (envelope_edge=1),turn the convective core off(core_edge=1).
!  this is done because core_edge is used in computing the nuclear timestep,
!  and putting it at the surface does strange things.
      if(envelope_edge.eq.1) core_edge = 1
      return
end subroutine findsh
