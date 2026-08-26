!----------------------------------------------------------------------
! wtime
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wtime.f; only variable names, source form, and comment
! style were updated.
!
! For rotating models, limit the timestep based on max change in
! omega. The user parameter max_domega_global (originally DTWIND)
! governs the maximum change allowed between models.
subroutine timestep_limit_omega(previous_timestep, num_points, omega, rotation_dt, &
     max_domega_frac)

      use star_info_lib, only: star, json
      implicit none

! previous_timestep: previous model timestep.
! num_points: number of model points.
! omega: run of angular velocity in the current model.
! rotation_dt: rotation-based timestep.
! max_domega_frac: maximum fractional change in omega.
      double precision, intent(in) :: previous_timestep
      integer, intent(in) :: num_points
      double precision, intent(in) :: omega(json)
      double precision, intent(out) :: rotation_dt
      double precision, intent(out) :: max_domega_frac
      integer :: start_index, i, max_index
      double precision :: test_domega, dt_factor, dt_factor_limit

      start_index = 1
      max_domega_frac = 2.0d0*abs(omega(start_index)-star%old_omega(start_index))/ &
           (omega(start_index)+star%old_omega(start_index))
      max_index = start_index
      do i = start_index+1,num_points
         test_domega=2.0d0*abs(omega(i)-star%old_omega(i))/(omega(i)+star%old_omega(i))
         if(test_domega.gt.max_domega_frac) then
            max_domega_frac = test_domega
            max_index = i
         endif
      end do
      dt_factor = max_domega_frac/star%job%max_domega_global
! if no change from previous model,set rotation_dt to timestep
! stored in the previous model.
      if (dt_factor.eq.0.d0)then
          rotation_dt = 1.0d20
          continue
          return
      endif
! restrict change in timestep to no more than a factor of atime(14)%
! up or down.
! mhp 10/14 use atime(13) as the global timestep limiter for the code
      dt_factor_limit = star%ctrl%atime(13)
      if (dt_factor.gt.dt_factor_limit) dt_factor=dt_factor_limit
      if (dt_factor.lt.1.0d0/dt_factor_limit) dt_factor=1.0d0/dt_factor_limit
      rotation_dt = previous_timestep/dt_factor
      return
end subroutine timestep_limit_omega
