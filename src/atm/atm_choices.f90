!----------------------------------------------------------------------
! atm_choices
!----------------------------------------------------------------------
! What each atm_choice is, in one place: its name (the run log's
! "surface boundary:" line), whether it is a table, which
! atm_get_surface_pt looks up and atm_init reads, or a T-tau relation,
! which the star layer's envelope integrates (t_tau_relations), and
! where a table hands the surface over to the gray relation. atm_lib
! re-exports the choices' names and the queries; the star layer and io
! ask them instead of testing choice numbers.
!     atm_choice_is_known      a choice atm_init accepts (0 to atm_last_choice)
!     atm_choice_is_tabulated  a table (3, 4, 5)
!     atm_choice_name          its name in the run log
!     has_gray_switch, kurucz_gray_switch_log_teff
!                              whether the star layer replaces the table
!                              by the gray relation above a log Teff, and
!                              the Kurucz table's (atm_lib's
!                              atm_leaves_table and atm_reenters_table add
!                              the Allard tables', which their reader
!                              finds)
module atm_choices
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: atm_choice_is_known, atm_choice_is_tabulated, atm_choice_name, &
           has_gray_switch

! The choices (the control atm_choice).
      integer, parameter, public :: atm_eddington_gray = 0, atm_krishna_swamy = 1, &
           atm_hsra = 2, atm_kurucz = 3, atm_allard = 4, atm_castelli_kurucz = 5
! the first and the last choice (atm_last_choice: callers loop over the
! choices with it)
      integer, parameter :: first_choice = atm_eddington_gray
      integer, parameter, public :: atm_last_choice = atm_castelli_kurucz

! The log Teff at and above which the star layer replaces the Kurucz
! table by the gray relation, and below which it returns to it (MHP
! 9/01; the table itself reaches log Teff 4.6). The Allard
! tables switch at their last log Teff; the Castelli/Kurucz table and the
! T-tau relations never.
      real(dp), parameter, public :: kurucz_gray_switch_log_teff = 3.95d0

      integer, parameter, public :: atm_choice_name_length = 44
      type :: choice_description
         character(len=atm_choice_name_length) :: name
         logical :: tabulated
         logical :: gray_switch
      end type choice_description
      type(choice_description), parameter :: choices(first_choice:atm_last_choice) = [ &
           choice_description('Eddington gray T(tau) relation', .false., .false.), &
           choice_description('Krishna-Swamy T(tau) relation', .false., .false.), &
           choice_description('Harvard-Smithsonian reference atmosphere', .false., .false.), &
           choice_description('Kurucz atmosphere tables', .true., .true.), &
           choice_description('Allard atmosphere tables', .true., .true.), &
           choice_description('Kurucz/Castelli atmosphere tables', .true., .false.)]

contains

! atm_init accepts the choice
pure logical function atm_choice_is_known(atm_choice)
      integer, intent(in) :: atm_choice
      atm_choice_is_known = atm_choice >= first_choice .and. atm_choice <= atm_last_choice
end function atm_choice_is_known

! the choice is a table (atm_get_surface_pt), not a T-tau relation
pure logical function atm_choice_is_tabulated(atm_choice)
      integer, intent(in) :: atm_choice
      atm_choice_is_tabulated = .false.
      if (atm_choice_is_known(atm_choice)) atm_choice_is_tabulated = choices(atm_choice)%tabulated
end function atm_choice_is_tabulated

! the star layer replaces the table by the gray relation above a log Teff
pure logical function has_gray_switch(atm_choice)
      integer, intent(in) :: atm_choice
      has_gray_switch = .false.
      if (atm_choice_is_known(atm_choice)) has_gray_switch = choices(atm_choice)%gray_switch
end function has_gray_switch

! the choice's name in the run log's banner ('' for an unknown choice)
pure function atm_choice_name(atm_choice) result(name)
      integer, intent(in) :: atm_choice
      character(len=atm_choice_name_length) :: name
      name = ''
      if (atm_choice_is_known(atm_choice)) name = choices(atm_choice)%name
end function atm_choice_name

!======================================================================
! Sources: the F77's KTTAU codes (0-5; 5, Castelli/Kurucz, JNT 6/14),
! PARMIN's run-log banner, and SURFBC's switch (MHP 9/01).

end module atm_choices
