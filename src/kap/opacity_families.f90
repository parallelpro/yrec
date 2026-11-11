!----------------------------------------------------------------------
! opacity_families
!----------------------------------------------------------------------
! The opacity table families, as data: each role's families in
! precedence order, and the family of each role that a run uses. This is
! the one place that orders them; kap_init loads, kap_update_surface_tables
! refreshes and kap_eval looks up only the families select_families
! returns.
!   interior   the atomic tables, used above molecular_opacity_logt_min:
!              OPAL 1995-format, else OPAL 1992, else LAOL 1989
!   molecular  the low-temperature tables, used below
!              molecular_opacity_logt_max: Alexander 2006, else
!              Alexander 1994 (the alex95 controls), else Kurucz 1990
!   helium core  the pure-Z table, toward which a helium-burning core's
!              interior opacity is blended
!   conduction the Potekhin conductive opacity
! A family that is selected but not first in its role's list is neither
! loaded nor used.
module opacity_families
      use opacity_tables, only: kap_settings_t
      use opacity_source, only: kap_source
      use opal95_opacity, only: opal95_source
      use opal92_opacity, only: opal92_source
      use laol_opacity, only: laol89_source, pure_z_source
      use alex06_opacity, only: alex06_source
      use alex94_opacity, only: alex94_source
      use kurucz_opacity, only: kurucz90_source
      use conductive_opacity, only: conductive_source
      implicit none
      private
      public :: interior_families, molecular_families, select_families

contains

! the interior families, in precedence order
function interior_families(settings) result(families)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: families(3)

      families(1) = opal95_source(settings)
      families(2) = opal92_source(settings)
      families(3) = laol89_source(settings)
end function interior_families

! the molecular families, in precedence order
function molecular_families(settings) result(families)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source) :: families(3)

      families(1) = alex06_source(settings)
      families(2) = alex94_source(settings)
      families(3) = kurucz90_source(settings)
end function molecular_families

! The families a run with these settings uses: the first selected of each
! role's list (a descriptor named 'none', not selected, when the role has
! none), the pure-Z table and the conductive opacity, with the windows
! of the molecular ramp: the molecular family below
! molecular_opacity_logt_max, the interior one above
! molecular_opacity_logt_min.
subroutine select_families(settings, interior, molecular, helium_core, conduction)
      type(kap_settings_t), intent(in) :: settings
      type(kap_source), intent(out) :: interior, molecular, helium_core, conduction

      interior = first_selected(interior_families(settings))
      interior%logt_min = settings%molecular_opacity_logt_min
      molecular = first_selected(molecular_families(settings))
      molecular%logt_max = settings%molecular_opacity_logt_max
      helium_core = pure_z_source(settings)
      conduction = conductive_source(settings)
end subroutine select_families

function first_selected(families) result(first)
      type(kap_source), intent(in) :: families(:)
      type(kap_source) :: first
      integer :: k

      do k = 1, size(families)
         if (families(k)%selected) then
            first = families(k)
            return
         end if
      end do
end function first_selected

end module opacity_families
