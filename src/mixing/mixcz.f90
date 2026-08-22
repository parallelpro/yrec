!----------------------------------------------------------------------
! mixcz
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mixcz.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Locates the boundaries of standard (non-overshoot) convection zones
! from convective_flag and homogenizes the composition within each
! one. Also renormalizes the composition to guard against small
! negative/overflowing abundances at the end.
!
! G Somers 6/14 originally added an IFSTCL ("first call") argument
! plus HR/HP/HD/HG/HS1 dummy arguments to support an active taucz
! (convective overturn timescale) calculation in this routine; G
! Somers 3/17 moved that calculation elsewhere (it is now passed in
! via common/ovrtrn/) and the taucz code below was commented out. KC
! 2025-05-31 correspondingly dropped those now-unused dummy arguments
! from the subroutine signature. The dead code is preserved as
! comments below, unmodified, for historical reference; none of it
! executes.
subroutine mixcz(composition, shell_mass, convective_flag, num_zones)

      use const_lib
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: shell_mass(json)
      logical, intent(inout) :: convective_flag(json)
      integer, intent(in) :: num_zones



! JVS 02/12 common blocks added for the calculation of taucz (now
! unused here -- see header note above; declared only to preserve
! layout).








      double precision :: species_sum(15)
      integer :: zone_bounds(24)
      integer :: j_idx, num_zones_plus1, zone_idx, zone_start, zone_end, &
           species_idx, inner_idx, num_species
      double precision :: weight_sum
      logical :: in_convection_zone

! COMPUTE STANDARD CONVECTION ZONES
      j_idx = 1
      in_convection_zone = .false.
      num_zones_plus1 = num_zones + 1
      convective_flag(num_zones_plus1) = .false.
      do zone_idx = 1, num_zones_plus1
         if (.not. (.not.convective_flag(zone_idx))) then
! CONVECTION ZONE
         if (in_convection_zone) cycle
! START OF CONVECTION ZONE
         in_convection_zone = .true.
         zone_start = zone_idx
         cycle
         end if
   10    if (.not.in_convection_zone) cycle
!   END OF CONVECTION ZONE
         in_convection_zone = .false.
         zone_bounds(j_idx) = zone_start
         zone_bounds(j_idx+1) = zone_idx - 1
         j_idx = j_idx + 2
         if (j_idx.lt.24) cycle
         exit
   11    continue
      end do
      if (zone_idx > num_zones_plus1) then
      zone_bounds(j_idx) = 0
      end if
   12 continue
      num_species = 11
      if (use_extended_composition) num_species = 15
! MIX ALL CONVECTIVE ZONES
      do j_idx = 1, 24, 2
         if (zone_bounds(j_idx).le.0) exit
         zone_start = zone_bounds(j_idx)
         zone_end = min0(num_zones, zone_bounds(j_idx+1))
         if (zone_start.ne.1 .and. zone_start.ge.zone_end) cycle
! INITIALIZE SUMS
         weight_sum = 0.0d0
         do species_idx = 1, num_species
            species_sum(species_idx) = 0.0d0
   40    continue
         end do
         do inner_idx = zone_start, zone_end
            weight_sum = weight_sum + shell_mass(inner_idx)
            do species_idx = 1, num_species
               species_sum(species_idx) = species_sum(species_idx) + &
                    composition(species_idx,inner_idx)*shell_mass(inner_idx)
   50       continue
            end do
   60    continue
         end do
         do species_idx = 1, num_species
            species_sum(species_idx) = species_sum(species_idx)/weight_sum
   70    continue
         end do
         do inner_idx = zone_start, zone_end
            do species_idx = 1, num_species
               composition(species_idx,inner_idx) = species_sum(species_idx)
   80       continue
            end do
   90    continue
         end do
  100 continue
      end do
  110 continue
! RENORMALIZE COMPOSITION IF NECESSARY
      do zone_idx = 1, num_zones
         composition(1,zone_idx) = dmax1(composition(1,zone_idx), 0.0d0)
         composition(3,zone_idx) = dmin1(composition(3,zone_idx), &
              1.0d0-composition(1,zone_idx))
         composition(9,zone_idx) = dmax1(composition(9,zone_idx), &
              0.99d-3*(composition(3,zone_idx)-star%env_comp%zenvm))
      end do

! G Somers 3/17, commented out this taucz calculation. It is now
! passed in in the OVRTRN common block.
!
! JVS 02/12 calculate the local convective overturn timescale at the
! base of the CZ. In older versions this was only done for rotating
! models; this makes it so taucz is calculated for all models.
! This code snagged from midmod
!
!  determine extent of surface convection zone.
!      lallcz = .false.
!      if(lcz(m))then
!  surface c.z. exists.  find lowest shell (imax), which is also the
!  uppermost zone considered for stability against rotationally induced mixing.
!         do 71 i = m-1,1,-1
!            if(.not.lcz(i)) goto 81
!   71    continue
!         lallcz = .true.
!         i = 0
!   81    imax = i + 1
!  hstop is the mass at the top of the c.z.
!  hsbot is the mass at the bottom of the c.z.
!         hstop = exp(cln*star%env_comp%stotal)
!         if(imax.gt.1)then
!            hsbot = 0.5d0*(hs1(imax)+hs1(imax-1))
!         else
!            hsbot = 0.0d0
!         endif
!  lczsur=t if a surface c.z.deep enough for angular momentum loss exists
!         if((hstop-hsbot)/cmsun.gt.0.0d0)then
!            lczsur = .true.
!         else
!            lczsur= .false.
!         endif
!      else
!  no surface c.z.
!         imax = m
!         lczsur = .false.
!      endif
!
!      if(lczsur)then
!         if(.not.lallcz)then
!            if(.not.lrot)then
!               hg(imax)=hs1(imax)*exp(cln*(cgl-2.0d0*hr(imax)))
!               hg(imax-1)=hs1(imax-1)*exp(cln*(cgl-2.0d0*hr(imax-1)))
!            endif
!  pinpoint rcz
!  g somers 6/14, check whether this run of mixcz occured before
!  or after midmod. if before, use sdel. if after, use the updated
!  variables delrm and delam.
!            if(ifstcl.eq..true.)then
!               dd2 = sdel(1,imax-1)-sdel(3,imax-1)
!               dd1 = sdel(1,imax)-sdel(3,imax)
!            else
!               dd2 = delrm(imax-1)-delam(imax-1)
!               dd1 = delrm(imax)-delam(imax)
!            endif
!  g somers end
!            fx = dd2/(dd2-dd1)
!  infer hp
!            envrl = hr(imax-1)+fx*(hr(imax)-hr(imax-1))
!            envr = exp(cln*envrl)
!            ps2 = exp(cln*(hp(imax)-hd(imax)))/hg(imax)
!            ps1 = exp(cln*(hp(imax-1)-hd(imax-1)))/hg(imax-1)
!            psca = ps1 + fx*(ps2-ps1)
!            rtestl = dlog10(envr+psca)
!  find v
!            do k = imax+1,m
!               if(hr(k).gt.rtestl)then
!                  fx = (rtestl-hr(k-1))/(hr(k)-hr(k-1))
!                  cvel = star%diag%svel(k-1)+fx*(star%diag%svel(k)-star%diag%svel(k-1))
!                  goto 85
!               endif
!            end do
!            cvel = star%diag%svel(m)
! 85         continue
!  define taucz
!            taucz = psca/cvel
!         else
!  infer hp
!            if(.not.lrot)then
!               hg(1)=hs1(1)*exp(cln*(cgl-2.0d0*hr(1)))
!            endif
!            psca2 = exp(cln*(hp(1)-hd(1)))/hg(1)
!            rtest2 = exp(cln*hr(1))
!            if(psca2.le.rtest2)then
!  hp < r at the first point.  assume v constant inside and hp = k/r for
!  slowly varying density and pressure near the center.
!               cvel = star%diag%svel(1)
!               psca = (psca2*rtest2)**0.5d0
!               taucz = psca/cvel
!            else
!               do k = 2,m
!                  psca1 = psca2
!                  rtest1 = rtest2
!                  if(.not.lrot)then
!                     hg(k)=hs1(k)*exp(cln*(cgl-2.0d0*hr(k)))
!                  endif
!                  psca2 = exp(cln*(hp(k)-hd(k)))/hg(k)
!                  rtest2 = exp(cln*hr(k))
!  find location where hp = r
!                  if(psca2.le.rtest2)then
!                     fx = (rtest1-psca1)/((psca2-rtest2)-(psca1-rtest1))
!  find v
!                     cvel = star%diag%svel(k-1)+fx*(star%diag%svel(k)-star%diag%svel(k-1))
!                     psca = psca1+fx*(psca2-psca1)
!  define taucz
!                     taucz = psca/cvel
!                     goto 95
!                  endif
!               end do
!               k = m
!               cvel = star%diag%svel(m)
!               psca = psca2
!               taucz = psca/cvel
! 95            continue
!            endif
!         endif
!      else
!         taucz = 0.0d0
!      endif

! end jvs

      return
end subroutine mixcz
