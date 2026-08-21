!----------------------------------------------------------------------
!  YREC version value assignment
!----------------------------------------------------------------------
#define YREC_VERSION '202602'

#ifndef GIT_HASH
#define GIT_HASH 'unknown'
#endif

subroutine setversion()

      use const_lib
      implicit none
! former common/version/: yrec_version_string (released version tag)
! and git_hash_string (short git commit hash + indicator if working
! tree was not clean) are now use-associated from const_lib.
      yrec_version_string = YREC_VERSION
      git_hash_string = GIT_HASH
end subroutine setversion
