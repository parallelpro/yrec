!----------------------------------------------------------------------
!  YREC version value assignment
!----------------------------------------------------------------------
#define YREC_VERSION '202602'

#ifndef GIT_HASH
#define GIT_HASH 'unknown'
#endif

subroutine setversion()

      implicit none
! yrec_version_string: released version tag.
! git_hash_string: short git commit hash + indicator if working tree
!                   was not clean.
      character*10 yrec_version_string
      character*20 git_hash_string
      common/version/ yrec_version_string, git_hash_string
      yrec_version_string = YREC_VERSION
      git_hash_string = GIT_HASH
end subroutine setversion
