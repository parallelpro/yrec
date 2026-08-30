!----------------------------------------------------------------------
!  YREC version value assignment
!----------------------------------------------------------------------
#define YREC_VERSION '202602'

#ifndef GIT_HASH
#define GIT_HASH 'unknown'
#endif

subroutine setversion()

      use phys_const_lib
      implicit none
! former common/version/: yrec_version_string (released version tag)
! and git_hash_string (short git commit hash + indicator if working
! tree was not clean) are now use-associated from const_lib.
      yrec_version_string = YREC_VERSION
      git_hash_string = GIT_HASH
! Builds outside a git checkout can define GIT_HASH as an empty (or
! '-dirty'-only) string, which the ifndef above cannot catch; an
! empty hash used to produce a fatal zero-width A00 edit descriptor
! in the run-log banner.
      if (len_trim(git_hash_string) == 0) git_hash_string = 'unknown'
end subroutine setversion
