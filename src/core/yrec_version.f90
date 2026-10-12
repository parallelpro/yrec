!----------------------------------------------------------------------
! yrec_version
!----------------------------------------------------------------------
! The code's version stamp: set_version_strings fills physical_constants'
! yrec_version_string (the release tag) and git_hash_string (the build's
! commit, from the Makefile's GIT_HASH; the file is preprocessed).
#define YREC_VERSION '202602'

#ifndef GIT_HASH
#define GIT_HASH 'unknown'
#endif

module yrec_version
      implicit none
      private
      public :: set_version_strings

contains

subroutine set_version_strings()
      use physical_constants
      implicit none
      yrec_version_string = YREC_VERSION
      git_hash_string = GIT_HASH
! A build outside a git checkout can define GIT_HASH as an empty (or
! '-dirty'-only) string, which the ifndef above cannot catch; an empty
! hash would make a zero-width A00 edit descriptor in the run-log banner.
      if (len_trim(git_hash_string) == 0) git_hash_string = 'unknown'
end subroutine set_version_strings

end module yrec_version
