!----------------------------------------------------------------------
! yrec_version
!----------------------------------------------------------------------
! The code's version stamp: yrec_version_string (the release tag) and
! git_hash_string (the build's commit, from the Makefile's GIT_HASH; the
! file is preprocessed), filled by set_version_strings at the start of
! run_yrec, printed in the controls read's banner (io/controls_reader.f90,
! io/derived_controls.f90) and in the history header. A new commit
! recompiles this file only: the module's interface does not change, so
! its users are not rebuilt.
!
! Sources
!   The two strings were physical_constants' until 2026-10-29, and the
!   F77's common/version/ (yrecver, githash) before that.
#define YREC_VERSION '202602'

#ifndef GIT_HASH
#define GIT_HASH 'unknown'
#endif

module yrec_version
      implicit none
      private
      public :: set_version_strings, yrec_version_string, git_hash_string

      character(len=10) :: yrec_version_string = ' '
      character(len=20) :: git_hash_string = ' '

contains

subroutine set_version_strings()
      yrec_version_string = YREC_VERSION
      git_hash_string = GIT_HASH
! A build outside a git checkout can define GIT_HASH as an empty (or
! '-dirty'-only) string, which the ifndef above cannot catch; an empty
! hash would make a zero-width A00 edit descriptor in the run-log banner.
      if (len_trim(git_hash_string) == 0) git_hash_string = 'unknown'
end subroutine set_version_strings

end module yrec_version
