!----------------------------------------------------------------------
! fermi_table_reader
!----------------------------------------------------------------------
! The Fermi-Dirac table's reader, read_fermi_table (eos_init calls it):
! the file's I/O, kept apart from the Yale EOS that interpolates the
! table (yale_eos prints nothing). The table lives in yale_eos_tables.
module fermi_table_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use yale_eos_tables, only: yale, fermi_n_x, fermi_n_quantities, fermi_n_cards, &
           fermi_n_bins, fermi_bins_per_unit_x
      use eos_messages, only: eos_report, message_length, unwritten_line
      implicit none
      private
      public :: read_fermi_table

contains

!----------------------------------------------------------------------
! read_fermi_table
!----------------------------------------------------------------------
! Reads the Fermi-Dirac integral table of the degenerate electron gas
! (FERMI.TAB: 43 abscissas, 43 degeneracies, 860 cards of 5 values)
! into yale, and bins the abscissas into fermi_table_x_lookup, the
! interval index compute_fully_ionized_eos looks up. A table that cannot
! be opened, a read that fails (a truncated or corrupt table), a card
! whose (abscissa, row) is not the next in the file's order (abscissa
! first), a non-increasing abscissa and abscissae that span more than
! the lookup's bins are ierr = 1, with a message: a card's own indices
! are not trusted as subscripts.
subroutine read_fermi_table(fermi_table_path, ierr)
      character(len=*), intent(in) :: fermi_table_path
      integer, intent(out) :: ierr
      character(len=message_length) :: eos_lines(4)
      integer :: fermi_unit, grid_idx, iden_idx, card_idx, ios
      integer :: bin_start, bin_width, bin_end
      real(dp) :: card_values(fermi_n_quantities)
      logical :: misplaced_card

      ierr = 0
      open(newunit=fermi_unit, file=fermi_table_path, form='FORMATTED', status='OLD', &
           iostat=ierr)
      if (ierr /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: cannot open the Fermi-Dirac table ', trim(fermi_table_path)
         call eos_report(eos_lines)
         ierr = 1
         return
      end if
      read(fermi_unit,150,iostat=ios)  (yale%fermi_table_x_grid(grid_idx), grid_idx = 1,fermi_n_x)
      if (ios == 0) read(fermi_unit,150,iostat=ios)  (yale%fermi_table_eta(grid_idx), grid_idx = 1,fermi_n_x)
  150 format(8x,9f8.4)
! the cards, one per (abscissa, row), the abscissa varying first
      misplaced_card = .false.
      card_idx = 0
      do while (ios == 0 .and. card_idx < fermi_n_cards)
         card_idx = card_idx + 1
         read(fermi_unit,160,iostat=ios) grid_idx, iden_idx, card_values
  160    format(2i5,5f12.7)
         if (ios /= 0) exit
         misplaced_card = grid_idx /= mod(card_idx - 1, fermi_n_x) + 1 .or. &
              iden_idx /= (card_idx - 1)/fermi_n_x + 1
         if (misplaced_card) exit
         yale%fermi_table_data(:,grid_idx,iden_idx) = card_values
      end do
      if (ios /= 0) then
         eos_lines = unwritten_line
         write(eos_lines,*) 'eos_init: the Fermi-Dirac table ends early or is corrupt: ', &
              trim(fermi_table_path)
         call eos_report(eos_lines)
         ierr = 1
      else if (misplaced_card) then
         eos_lines = unwritten_line
         write(eos_lines,'(a,i0,a,i0,a,i0,a,a)') 'read_fermi_table: card ', card_idx, &
              ' is for abscissa ', grid_idx, ' and row ', iden_idx, &
              ', not the next in the table: ', trim(fermi_table_path)
         call eos_report(eos_lines)
         ierr = 1
      else
! the X lookup: each abscissa's bins name it; a non-increasing abscissa,
! or abscissae that span more than the bins, are a corrupt table
         bin_start = 1
         do grid_idx = 1, fermi_n_x - 2
            if (yale%fermi_table_x_grid(grid_idx+1) <= yale%fermi_table_x_grid(grid_idx)) then
               eos_lines = unwritten_line
               write(eos_lines,'(a,i0,a,a)') 'read_fermi_table: abscissa ', grid_idx, &
                    ' does not increase: ', trim(fermi_table_path)
               call eos_report(eos_lines)
               ierr = 1
               exit
            end if
            bin_width = int((yale%fermi_table_x_grid(grid_idx+1) - yale%fermi_table_x_grid(grid_idx))* &
                 fermi_bins_per_unit_x + 0.10d0)
            bin_end = bin_start + bin_width - 1
            if (grid_idx == fermi_n_x - 2) bin_end = fermi_n_bins
            if (bin_end > fermi_n_bins) then
               eos_lines = unwritten_line
               write(eos_lines,'(a,i0,a,a)') 'read_fermi_table: the abscissae span more than the ', &
                    fermi_n_bins, ' bins of the X lookup: ', trim(fermi_table_path)
               call eos_report(eos_lines)
               ierr = 1
               exit
            end if
            do iden_idx = bin_start, bin_end
               yale%fermi_table_x_lookup(iden_idx) = grid_idx
            end do
            bin_start = bin_end + 1
         end do
      end if
      close(fermi_unit)
end subroutine read_fermi_table

end module fermi_table_reader
