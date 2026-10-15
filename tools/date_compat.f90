! SPDX-License-Identifier: GPL-2.0-or-later
! Compatibility implementation for the legacy Vanderbilt generator.
subroutine idate(idmy)
  implicit none
  integer, intent(out) :: idmy(3)
  integer :: values(8)
  call date_and_time(values=values)
  idmy = [values(3), values(2), values(1)]
end subroutine idate
