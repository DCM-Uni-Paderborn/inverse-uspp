! SPDX-License-Identifier: GPL-2.0-or-later
! Use Quantum ESPRESSO's readers, augmentation expansion and UPF v2 writer.
! Based on the conversion sequence in QE 6.3MaX upftools/upf2upf2.f90.
program convert_gbrv
  use pseudo_types, only: pseudo_upf, nullify_pseudo_upf
  use radial_grids, only: radial_grid_type, nullify_radial_grid
  use read_upf_v1_module, only: read_upf_v1
  use upf_to_internal, only: set_pseudo_upf
  use write_upf_module, only: write_upf
  implicit none
  type(pseudo_upf) :: upf
  type(radial_grid_type), target :: grid
  character(2048) :: filename
  integer :: unit,status
  integer, external :: atomic_number
  call get_command_argument(1,filename)
  if(len_trim(filename)==0) error stop 'Supply an original GBRV UPF v1 file'
  call nullify_pseudo_upf(upf)
  call nullify_radial_grid(grid)
  open(newunit=unit,file=trim(filename),status='old')
  call read_upf_v1(unit,upf,grid,status)
  close(unit)
  if(status/=0) error stop 'Cannot read original UPF v1'
  if(.not.upf%tvanp.or.upf%tpawp.or.upf%has_so) error stop 'Expected scalar USPP'
  ! This is QE's own polynomial-to-angular-table conversion, without resampling.
  call set_pseudo_upf(1,upf)
  upf%q_with_l=.true.
  upf%nqf=0
  upf%nv='2.0.1'
  upf%rmax=upf%r(upf%mesh)
  upf%zmesh=atomic_number(upf%psd)
  upf%rel='scalar'
  if(upf%r(1)>1.d-16) then
    upf%dx=log(upf%rmax/upf%r(1))/(upf%mesh-1)
    upf%xmin=log(upf%r(1)*upf%zmesh)
  endif
  if(.not.associated(upf%epseu)) then
    allocate(upf%epseu(upf%nwfc))
    upf%epseu=0
  endif
  allocate(upf%nchi(upf%nwfc),upf%rcut_chi(upf%nwfc),upf%rcutus_chi(upf%nwfc))
  upf%nchi=0
  if(associated(upf%nn)) upf%nchi=upf%nn(1:upf%nwfc)
  upf%rcut_chi=upf%rcut(1:upf%nwfc)
  upf%rcutus_chi=upf%rcutus(1:upf%nwfc)
  call write_upf(filename=trim(filename)//'.UPF',upf=upf,schema='v2')
end program convert_gbrv
