!> Program to track a Bmad beam and write projected RMS phase-space values
!> at every saved point along the beamline.
!>
!> The saved RMS values are computed from the bunch covariance matrix
!> \f$\Sigma\f$ returned by `calc_bunch_params` via the internal
!> `bunch_track_struct` save path:
!> \f[
!>   \mathrm{rms}(u_i) = \sqrt{\Sigma_{ii}},
!> \f]
!> with \f$u = (x, p_x, y, p_y, z, p_z)\f$ in Bmad's beam coordinates.
program bmad_beam_rms_track

use bmad
use beam_mod
use beam_file_io

implicit none

type (lat_struct)               :: lat
type (beam_struct), target      :: beam
type (beam_init_struct)         :: beam_init
type (ele_struct), pointer      :: ele1, ele2
type (bunch_track_struct), allocatable :: bunch_tracks(:)

character(256) :: in_filename
character(256) :: lat_filename
character(256) :: outfile_name
character(256) :: rms_outfile_name

integer :: n_arg, ran_seed
integer :: i_bunch
logical :: err
real(rp) :: ds_save

namelist /beam_track_params/ lat_filename, beam_init, outfile_name, rms_outfile_name, ran_seed, ds_save

call read_input_namelist(in_filename, lat_filename, beam_init, outfile_name, rms_outfile_name, ran_seed, ds_save)

call ran_seed_put(ran_seed)
call bmad_parser(lat_filename, lat)

ele1 => lat%ele(0)
call init_beam_distribution(ele1, lat%param, beam_init, beam)

allocate(bunch_tracks(size(beam%bunch)))
do i_bunch = 1, size(bunch_tracks)
  bunch_tracks(i_bunch)%ds_save = ds_save
  bunch_tracks(i_bunch)%n_pt = -1
end do

ele2 => lat%ele(lat%n_ele_track)
call track_beam(lat, beam, ele1, ele2, err, bunch_tracks = bunch_tracks)
if (err) stop 1

call write_rms_track_file(rms_outfile_name, lat, bunch_tracks)
call write_beam_file(outfile_name, beam, lat = lat)

contains

  !> Read the beam tracking namelist with defaults compatible with
  !> `beam_track_example`.
  subroutine read_input_namelist(in_file, lat_file, beam_init_out, beam_file, rms_file, seed, ds_save_out)
    character(*), intent(out)        :: in_file
    character(*), intent(out)        :: lat_file
    type (beam_init_struct), intent(out) :: beam_init_out
    character(*), intent(out)        :: beam_file
    character(*), intent(out)        :: rms_file
    integer, intent(out)             :: seed
    real(rp), intent(out)            :: ds_save_out

    integer :: iu
    integer :: n_arg_local

    n_arg_local = command_argument_count()
    if (n_arg_local > 1) then
      print *, 'Usage: bmad_beam_rms_track <input_file>'
      print *, 'Default: <input_file> = beam_track_example.in'
      stop 1
    end if

    in_file = 'beam_track_example.in'
    if (n_arg_local == 1) call get_command_argument(1, in_file)

    lat_file = 'lat.bmad'
    beam_file = 'beam.dat'
    rms_file = 'beam_rms.dat'
    beam_init_out%n_particle = 1
    beam_init_out%n_bunch = 1
    seed = 0
    ds_save_out = 0.0_rp

    open(newunit = iu, file = trim(in_file), status = 'old')
    read(iu, nml = beam_track_params)
    close(iu)
  end subroutine read_input_namelist

  !> Write a plain-text table of bunch RMS values.
  subroutine write_rms_track_file(file_name, lat_in, bunch_tracks_in)
    character(*), intent(in)                    :: file_name
    type (lat_struct), intent(in)               :: lat_in
    type (bunch_track_struct), intent(in)       :: bunch_tracks_in(:)

    integer :: iu
    integer :: ix_bunch, ix_pt, ix_ele
    character(40) :: ele_name
    character(20) :: ele_type
    character(16) :: loc_name

    open(newunit = iu, file = trim(file_name), status = 'replace')
    write(iu, '(a)') '# Bmad bunch RMS track'
    write(iu, '(a)') '# rms_i = sqrt(max(0, sigma_ii)) for (x, px, y, py, z, pz)'
    write(iu, '(a)') '# bunch s t ix_ele ele_name ele_type location n_live rms_x rms_px rms_y rms_py rms_z rms_pz'

    do ix_bunch = 1, size(bunch_tracks_in)
      if (.not. allocated(bunch_tracks_in(ix_bunch)%pt)) cycle
      if (bunch_tracks_in(ix_bunch)%n_pt < 0) cycle

      do ix_pt = 0, bunch_tracks_in(ix_bunch)%n_pt
        ix_ele = bunch_tracks_in(ix_bunch)%pt(ix_pt)%ix_ele
        ele_name = ''
        ele_type = ''
        if (ix_ele >= 0 .and. ix_ele <= lat_in%n_ele_track) then
          ele_name = trim(lat_in%ele(ix_ele)%name)
          ele_type = trim(key_name(lat_in%ele(ix_ele)%key))
        end if
        loc_name = trim(location_name(max(0, min(ubound(location_name, 1), bunch_tracks_in(ix_bunch)%pt(ix_pt)%location))))

        write(iu, '(i4,1x,f14.8,1x,es18.10,1x,i6,1x,a40,1x,a20,1x,a16,1x,i8,6(1x,es18.10))') &
          ix_bunch, &
          bunch_tracks_in(ix_bunch)%pt(ix_pt)%s, &
          bunch_tracks_in(ix_bunch)%pt(ix_pt)%t, &
          ix_ele, &
          trim(ele_name), &
          trim(ele_type), &
          trim(loc_name), &
          bunch_tracks_in(ix_bunch)%pt(ix_pt)%n_particle_live, &
          rms_from_sigma(bunch_tracks_in(ix_bunch)%pt(ix_pt)%sigma(1,1)), &
          rms_from_sigma(bunch_tracks_in(ix_bunch)%pt(ix_pt)%sigma(2,2)), &
          rms_from_sigma(bunch_tracks_in(ix_bunch)%pt(ix_pt)%sigma(3,3)), &
          rms_from_sigma(bunch_tracks_in(ix_bunch)%pt(ix_pt)%sigma(4,4)), &
          rms_from_sigma(bunch_tracks_in(ix_bunch)%pt(ix_pt)%sigma(5,5)), &
          rms_from_sigma(bunch_tracks_in(ix_bunch)%pt(ix_pt)%sigma(6,6))
      end do
    end do

    close(iu)
  end subroutine write_rms_track_file

  !> Guard the RMS extraction against tiny negative roundoff.
  pure real(rp) function rms_from_sigma(sigma_ii) result(rms_value)
    real(rp), intent(in) :: sigma_ii
    rms_value = sqrt(max(0.0_rp, sigma_ii))
  end function rms_from_sigma

end program bmad_beam_rms_track
