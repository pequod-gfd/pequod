program test_irregular_geometry_reader
  use qg_geometry, only: irregular_geometry, read_geometry
  use qg_parameters, only: real_8
  use qg_sha256, only: sha256, hash_line, hash_hex
  implicit none
  type(irregular_geometry) :: geometry
  type(sha256) :: hash
  character(1024) :: path, value, error
  integer :: ni, nj
  real(real_8) :: h
  integer, parameter :: wide = selected_int_kind(18)

  if (command_argument_count() == 0) then
    if (hash_hex(hash) /= 'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855') stop 1
    call hash_line(hash, 'abc')
    if (hash_hex(hash) /= 'edeaaff3f1774ad2888673770c6d64097e391bc362d7d6fb34982ddf0efd18cb') stop 1
    print '(a)', '[SHA256 known digests] : PASS'
    stop
  end if

  if (command_argument_count() /= 4) stop 2
  call get_command_argument(1,path)
  call get_command_argument(2,value)
  read(value,*) ni
  call get_command_argument(3,value)
  read(value,*) nj
  call get_command_argument(4,value)
  read(value,*) h
  call read_geometry(trim(path),ni,nj,h,geometry,error)
  if (len_trim(error) /= 0) then
    print '(a)', trim(error)
    stop 1
  end if
  print '(a,2(1x,i0))', geometry%digest, size(geometry%loops), sum(int(popcnt(geometry%wedge),wide))
end program test_irregular_geometry_reader
