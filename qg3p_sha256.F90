! Streaming SHA-256 for geometry identities. Use 64-bit intermediates so
! 32-bit word addition wraps explicitly without signed integer overflow.
module qg_sha256
  implicit none
  private
  public :: sha256, hash_line, hash_hex

  integer, parameter :: wide = selected_int_kind(18)
  integer(wide), parameter :: word_mask = int(z'FFFFFFFF', wide)
  integer(wide), parameter :: initial(8) = [ &
    int(z'6A09E667',wide), int(z'BB67AE85',wide), int(z'3C6EF372',wide), int(z'A54FF53A',wide), &
    int(z'510E527F',wide), int(z'9B05688C',wide), int(z'1F83D9AB',wide), int(z'5BE0CD19',wide) ]
  integer(wide), parameter :: constants(64) = [ &
    int(z'428A2F98',wide), int(z'71374491',wide), int(z'B5C0FBCF',wide), int(z'E9B5DBA5',wide), &
    int(z'3956C25B',wide), int(z'59F111F1',wide), int(z'923F82A4',wide), int(z'AB1C5ED5',wide), &
    int(z'D807AA98',wide), int(z'12835B01',wide), int(z'243185BE',wide), int(z'550C7DC3',wide), &
    int(z'72BE5D74',wide), int(z'80DEB1FE',wide), int(z'9BDC06A7',wide), int(z'C19BF174',wide), &
    int(z'E49B69C1',wide), int(z'EFBE4786',wide), int(z'0FC19DC6',wide), int(z'240CA1CC',wide), &
    int(z'2DE92C6F',wide), int(z'4A7484AA',wide), int(z'5CB0A9DC',wide), int(z'76F988DA',wide), &
    int(z'983E5152',wide), int(z'A831C66D',wide), int(z'B00327C8',wide), int(z'BF597FC7',wide), &
    int(z'C6E00BF3',wide), int(z'D5A79147',wide), int(z'06CA6351',wide), int(z'14292967',wide), &
    int(z'27B70A85',wide), int(z'2E1B2138',wide), int(z'4D2C6DFC',wide), int(z'53380D13',wide), &
    int(z'650A7354',wide), int(z'766A0ABB',wide), int(z'81C2C92E',wide), int(z'92722C85',wide), &
    int(z'A2BFE8A1',wide), int(z'A81A664B',wide), int(z'C24B8B70',wide), int(z'C76C51A3',wide), &
    int(z'D192E819',wide), int(z'D6990624',wide), int(z'F40E3585',wide), int(z'106AA070',wide), &
    int(z'19A4C116',wide), int(z'1E376C08',wide), int(z'2748774C',wide), int(z'34B0BCB5',wide), &
    int(z'391C0CB3',wide), int(z'4ED8AA4A',wide), int(z'5B9CCA4F',wide), int(z'682E6FF3',wide), &
    int(z'748F82EE',wide), int(z'78A5636F',wide), int(z'84C87814',wide), int(z'8CC70208',wide), &
    int(z'90BEFFFA',wide), int(z'A4506CEB',wide), int(z'BEF9A3F7',wide), int(z'C67178F2',wide) ]

  type sha256
    private
    integer(wide) :: digest(8) = initial
    integer(wide) :: bytes = 0_wide
    integer(wide) :: block(64)
    integer :: used = 0
  end type sha256

contains

  pure function rotate(word, n) result(value)
    integer(wide), intent(in) :: word
    integer, intent(in) :: n
    integer(wide) :: value

    value = ishftc(word, -n, 32)
  end function rotate

  subroutine compress(hash)
    type(sha256), intent(inout) :: hash
    integer(wide) :: schedule(64), state(8), a, b, e, f, g, first, second, s0, s1, majority, choice
    integer :: i, j

    do i = 1,16
      j = 4*(i-1)
      schedule(i) = ishft(hash%block(j+1),24) + ishft(hash%block(j+2),16) + &
                    ishft(hash%block(j+3),8) + hash%block(j+4)
    end do
    do i = 17,64
      a = schedule(i-15)
      b = schedule(i-2)
      s0 = ieor(ieor(rotate(a,7),rotate(a,18)),ishft(a,-3))
      s1 = ieor(ieor(rotate(b,17),rotate(b,19)),ishft(b,-10))
      schedule(i) = iand(schedule(i-16) + s0 + schedule(i-7) + s1, word_mask)
    end do
    state = hash%digest
    do i = 1,64
      a = state(1); b = state(2)
      e = state(5); f = state(6); g = state(7)
      choice = ieor(iand(e,f),iand(not(e),g))
      majority = ieor(iand(a,b),iand(ieor(a,b),state(3)))
      s0 = ieor(ieor(rotate(a,2),rotate(a,13)),rotate(a,22))
      s1 = ieor(ieor(rotate(e,6),rotate(e,11)),rotate(e,25))
      first = state(8) + s1 + choice + constants(i) + schedule(i)
      second = s0 + majority
      state(2:8) = state(1:7)
      state(1) = iand(first + second, word_mask)
      state(5) = iand(state(5) + first, word_mask)
    end do
    hash%digest = iand(hash%digest + state, word_mask)
  end subroutine compress

  subroutine hash_byte(hash, byte)
    type(sha256), intent(inout) :: hash
    integer, intent(in) :: byte

    hash%used = hash%used + 1
    hash%bytes = hash%bytes + 1
    hash%block(hash%used) = byte
    if (hash%used == 64) then
      call compress(hash)
      hash%used = 0
    end if
  end subroutine hash_byte

  subroutine hash_line(hash, line)
    type(sha256), intent(inout) :: hash
    character(*), intent(in) :: line
    integer :: i

    do i = 1,len_trim(line)
      call hash_byte(hash, iachar(line(i:i)))
    end do
    call hash_byte(hash, 10)
  end subroutine hash_line

  function hash_hex(hash) result(hex)
    type(sha256), intent(in) :: hash
    character(64) :: hex
    type(sha256) :: padded
    integer(wide) :: bits
    integer :: i

    padded = hash
    bits = hash%bytes * 8
    call hash_byte(padded, 128)
    do while (padded%used /= 56)
      call hash_byte(padded, 0)
    end do
    do i = 56,0,-8
      call hash_byte(padded, int(iand(ishft(bits,-i),255_wide)))
    end do
    write(hex, '(8z8.8)') padded%digest
    do i = 1,len(hex)
      if (hex(i:i) >= 'A' .and. hex(i:i) <= 'F') hex(i:i) = achar(iachar(hex(i:i)) + 32)
    end do
  end function hash_hex

end module qg_sha256
