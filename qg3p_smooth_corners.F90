! Polynomial diagnostics on smooth corner data. Singular wall fluxes require
! the separate coupled weak operator; this point plan is not that closure.
module qg_smooth_corners
  use lapack_interfaces, only: dgesvd
  use qg_diagnostics, only: BC_FREE_SLIP, BC_PARTIAL_SLIP, BC_NO_SLIP
  use qg_irregular_parallel, only: distributed_geometry, donor_plan, build_donor_plan, gather_donors, &
                                   half_coordinates, make_sector_moments, cell_mean_moments, broadcast_error
  use qg_geometry, only: fluid_node
  use qg_moments, only: n_moments, x_power, y_power, monomial_index, segment_moments, point_moments
  use qg_reconstruction, only: visible_donor, minimum_norm_weights
  use qg_parameters, only: real_8
  implicit none
  private
  public :: smooth_corner_plan, setup_smooth_corners, reconstruct_smooth_corners, smooth_corner_weights
  public :: corner_psi_mean, corner_omega_point, corner_omega_mean, corner_gamma_first, corner_gamma_second

  integer, parameter :: terms = n_moments-1, radius_limit = 4, max_sites = (2*radius_limit+1)**2
  integer, parameter :: corner_psi_mean=1, corner_omega_point=2, corner_omega_mean=3, corner_gamma_first=4, corner_gamma_second=5
  integer, parameter :: quantities = 5
  real(real_8), parameter :: ray(2,0:7) = reshape([1,0,1,1,0,1,-1,1,-1,0,-1,-1,0,-1,1,-1],[2,8])

  type smooth_corner_row
    integer :: i,j,component,start,count
    real(real_8), allocatable :: weights(:,:)
  end type smooth_corner_row

  type smooth_corner_plan
    private
    type(donor_plan) :: donors
    type(smooth_corner_row), allocatable :: rows(:)
    integer :: donor_count
  end type smooth_corner_plan

contains

  subroutine smooth_corner_weights(mask, slip, length, points, weights, accepted)
    integer, intent(in) :: mask, slip
    real(real_8), intent(in) :: length, points(:,:)
    real(real_8), intent(out) :: weights(:,:)
    logical, intent(out) :: accepted
    integer :: nr,r,first,last,m,k,rank,ndof,info,n,branch
    real(real_8) :: t(2,2), outward(2,2), unit_length, ell_scale, value(n_moments), means(n_moments), line(n_moments)
    real(real_8) :: row(terms), derivative(terms), second(terms), constraint(16,terms), vt(terms,terms)
    real(real_8) :: singular(terms), work(256), u(1,1), basis(terms,terms), func(terms,quantities)
    real(real_8) :: design(size(points,2),terms), distance_weight(size(points,2)), scale
    integer, parameter :: degree(terms) = x_power(2:)+y_power(2:)

    weights = 0.0D0
    accepted = .false.
    if (popcnt(mask) == 0 .or. popcnt(mask) == 8) return
    do first = 0,7
      if (btest(mask,first) .and. .not. btest(mask,modulo(first-1,8))) exit
    end do
    last = modulo(first+popcnt(mask),8)
    t(:,1) = ray(:,first)
    t(:,2) = ray(:,last)
    outward(:,1) = [t(2,1),-t(1,1)]/norm2(t(:,1))
    outward(:,2) = [-t(2,2),t(1,2)]/norm2(t(:,2))
    nr = 0
    ell_scale = max(1.0D0,length)
    do branch = 1,2
      value = point_moments(t(:,branch))
      do k = 1,4
        call append(where_degree(value(2:),k))
      end do
      derivative = normal_derivative(value,-outward(:,branch))
      second = normal_second(value,-outward(:,branch))
      if (slip == BC_NO_SLIP) then
        do k = 0,3
          call append(where_degree(derivative,k+1))
        end do
      else
        do k = 0,2
          row = where_degree(second,k+2)
          if (slip == BC_PARTIAL_SLIP) row = (length/ell_scale)*row-where_degree(derivative,k+1)/ell_scale
          call append(row)
        end do
      end if
    end do
    call dgesvd('N','A',nr,terms,constraint,16,singular,u,1,vt,terms,work,size(work),info)
    if (info /= 0) return
    rank = count(singular(:min(nr,terms)) > 64*epsilon(1.0D0)*max(nr,terms)*singular(1))
    ndof = terms-rank
    if (ndof == 0) then
      accepted = .true.
      return
    end if
    n = size(points,2)
    if (n < ndof) return
    basis(:,:ndof) = transpose(vt(rank+1:,:))
    do k = 1,n
      value = point_moments(points(:,k))
      design(k,:) = value(2:)
      distance_weight(k) = 1.0D0/(1.0D0+sum(points(:,k)**2))
    end do
    scale = norm2(design*spread(distance_weight,2,terms))
    means = cell_mean_moments(mask,make_sector_moments())
    func(:,corner_psi_mean) = means(2:)
    func(:,corner_omega_point) = laplacian_row(point_moments([0.0D0,0.0D0]))
    func(:,corner_omega_mean) = laplacian_row(means)
    do branch = 1,2
      line = segment_moments([0.0D0,0.0D0],0.5D0*t(:,branch))
      unit_length = norm2(t(:,branch))/2
      func(:,3+branch) = unit_length*normal_derivative(line,outward(:,branch))
    end do
    do k = 1,quantities
      call minimum_norm_weights(matmul(design,basis(:,:ndof)),distance_weight, &
           matmul(transpose(basis(:,:ndof)),func(:,k)),weights(:,k),accepted,rank_scale=scale)
      if (.not. accepted) return
    end do

  contains
    pure function where_degree(values,d) result(selected)
      real(real_8), intent(in) :: values(terms)
      integer, intent(in) :: d
      real(real_8) :: selected(terms)
      selected = merge(values,0.0D0,degree==d)
    end function where_degree

    subroutine append(values)
      real(real_8), intent(in) :: values(terms)
      real(real_8) :: row_norm

      row_norm = norm2(values)
      if (row_norm == 0.0D0) return
      nr = nr+1
      constraint(nr,:) = values/row_norm
    end subroutine append
  end subroutine smooth_corner_weights

  pure function normal_derivative(values,normal) result(row)
    real(real_8), intent(in) :: values(n_moments),normal(2)
    real(real_8) :: row(terms)
    integer :: m,a,b

    row = 0.0D0
    do m = 2,n_moments
      a = x_power(m)
      b = y_power(m)
      if (a > 0) row(m-1) = a*normal(1)*values(monomial_index(a-1,b))
      if (b > 0) row(m-1) = row(m-1)+b*normal(2)*values(monomial_index(a,b-1))
    end do
  end function normal_derivative

  pure function normal_second(values,normal) result(row)
    real(real_8), intent(in) :: values(n_moments),normal(2)
    real(real_8) :: row(terms)
    integer :: m,a,b

    row = 0.0D0
    do m = 4,n_moments
      a = x_power(m)
      b = y_power(m)
      if (a >= 2) row(m-1) = a*(a-1)*normal(1)**2*values(monomial_index(a-2,b))
      if (b >= 2) row(m-1) = row(m-1)+b*(b-1)*normal(2)**2*values(monomial_index(a,b-2))
      if (a > 0 .and. b > 0) row(m-1) = row(m-1)+2*a*b*product(normal)*values(monomial_index(a-1,b-1))
    end do
  end function normal_second

  pure function laplacian_row(values) result(row)
    real(real_8), intent(in) :: values(n_moments)
    real(real_8) :: row(terms)
    integer :: m,a,b

    row = 0.0D0
    do m = 4,n_moments
      a = x_power(m)
      b = y_power(m)
      if (a >= 2) row(m-1) = a*(a-1)*values(monomial_index(a-2,b))
      if (b >= 2) row(m-1) = row(m-1)+b*(b-1)*values(monomial_index(a,b-2))
    end do
  end function laplacian_row

  subroutine setup_smooth_corners(geom,slip,length,plan,error)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: slip
    real(real_8), intent(in) :: length
    type(smooth_corner_plan), intent(out) :: plan
    character(*), intent(out) :: error
    type(smooth_corner_row), allocatable :: rows(:)
    integer, allocatable :: requested(:,:)
    integer :: n,nd,edge,gi,gj,i,j,k,r,mask,origin(2),site(2,max_sites),ni,nj
    real(real_8) :: coordinates(2,max_sites), weights(max_sites,quantities)
    logical :: accepted

    error = ''
    if (slip < min(BC_FREE_SLIP,BC_NO_SLIP,BC_PARTIAL_SLIP) .or. slip > max(BC_FREE_SLIP,BC_NO_SLIP,BC_PARTIAL_SLIP)) &
      error = 'unsupported smooth corner slip selector'
    if (slip == BC_PARTIAL_SLIP .and. length <= 0.0D0) error = 'partial-slip smooth corner requires a positive length'
    allocate(rows(size(geom%walls)),requested(2,max_sites*size(geom%walls)))
    ni = size(geom%owner,1)
    nj = size(geom%owner,2)
    n = 0
    nd = 0
    do edge = 1,size(geom%walls)
      if (len_trim(error) /= 0) exit
      origin = geom%walls(edge)%endpoints(:,1)
      gi = ni-origin(2)/2
      gj = origin(1)/2+1
      mask = geom%global%wedge(gi,gj)
      if (popcnt(mask) == 4) cycle
      do r = 1,radius_limit
        k = 0
        do j = max(1,gj-r),min(nj,gj+r)
          do i = max(1,gi-r),min(ni,gi+r)
            if (geom%global%kind(i,j) /= fluid_node) cycle
            if (.not. visible_donor(geom,2*origin,2*half_coordinates(i,j,ni))) cycle
            k = k+1
            site(:,k) = [i,j]
            coordinates(:,k) = [j-gj,gi-i]
          end do
        end do
        call smooth_corner_weights(mask,slip,length,coordinates(:,:k),weights(:k,:),accepted)
        if (accepted) exit
      end do
      if (.not. accepted) then
        write(error,'(a,2(1x,i0))') 'insufficient strict-fluid donors for smooth corner at i, j',gi,gj
        exit
      end if
      n = n+1
      rows(n)%i = gi-geom%i0
      rows(n)%j = gj-geom%j0
      rows(n)%component = geom%walls(edge)%component
      rows(n)%start = nd
      rows(n)%count = k
      rows(n)%weights = weights(:k,:)
      requested(:,nd+1:nd+k) = site(:,:k)
      nd = nd+k
    end do
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    plan%rows = rows(:n)
    plan%donor_count = nd
    call build_donor_plan(geom,requested(:,:nd),plan%donors,error)
  end subroutine setup_smooth_corners

  ! Write only smooth corner sites. Column one is psi cell mean, then
  ! point/mean vorticity and the two oriented half-wall integrals.
  subroutine reconstruct_smooth_corners(plan,psi,wall,result)
    type(smooth_corner_plan), intent(in) :: plan
    real(real_8), intent(in) :: psi(:,:,:),wall(0:,:)
    real(real_8), intent(inout) :: result(:,:,:,:) ! (i,j,layer,quantity)
    real(real_8) :: donors(size(psi,3),plan%donor_count), difference(size(psi,3)), c(size(psi,3))
    integer :: n,k,m

    call gather_donors(plan%donors,psi,donors)
    do n = 1,size(plan%rows)
      associate(row => plan%rows(n))
        c = wall(row%component,:)
        result(row%i,row%j,:,:) = 0.0D0
        result(row%i,row%j,:,corner_psi_mean) = c
        do k = 1,row%count
          difference = donors(:,row%start+k)-c
          do m = 1,quantities
            result(row%i,row%j,:,m) = result(row%i,row%j,:,m)+row%weights(k,m)*difference
          end do
        end do
      end associate
    end do
  end subroutine reconstruct_smooth_corners

end module qg_smooth_corners
