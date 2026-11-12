! Local polynomial functionals of wet-cell means or accessible nodal data.
module qg_reconstruction
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use lapack_interfaces, only: dgesvd
  use qg_geometry, only: solid_node, wall_node
  use qg_irregular_parallel, only: distributed_geometry, donor_plan, build_donor_plan, gather_donors, &
                                   make_sector_moments, cell_mean_moments, half_coordinates, broadcast_error
  use qg_moments, only: n_moments, point_moments, translate_moments
  use qg_parameters, only: real_8
  implicit none
  private
  public :: reconstruction_plan, cubic_terms, build_reconstruction, evaluate_reconstruction
  public :: minimum_norm_weights, visible_donor

  integer, parameter :: cubic_terms = 10, max_radius = 4, max_sites = (2*max_radius+1)**2
  integer, parameter :: wide = selected_int_kind(18)
  integer, parameter :: di(0:7) = [0,-1,-1,-1,0,1,1,1], dj(0:7) = [1,1,0,-1,-1,-1,0,1]

  type reconstruction_row
    integer, allocatable :: sample(:) ! positive donor index or -(wall component+1)
    real(real_8), allocatable :: weight(:)
    real(real_8) :: constant
  end type reconstruction_row

  type reconstruction_plan
    private
    type(reconstruction_row), allocatable :: rows(:)
    type(donor_plan) :: donors
    integer :: donor_count
  end type reconstruction_plan

contains

  ! Return one linear functional exact on the chosen polynomial space. The
  ! absolute rank floor uses the unrestricted design scale where supplied.
  subroutine minimum_norm_weights(design, distance_weight, target, weights, accepted, rank_scale)
    real(real_8), intent(in) :: design(:,:), distance_weight(:), target(:)
    real(real_8), intent(out) :: weights(:)
    logical, intent(out) :: accepted
    real(real_8), optional, intent(in) :: rank_scale
    integer :: m,n,info
    real(real_8) :: a(size(design,1),size(design,2)), singular(size(design,2)), &
                    u(size(design,1),size(design,2)), vt(size(design,2),size(design,2)), &
                    work(8*max(size(design,1),size(design,2))), scale

    m = size(design,1)
    n = size(design,2)
    accepted = .false.
    if (m < n .or. n == 0) return
    a = design*spread(distance_weight,2,n)
    scale = norm2(a)
    if (present(rank_scale)) scale = max(scale,rank_scale)
    call dgesvd('S','A',m,n,a,m,singular,u,m,vt,n,work,size(work),info)
    if (info /= 0 .or. .not. all(ieee_is_finite(singular))) return
    if (singular(n) <= max(singular(1)*1.0D-6,64*epsilon(scale)*max(m,n)*scale)) return
    weights = distance_weight*matmul(u,matmul(vt,target)/singular)
    accepted = maxval(abs(matmul(transpose(design),weights)-target)) < 1.0D-10*max(1.0D0,maxval(abs(target)))
  end subroutine minimum_norm_weights

  ! Origins are in quarter-grid physical coordinates. Each sample is either
  ! a strict point (with exact C for wall samples) or a positive-area mean.
  subroutine build_reconstruction(geom, origins, targets, from_means, plan, error)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: origins(:,:)
    real(real_8), intent(in) :: targets(:,:)
    logical, intent(in) :: from_means
    type(reconstruction_plan), intent(out) :: plan
    character(*), intent(out) :: error
    integer, allocatable :: donor_index(:,:), requested(:,:)
    integer :: site,nsite,n,ni,nj,i,j,k,nd,point(2,max_sites), radius, ic,jc, kind, mask, p(2)
    real(real_8) :: design(max_sites,cubic_terms), weight(max_sites), attenuation(max_sites), offset(2)
    real(real_8) :: moment(n_moments), sector(n_moments,0:7)
    logical :: accepted

    nsite = size(origins,2)
    ni = size(geom%owner,1)
    nj = size(geom%owner,2)
    allocate(plan%rows(nsite),donor_index(ni,nj),requested(2,ni*nj))
    donor_index = 0
    nd = 0
    sector = make_sector_moments()
    error = ''
    if (size(origins,1) /= 2 .or. any(shape(targets) /= [cubic_terms,nsite])) &
      error = 'cubic reconstruction has incompatible site/target shape'
    do site = 1,nsite
      if (len_trim(error) /= 0) exit
      ic = ni-origins(2,site)/4
      jc = 1+origins(1,site)/4
      accepted = .false.
      do radius = 2,max_radius
        n = 0
        do i = max(1,ic-radius),min(ni,ic+radius)
          do j = max(1,jc-radius),min(nj,jc+radius)
            mask = geom%global%wedge(i,j)
            if (mask == 0) cycle
            p = 2*half_coordinates(i,j,ni)
            if (.not. visible_donor(geom,origins(:,site),p)) cycle
            n = n+1
            point(:,n) = [i,j]
            offset = 0.25D0*(p-origins(:,site))
            attenuation(n) = 1.0D0/(1.0D0+sum(offset**2))
            if (from_means) then
              moment = translate_moments(cell_mean_moments(mask,sector),offset)
            else
              moment = point_moments(offset)
            end if
            design(n,:) = moment(:cubic_terms)
          end do
        end do
        call minimum_norm_weights(design(:n,:),attenuation(:n),targets(:,site),weight(:n),accepted)
        if (accepted) exit
      end do
      if (.not. accepted) then
        write(error,'(a,2(1x,i0))') 'insufficient accessible cubic donors at quarter-grid x, y', origins(:,site)
        exit
      end if
      associate(row => plan%rows(site))
        row%constant = targets(1,site)
        row%weight = weight(:n)
        allocate(row%sample(n))
        do k = 1,n
          i = point(1,k)
          j = point(2,k)
          kind = geom%global%kind(i,j)
          if (kind == wall_node .and. .not. from_means) then
            row%sample(k) = -geom%global%component(i,j)-1
          else
            if (donor_index(i,j) == 0) then
              nd = nd+1
              donor_index(i,j) = nd
              requested(:,nd) = [i,j]
            end if
            row%sample(k) = donor_index(i,j)
          end if
        end do
      end associate
    end do
    call broadcast_error(geom%comm,geom%rank,geom%ranks,error)
    if (len_trim(error) /= 0) return
    plan%donor_count = nd
    call build_donor_plan(geom,requested(:,:nd),plan%donors,error,cell_means=from_means)
  end subroutine build_reconstruction

  subroutine evaluate_reconstruction(plan, field, values, wall_values)
    type(reconstruction_plan), intent(in) :: plan
    real(real_8), intent(in) :: field(:,:,:)
    real(real_8), intent(out) :: values(:,:) ! (layers, sites)
    real(real_8), optional, intent(in) :: wall_values(0:,:)
    real(real_8) :: donors(size(field,3),plan%donor_count), reference(size(field,3)), sample(size(field,3))
    integer :: i,k

    call gather_donors(plan%donors,field,donors)
    do i = 1,size(plan%rows)
      associate(row => plan%rows(i))
        reference = sample_at(row%sample(1))
        values(:,i) = row%constant*reference
        do k = 2,size(row%sample)
          sample = sample_at(row%sample(k))
          values(:,i) = values(:,i)+row%weight(k)*(sample-reference)
        end do
      end associate
    end do

  contains
    function sample_at(index) result(sample)
      integer, intent(in) :: index
      real(real_8) :: sample(size(field,3))
      if (index > 0) then
        sample = donors(:,index)
      else
        sample = wall_values(-index-1,:)
      end if
    end function sample_at
  end subroutine evaluate_reconstruction

  ! The segment may meet a wall at either endpoint. Interior intersections
  ! and overlaps are rejected. Search only wall rays in its bounding box.
  pure logical function visible_donor(geom, start, finish) result(visible)
    type(distributed_geometry), intent(in) :: geom
    integer, intent(in) :: start(2), finish(2)
    integer :: i,j,ilo,ihi,jlo,jhi,mask,r,a(2),b(2),ni

    visible = .true.
    if (all(start == finish)) return
    ni = size(geom%owner,1)
    ilo = max(1,ni-max(start(2),finish(2))/4-1)
    ihi = min(ni,ni-min(start(2),finish(2))/4+1)
    jlo = max(1,min(start(1),finish(1))/4)
    jhi = min(size(geom%owner,2),max(start(1),finish(1))/4+2)
    do i = ilo,ihi
      do j = jlo,jhi
        if (geom%global%kind(i,j) /= wall_node) cycle
        mask = geom%global%wedge(i,j)
        do r = 0,7
          if (btest(mask,r) .and. .not. btest(mask,modulo(r-1,8))) exit
        end do
        a = 2*half_coordinates(i,j,ni)
        b = a+4*[dj(r),-di(r)]
        if (blocks(a,b,start,finish)) then
          visible = .false.
          return
        end if
      end do
    end do
  end function visible_donor

  pure logical function blocks(a,b,p,q)
    integer, intent(in) :: a(2),b(2),p(2),q(2)
    integer(wide) :: c,d,e,f
    integer :: axis

    c = cross(q-p,a-p)
    d = cross(q-p,b-p)
    e = cross(b-a,p-a)
    f = cross(b-a,q-a)
    if (c == 0 .and. d == 0) then
      axis = maxloc(abs(q-p),dim=1)
      blocks = max(min(p(axis),q(axis)),min(a(axis),b(axis))) < min(max(p(axis),q(axis)),max(a(axis),b(axis)))
    else
      blocks = opposite(c,d) .and. opposite(e,f)
      if (c == 0) blocks = blocks .or. strictly_between(a,p,q)
      if (d == 0) blocks = blocks .or. strictly_between(b,p,q)
    end if
  end function blocks

  pure logical function strictly_between(a,p,q)
    integer, intent(in) :: a(2),p(2),q(2)
    integer :: axis

    axis = maxloc(abs(q-p),dim=1)
    strictly_between = min(p(axis),q(axis)) < a(axis) .and. a(axis) < max(p(axis),q(axis))
  end function strictly_between

  pure logical function opposite(a,b)
    integer(wide), intent(in) :: a,b
    opposite = (a < 0 .and. b > 0) .or. (a > 0 .and. b < 0)
  end function opposite

  pure integer(wide) function cross(a,b)
    integer, intent(in) :: a(2),b(2)
    cross = int(a(1),wide)*b(2)-int(a(2),wide)*b(1)
  end function cross

end module qg_reconstruction
