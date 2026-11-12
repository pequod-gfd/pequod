#include "config.h"

module qg_irregular_diagnostics
  use qg_diagnostics, only: stratification
  use qg_geometry, only: wall_node
  use qg_irregular_parallel, only: distributed_geometry
  use qg_parallel, only: pabort
  use qg_parameters, only: real_8
  implicit none
  include 'mpif.h'
  private
  public :: irregular_energy_totals, set_wall_pv

contains

  ! Layer kinetic and potential energy totals weighted by actual fluid area.
  ! Derivatives use native grid units, as in the rectangle diagnostic.
  ! psi has one populated halo ring wherever the neighboring node is wet.
  subroutine irregular_energy_totals(geom, strat, psi, area, kinetic, potential)
    type(distributed_geometry), intent(in) :: geom
    type(stratification), intent(in) :: strat
    real(real_8), intent(in) :: psi(:,:,:)
    real(real_8), intent(out) :: area, kinetic(strat%layers), potential
    integer, parameter :: ci(4) = [1,1,0,0], cj(4) = [0,1,1,0] ! SW, SE, NE, NW
    integer :: i,j,k,n,gi,gj,ierr,c,wet,next,prev,wedge_count(4),ni,nj
    real(real_8) :: dx,dy,value(4),depth(strat%layers),weight,local(strat%layers+2),total(strat%layers+2)

    n = strat%layers
    if (any(shape(psi) /= [geom%ni+2,geom%nj+2,n])) then
      labort('energy diagnostic requires owned streamfunction and one halo ring')
    end if
    ni = size(geom%global%wedge,1)
    nj = size(geom%global%wedge,2)
    depth = strat%h/sum(strat%h)
    local = 0.0D0
    do j = 1,geom%nj
      gj = geom%j0+j
      do i = 1,geom%ni
        gi = geom%i0+i
        ! Lumped dual-cell potential energy uses the same quadrature as mass.
        if (geom%area(i,j) > 0.0D0) then
          local(n+2) = local(n+2)+geom%area(i,j)
          do k = 1,n-1
            weight = 0.5D0*geom%area(i,j)*depth(k)*strat%s2(k)
            local(n+1) = local(n+1)+weight*(psi(i+1,j+1,k)-psi(i+1,j+1,k+1))**2
          end do
        end if

        if (gi == ni .or. gj == nj) cycle
        do c = 1,4
          wedge_count(c) = popcnt(ibits(geom%global%wedge(gi+ci(c),gj+cj(c)),2*(c-1),2))
        end do
        select case(sum(wedge_count))
          case(0)
            cycle
          case(8)
            ! Full primal square: original cell-centre gradient.
            do k = 1,n
              do c = 1,4
                value(c) = psi(i+ci(c)+1,j+cj(c)+1,k)
              end do
              dx = 0.5D0*(value(2)+value(3)-value(1)-value(4))
              dy = 0.5D0*(value(3)+value(4)-value(1)-value(2))
              local(k) = local(k)+0.5D0*geom%h**2*depth(k)*(dx*dx+dy*dy)
            end do
          case(4)
            ! Wet half square: a P1 triangle, with exact integrated gradient.
            wet = maxloc(wedge_count,dim=1)
            next = mod(wet,4)+1
            prev = modulo(wet-2,4)+1
            do k = 1,n
              dx = psi(i+ci(next)+1,j+cj(next)+1,k)-psi(i+ci(wet)+1,j+cj(wet)+1,k)
              dy = psi(i+ci(prev)+1,j+cj(prev)+1,k)-psi(i+ci(wet)+1,j+cj(wet)+1,k)
              local(k) = local(k)+0.25D0*geom%h**2*depth(k)*(dx*dx+dy*dy)
            end do
          case default
            labort('invalid wet primal square in energy diagnostic')
        end select
      end do
    end do
    call mpi_allreduce(local,total,size(local),MPI_DOUBLE_PRECISION,MPI_SUM,geom%comm,ierr)
    kinetic = total(:n)
    potential = total(n+1)
    area = total(n+2)
  end subroutine irregular_energy_totals

  ! Relative wall PV is omega + S psi. It remains a wall point diagnostic.
  subroutine set_wall_pv(geom, strat, wall_values, omega, zeta)
    type(distributed_geometry), intent(in) :: geom
    type(stratification), intent(in) :: strat
    real(real_8), intent(in) :: wall_values(0:,:), omega(:,:,:)
    real(real_8), intent(inout) :: zeta(:,:,:)
    real(real_8) :: stretching(strat%layers,0:ubound(wall_values,1))
    integer :: i,j,component

    do component = 0,ubound(wall_values,1)
      stretching(:,component) = matmul(strat%s_mat,wall_values(component,:))
    end do
    do j = 1,geom%nj
      do i = 1,geom%ni
        if (geom%global%kind(geom%i0+i,geom%j0+j) /= wall_node) cycle
        component = geom%global%component(geom%i0+i,geom%j0+j)
        zeta(i,j,:) = omega(i,j,:)+stretching(:,component)
      end do
    end do
  end subroutine set_wall_pv

end module qg_irregular_diagnostics
