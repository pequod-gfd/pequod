!!NP 12

#include "config.h"

program test_helmholtz_12

  use qg_laplace_bc
  use qg_parallel
  use qg_parameters, only : pi, real_8
  use qg_helmholtz
  use test_tools
  
  implicit none
  
  integer :: g_ii, g_jj, i, j, n, p
  integer, dimension(6) :: p_iis = (/1, 2, 3, 4, 6, 12/), p_jjs = (/12, 6, 4, 3, 2, 1/)
  integer, dimension(7), parameter :: g_iis = (/33, 65, 129, 257, 513, 1025, 2049/)
  integer, dimension(:), allocatable :: iis, jjs
  real(kind = real_8) :: x, y
  real(kind = real_8), dimension(size(g_iis)) :: m, o
  real(kind = real_8), dimension(:, :, :), allocatable :: phi
  type(decomposition) :: decomp
  type(helm_cache) :: helm_c

  call parallel_init()
  
  n = nprocs()
  call preport_test("nprocs == 12", n /= 12)
  
  do p = 1, size(p_iis)
    allocate(iis(p_iis(p)), jjs(p_jjs(p)))
    
    do n = 1, size(g_iis)
      g_ii = g_iis(n)
      g_jj = g_iis(n)
      
      iis(:) = g_ii / size(iis)
      jjs(:) = g_jj / size(jjs)
      do i = 1, g_ii - iis(1) * size(iis)
        iis(size(iis) - i + 1) = iis(size(iis) - i + 1) + 1
      end do
      assert(sum(iis) == g_ii)
      do i = 1, g_jj - jjs(1) * size(jjs)
        jjs(size(jjs) - i + 1) = jjs(size(jjs) - i + 1) + 1
      end do
      assert(sum(jjs) == g_jj)
      
      call allocate(decomp, iis = iis, jjs = jjs, layers = 1)
      call allocate(helm_c, decomp, ss = (/0.0D0/))
      
      allocate(phi(decomp%ii, decomp%jj, decomp%layers))
      
      phi = huge(0.0D0)
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          phi(i, j, 1) = -pi * pi * dsin(pi * x) * dsin(pi * y)
        end do
      end do
      phi = phi * ((1.0D0 / real(decomp%g_ii_int, kind = real_8) ** 2) + (1.0D0 / real(decomp%g_jj_int, kind = real_8) ** 2))
      
      call helm(helm_c, phi)
      call apply_dirichlet_bc(decomp, 0.0D0, phi(:, :, 1))

      m(n) = 0.0D0
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          m(n) = max(m(n), dabs(phi(i, j, 1) - dsin(pi * x) * dsin(pi * y)))
        end do
      end do
      call allmax(m(n))
      
      deallocate(phi)
      
      call deallocate(helm_c)
      call deallocate(decomp)
    end do
    
    do n = 2, size(g_iis)
      o(n) = -dlog10(m(n) / m(n - 1)) / dlog10(real(g_iis(n), kind = real_8) / real(g_iis(n - 1), kind = real_8))
    end do
    if(procno() == 1) write(6, *) o(2:)
    call preport_test("helm Poisson equation " // int2str(size(iis)) // "x" // int2str(size(jjs)) // " square domain convergence", any(o(2:) < 2.0D0))
    
    do n = 1, size(g_iis)
      g_ii = g_iis(n)
      g_jj = g_iis(n)
      
      iis(:) = g_ii / size(iis)
      jjs(:) = g_jj / size(jjs)
      do i = 1, g_ii - iis(1) * size(iis)
        iis(size(iis) - i + 1) = iis(size(iis) - i + 1) + 1
      end do
      assert(sum(iis) == g_ii)
      do i = 1, g_jj - jjs(1) * size(jjs)
        jjs(size(jjs) - i + 1) = jjs(size(jjs) - i + 1) + 1
      end do
      assert(sum(jjs) == g_jj)
      
      call allocate(decomp, iis = iis, jjs = jjs, layers = 1)
      call allocate(helm_c, decomp, ss = (/0.5D0/))
      
      allocate(phi(decomp%ii, decomp%jj, decomp%layers))
      
      phi = huge(0.0D0)
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          phi(i, j, 1) = - dsin(pi * x) * dsin(pi * y)
        end do
      end do
      phi = phi * ((pi * pi / real(decomp%g_ii_int, kind = real_8) ** 2) + (pi * pi / real(decomp%g_jj_int, kind = real_8) ** 2) + 0.5D0)
      
      call helm(helm_c, phi)
      call apply_dirichlet_bc(decomp, 0.0D0, phi(:, :, 1))

      m(n) = 0.0D0
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          m(n) = max(m(n), dabs(phi(i, j, 1) - dsin(pi * x) * dsin(pi * y)))
        end do
      end do
      call allmax(m(n))
      
      deallocate(phi)
      
      call deallocate(helm_c)
      call deallocate(decomp)
    end do
    
    do n = 2, size(g_iis)
      o(n) = -dlog10(m(n) / m(n - 1)) / dlog10(real(g_iis(n), kind = real_8) / real(g_iis(n - 1), kind = real_8))
    end do
    if(procno() == 1) write(6, *) o(2:)
    call preport_test("helm Helmholtz equation " // int2str(size(iis)) // "x" // int2str(size(jjs)) // " square domain convergence", any(o(2:) < 3.9D0))
    
    do n = 1, size(g_iis)
      g_ii = g_iis(n)
      g_jj = (g_iis(n) - 1) * 2 + 1
      
      iis(:) = g_ii / size(iis)
      jjs(:) = g_jj / size(jjs)
      do i = 1, g_ii - iis(1) * size(iis)
        iis(size(iis) - i + 1) = iis(size(iis) - i + 1) + 1
      end do
      assert(sum(iis) == g_ii)
      do i = 1, g_jj - jjs(1) * size(jjs)
        jjs(size(jjs) - i + 1) = jjs(size(jjs) - i + 1) + 1
      end do
      assert(sum(jjs) == g_jj)
      
      call allocate(decomp, iis = iis, jjs = jjs, layers = 1)
      call allocate(helm_c, decomp, ss = (/0.0D0/))
      
      allocate(phi(decomp%ii, decomp%jj, decomp%layers))
      
      phi = huge(0.0D0)
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          phi(i, j, 1) = -pi * pi * dsin(pi * x) * dsin(pi * y)
        end do
      end do
      phi = phi * ((1.0D0 / real(decomp%g_ii_int, kind = real_8) ** 2) + (1.0D0 / real(decomp%g_jj_int, kind = real_8) ** 2))
      
      call helm(helm_c, phi)
      call apply_dirichlet_bc(decomp, 0.0D0, phi(:, :, 1))

      m(n) = 0.0D0
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          m(n) = max(m(n), dabs(phi(i, j, 1) - dsin(pi * x) * dsin(pi * y)))
        end do
      end do
      call allmax(m(n))
      
      deallocate(phi)
      
      call deallocate(helm_c)
      call deallocate(decomp)
    end do
    
    do n = 2, size(g_iis)
      o(n) = -dlog10(m(n) / m(n - 1)) / dlog10(real(g_iis(n), kind = real_8) / real(g_iis(n - 1), kind = real_8))
    end do
    if(procno() == 1) write(6, *) o(2:)
    call preport_test("helm Poisson equation " // int2str(size(iis)) // "x" // int2str(size(jjs)) // " non-square domain convergence", any(o(2:) < 2.0D0))
    
    do n = 1, size(g_iis)
      g_ii = g_iis(n)
      g_jj = (g_iis(n) - 1) * 2 + 1
      
      iis(:) = g_ii / size(iis)
      jjs(:) = g_jj / size(jjs)
      do i = 1, g_ii - iis(1) * size(iis)
        iis(size(iis) - i + 1) = iis(size(iis) - i + 1) + 1
      end do
      assert(sum(iis) == g_ii)
      do i = 1, g_jj - jjs(1) * size(jjs)
        jjs(size(jjs) - i + 1) = jjs(size(jjs) - i + 1) + 1
      end do
      assert(sum(jjs) == g_jj)
      
      call allocate(decomp, iis = iis, jjs = jjs, layers = 1)
      call allocate(helm_c, decomp, ss = (/0.5D0/))
      
      allocate(phi(decomp%ii, decomp%jj, decomp%layers))
      
      phi = huge(0.0D0)
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          phi(i, j, 1) = - dsin(pi * x) * dsin(pi * y)
        end do
      end do
      phi = phi * ((pi * pi / real(decomp%g_ii_int, kind = real_8) ** 2) + (pi * pi / real(decomp%g_jj_int, kind = real_8) ** 2) + 0.5D0)
      
      call helm(helm_c, phi)
      call apply_dirichlet_bc(decomp, 0.0D0, phi(:, :, 1))

      m(n) = 0.0D0
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          m(n) = max(m(n), dabs(phi(i, j, 1) - dsin(pi * x) * dsin(pi * y)))
        end do
      end do
      call allmax(m(n))
      
      deallocate(phi)
      
      call deallocate(helm_c)
      call deallocate(decomp)
    end do
    
    do n = 2, size(g_iis)
      o(n) = -dlog10(m(n) / m(n - 1)) / dlog10(real(g_iis(n), kind = real_8) / real(g_iis(n - 1), kind = real_8))
    end do
    if(procno() == 1) write(6, *) o(2:)
    call preport_test("helm Helmholtz equation " // int2str(size(iis)) // "x" // int2str(size(jjs)) // " square domain convergence", any(o(2:) < 3.9D0))
    
    do n = 1, size(g_iis)
      g_ii = g_iis(n)
      g_jj = (g_iis(n) - 1) * 2
      
      iis(:) = g_ii / size(iis)
      jjs(:) = g_jj / size(jjs)
      do i = 1, g_ii - iis(1) * size(iis)
        iis(size(iis) - i + 1) = iis(size(iis) - i + 1) + 1
      end do
      assert(sum(iis) == g_ii)
      do i = 1, g_jj - jjs(1) * size(jjs)
        jjs(size(jjs) - i + 1) = jjs(size(jjs) - i + 1) + 1
      end do
      assert(sum(jjs) == g_jj)
      
      call allocate(decomp, iis = iis, jjs = jjs, layers = 1, periodic_j = .true.)
      call allocate(helm_c, decomp, ss = (/0.0D0/))
      
      allocate(phi(decomp%ii, decomp%jj, decomp%layers))
      
      phi = huge(0.0D0)
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          phi(i, j, 1) = -pi * pi * dcos(2.0D0 * pi * x) * dsin(pi * y)
        end do
      end do
      phi = phi * ((1.0D0 / real(decomp%g_ii_int, kind = real_8) ** 2) + (4.0D0 / real(decomp%g_jj_int, kind = real_8) ** 2))
      
      call helm(helm_c, phi)
      call apply_dirichlet_bc(decomp, 0.0D0, phi(:, :, 1))

      m(n) = 0.0D0
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          m(n) = max(m(n), dabs(phi(i, j, 1) - dcos(2.0D0 * pi * x) * dsin(pi * y)))
        end do
      end do
      call allmax(m(n))
      
      deallocate(phi)
      
      call deallocate(helm_c)
      call deallocate(decomp)
    end do
    
    do n = 2, size(g_iis)
      o(n) = -dlog10(m(n) / m(n - 1)) / dlog10(real(g_iis(n), kind = real_8) / real(g_iis(n - 1), kind = real_8))
    end do
    if(procno() == 1) write(6, *) o(2:)
    call preport_test("helm Poisson equation " // int2str(size(iis)) // "x" // int2str(size(jjs)) // " non-square domain convergence, periodic", any(o(2:) < 2.0D0))
    
    do n = 1, size(g_iis)
      g_ii = g_iis(n)
      g_jj = (g_iis(n) - 1) * 2
      
      iis(:) = g_ii / size(iis)
      jjs(:) = g_jj / size(jjs)
      do i = 1, g_ii - iis(1) * size(iis)
        iis(size(iis) - i + 1) = iis(size(iis) - i + 1) + 1
      end do
      assert(sum(iis) == g_ii)
      do i = 1, g_jj - jjs(1) * size(jjs)
        jjs(size(jjs) - i + 1) = jjs(size(jjs) - i + 1) + 1
      end do
      assert(sum(jjs) == g_jj)
      
      call allocate(decomp, iis = iis, jjs = jjs, layers = 1, periodic_j = .true.)
      call allocate(helm_c, decomp, ss = (/0.5D0/))
      
      allocate(phi(decomp%ii, decomp%jj, decomp%layers))
      
      phi = huge(0.0D0)
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          phi(i, j, 1) = - dcos(2.0D0 * pi * x) * dsin(pi * y)
        end do
      end do
      phi = phi * ((pi * pi / real(decomp%g_ii_int, kind = real_8) ** 2) + (4.0D0 * pi * pi / real(decomp%g_jj_int, kind = real_8) ** 2) + 0.5D0)
      
      call helm(helm_c, phi)
      call apply_dirichlet_bc(decomp, 0.0D0, phi(:, :, 1))

      m(n) = 0.0D0
      do i = 1, decomp%ii
        do j = 1, decomp%jj
          x = real(j + decomp%jj_0 - 1, kind = real_8) / real(decomp%g_jj_int, kind = real_8) 
          y = real(decomp%g_ii - i - decomp%ii_0, kind = real_8) / real(decomp%g_ii_int, kind = real_8) 
          m(n) = max(m(n), dabs(phi(i, j, 1) - dcos(2.0D0 * pi * x) * dsin(pi * y)))
        end do
      end do
      call allmax(m(n))
      
      deallocate(phi)
      
      call deallocate(helm_c)
      call deallocate(decomp)
    end do
    
    do n = 2, size(g_iis)
      o(n) = -dlog10(m(n) / m(n - 1)) / dlog10(real(g_iis(n), kind = real_8) / real(g_iis(n - 1), kind = real_8))
    end do
    if(procno() == 1) write(6, *) o(2:)
    call preport_test("helm Helmholtz equation equation " // int2str(size(iis)) // "x" // int2str(size(jjs)) // " non-square domain convergence, periodic", any(o(2:) < 3.9D0))
    
    deallocate(iis, jjs)
  end do
  
  call parallel_finalise()

end program test_helmholtz_12
