!!$*********************************************************************************************************
!!$  Differential equation to solve:                                                                       *
!!$  f'''(eta) + f(eta)f''(eta) = 0                                                                        *
!!$  where eta is the characteristic boundary layer variable.                                              * 
!!$                                                                                                        *
!!$  Boundary Conditions                                                                                   *
!!$  f(eta=0) = 0; f'(eta=0) = 0                                                                           *
!!$  f'(eta=infinity) = 1                                                                                  *
!!$                                                                                                        *
!!$  To Solve:                                                                                             *
!!$  G = f'                                                                                                *
!!$  H = G' = f''                                                                                          *
!!$  H' = -(1/2)f*H                                                                                        *
!!$                                                                                                        *
!!$  with the same boundary conditions as given above.                                                     *
!!$                                                                                                        *
!!$  4th order Runge Kutta with the shooting technique implemented to obtain solution.                     *
!!$  Shooting technique is used to convert the Boundary Value Problem into an Initial Value Problem.       *
!!$                                                                                                        *
!!$  This is the blasius boundary layer equation.                                                          *
!!$*********************************************************************************************************

module blasius

  implicit none

contains

  subroutine blasius_f_curves(n, f, G, H, H_prime, eta, delta, epsilon)

    implicit none

    integer :: n
    double precision, dimension(0:n) :: f, G, H, H_prime, eta
    double precision :: H1, H2, Gn, Htemp, f_temp, G_temp, H_temp, H_prime_temp
    double precision :: k1,k2,k3,k4
    double precision :: l1,l2,l3,l4
    double precision :: m1,m2,m3,m4
    double precision :: delta, epsilon
    double precision, parameter :: Ginfinity = 2.0d0
    integer :: i, iter_count

    H1 = 1.0d0
    H2 = 1.50d0

    iter_count = 1

    do i = 0, n
       eta(i) = i*delta
    end do

    do while ( abs(G(n) - Ginfinity) .gt. epsilon )

       f(0) = 0.0d0
       G(0) = 0.0d0

       if (iter_count .eq. 1) then
          H(0) = H1
       else if (iter_count .eq. 2) then
          H(0) = H2
       end if

       !G(n) = 1.0d0
       H_prime = -f(1)*H(1)

       do i = 0, n-1

          k1 = delta*G(i)
          l1 = delta*H(i)
          m1 = delta*H_prime(i)

          f_temp = f(i) + 0.50d0*k1
          G_temp = G(i) + 0.50d0*l1
          H_temp = H(i) + 0.50d0*m1
          H_prime_temp = -f_temp*H_temp

          k2 = delta*G_temp
          l2 = delta*H_temp
          m2 = delta*H_prime_temp

          f_temp = f(i) + 0.50d0*k2
          G_temp = G(i) + 0.50d0*l2
          H_temp = H(i) + 0.50d0*m2
          H_prime_temp = -f_temp*H_temp

          k3 = delta*G_temp
          l3 = delta*H_temp
          m3 = delta*H_prime_temp

          f_temp = f(i) + 0.50d0*k3
          G_temp = G(i) + 0.50d0*l3
          H_temp = H(i) + 0.50d0*m3
          H_prime_temp = -f_temp*H_temp

          k4 = delta*G_temp
          l4 = delta*H_temp
          m4 = delta*H_prime_temp

          f(i+1) = f(i) + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
          G(i+1) = G(i) + (1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0d0*l3 + l4)
          H(i+1) = H(i) + (1.0d0/6.0d0)*(m1 + 2.0d0*m2 + 2.0d0*m3 + m4)
          H_prime(i+1) = -f(i+1)*H(i+1)

       end do

       if ( iter_count .eq. 1 ) then
          Gn = G(n)
       end if

       if ( iter_count .eq. 2 ) then

          Htemp = H1 + ((H2 - H1)/(G(n) - Gn))*(Ginfinity-Gn)
          if ( abs(G(n) - Ginfinity) .lt. abs(Gn - Ginfinity) ) then
             H1 = H2
          end if
          H2 = Htemp

       end if

       if ( iter_count .eq. 1) then
          iter_count = iter_count + 1
       else if ( iter_count .eq. 2  ) then
          iter_count = 1
       end if

    end do

  end subroutine blasius_f_curves

  subroutine trapezoidal_integration(eta, func, lower_limit, upper_limit, delta, n, integral)

    implicit none

    integer :: n, i, lower_i, upper_i
    double precision :: lower_limit, upper_limit, lower_blasius_f, upper_blasius_f, delta, integral
    double precision, dimension(0:n) :: eta, func

    do i = 0, n

       if ( eta(i) .eq. lower_limit ) then
          lower_i = i
       else if ( eta(i) .eq. upper_limit ) then
          upper_i = i
       end if

    end do

    integral = 0.0d0

    do i = lower_i, upper_i

       if ( (i .eq. lower_i) .or. (i .eq. upper_i) ) then
          integral = integral + func(i)
       else
          integral = integral + 2.0d0*func(i)
       end if

    end do

    integral = integral*delta*0.50d0

  end subroutine trapezoidal_integration

  subroutine calc_r_function(eta, H, r_function, n, delta, Pr)

    implicit none

    integer :: n, i, j, k
    double precision, dimension(0:n) :: r_function, func_to_integrate, eta, H
    double precision :: integral, delta, temp_var1, temp_var2, temp_var3, temp_var4, temp_var5, Pr

    do i = 0, n

       func_to_integrate = H**2
       call trapezoidal_integration(eta, func_to_integrate, eta(i), eta(n), delta, n, integral)
       temp_var1 = integral

       func_to_integrate = H**(Pr)
       call trapezoidal_integration(eta, func_to_integrate, eta(i), eta(n), delta, n, integral)
       temp_var2 = integral

       temp_var4 = temp_var1 + temp_var2*H(0)**(2.0d0-Pr)

       temp_var1 = 0.0d0
       temp_var2 = 0.0d0
       temp_var3 = 0.0d0

       do j = 1, i
          temp_var1 = temp_var1 + H(j)**(2.0d0-Pr)
       end do
       temp_var1 = temp_var1*H(i)**Pr

       do j = 1, n
          temp_var2 = temp_var2 + H(j)**(2.0d0-Pr)
       end do
       temp_var2 = temp_var2*H(n)**Pr

       do k = i, n
          do j = 1, k
             temp_var3 = temp_var3 + 2.0d0*(H(k)**Pr)*H(j)**(2.0d0-Pr)
          end do
       end do

       temp_var5 = 2.0d0*(temp_var1 + temp_var2 + temp_var3)

       r_function(i) = (Pr*delta**2/8.0d0)*(temp_var4 + temp_var5)
       !write(*,*) "iteration number", i

    end do

  end subroutine calc_r_function

  subroutine calc_y_functions(n, eta, y, Pr, delta)

    implicit none

    integer :: n, i, status_int
    double precision :: y_far_field, yp_far_field, ydp_far_field, Pr, integral, delta
    double precision, dimension(0:n) :: eta, y, yp, ydp
    double precision :: k1,k2,k3,k4
    double precision :: l1,l2,l3,l4
    double precision :: y_temp, yp_temp, ydp_temp
    double precision, parameter :: eta_0 = 4.10d0, eta_n = 10.0d0, eta_size = 1000, delta_ = 0.010d0
    integer, parameter :: func_size = int((eta_n-eta_0)/delta_), eta_start = nint(eta_0/delta_)
    double precision, dimension(1:func_size) :: exp_func, func_to_integrate, x

    !write(*,*) func_size, eta_start

    do i = 1, func_size
       x(i) = dsqrt(Pr)*(eta(eta_start-1+i) - 0.860380d0)
       exp_func(i) = exp(-x(i)**2)
    end do

    func_to_integrate = exp_func
    call trapezoidal_integration(x, func_to_integrate, x(1), x(func_size), delta, func_size, integral)
    y_far_field = integral

    yp_far_field = -exp(-x(1)**2)

    ydp_far_field = -2.0d0*x(1)*exp(-1.0d0*x(1)**2)

    !write(*,*) "eta value", eta_0, y_far_field, yp_far_field, ydp_far_field

    ! Copying y, yp and ydp values from far field to array
    do i = func_size, eta_size
       y(i) = y_far_field
       yp(i) = yp_far_field
       ydp(i) = yp_far_field
    end do

    write(*,*) y(n)

    !call exit(status_int)

    do i = func_size, 1, -1

       k1 = -delta*yp(i)
       l1 = -delta*ydp(i)

       y_temp = y(i) + 0.50d0*k1
       yp_temp = yp(i) + 0.50d0*l1
       ydp_temp = -2.0d0*dsqrt(Pr)*(eta(i)+0.50d0*delta - 0.860380d0)*yp_temp

       k2 = -delta*yp_temp
       l2 = -delta*ydp_temp

       y_temp = y(i) + 0.50d0*k2
       yp_temp = yp(i) + 0.50d0*l2
       ydp_temp = -2.0d0*dsqrt(Pr)*(eta(i)+0.50d0*delta - 0.860380d0)*yp_temp

       k3 = -delta*yp_temp
       l3 = -delta*ydp_temp

       y_temp = y(i) + 0.50d0*k3
       yp_temp = yp(i) + 0.50d0*l3
       ydp_temp = -2.0d0*dsqrt(Pr)*(eta(i)+0.50d0*delta - 0.860380d0)*yp_temp

       k4 = -delta*yp_temp
       l4 = -delta*ydp_temp

       y(i-1) = y(i) + (1.0d0/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
       yp(i-1) = yp(i) + (1.0d0/6.0d0)*(l1 + 2.0d0*l2 + 2.0d0*l3 + l4)
       ydp(i-1) = -2.0d0*dsqrt(Pr)*(eta(i-1) - 0.860380d0)*yp(i-1)

    end do

    ! Scaling y by its initial value such that y(0) = 1
    y = y/y(0)

  end subroutine calc_y_functions

  subroutine calc_y_functions_alternate(n, eta, y, Pr, delta, fdp)

    integer :: n, i, status_int
    double precision :: Pr, integral, delta, temp_var1, temp_var2
    double precision, dimension(0:n) :: eta, y, func_to_integrate, fdp

    func_to_integrate = fdp**Pr
    call trapezoidal_integration(eta, func_to_integrate, eta(0), eta(n), delta, n, integral)
    temp_var2 = integral

    do i = 0, n

       func_to_integrate = fdp**Pr
       call trapezoidal_integration(eta, func_to_integrate, eta(i), eta(n), delta, n, integral)
       temp_var1 = integral

       y(i) = temp_var1/temp_var2

    end do

  end subroutine calc_y_functions_alternate

end module blasius

program compressible_laminar_boundary_layer

  use blasius

  implicit none

  integer, parameter :: n = 100, nx = 1000, ny = 2*n, nz = 4
  double precision, parameter ::  Lx = 4.0d0, Ly = 2.0d0, Lz = 2.0d0
  integer :: i, j, k, status_int
  double precision, parameter :: dx = Lx/nx, dy = Ly/ny, dz = Lz/nz
  double precision, dimension(0:nx, 0:ny, 0:nz) :: u, T, eta_deriv_y, rho, psi_deriv_y
  double precision, dimension(0:nx) :: x
  double precision, dimension(0:ny) :: y
  double precision, dimension(0:nz) :: z
  double precision, dimension(0:n) :: f, G, H, H_prime, eta, r_function, func_to_integrate,test_function, y0, eta_rhs
  double precision, dimension(0:n) :: y0_average, r_function_average, y0_average_deriv, r_function_average_deriv
  double precision, parameter :: delta = 0.10d0, epsilon = 0.000000010d0
  double precision, parameter :: Pr = 0.720d0, mach_no = 2.0d0, gamma = 1.40d0, t0 = 300.0d0, R = 287.150d0, tw = 1000.0d0
  double precision, parameter :: c0 = dsqrt(gamma*R*t0), u_inf = c0*mach_no, p_inf = 101325.0d0, rho_inf = p_inf/(R*t0)
  double precision, parameter :: rho_wall = p_inf/(R*tw)
  double precision, parameter :: dynamic_visc = 0.000019830d0, kinematic_visc = dynamic_visc/rho_inf
  double precision, parameter :: sutherland_constant = 375.3720d0
  double precision, parameter :: C_const = dsqrt(tw/t0)*(t0 + sutherland_constant)/(tw + sutherland_constant)
  double precision, parameter :: transform_const = 0.50d0*dsqrt(u_inf/(kinematic_visc*C_const)), blayer_thickness = 2.0d0*dsqrt(kinematic_visc*C_const/u_inf)
  double precision :: integral, te, a0, temp, y_limit
  integer, dimension(0:nx) :: tempintarray

  status_int = 0

  write(*,*) "boundary layer thickness", blayer_thickness

!!$  do i = 0, n
!!$     test_function(i) = (i*delta)**2
!!$  end do
!!$  call trapezoidal_integration(eta, test_function, eta(0), eta(n), delta, n, integral)
!!$  write(*,*) (n*delta)**3/3.0d0 - integral

!!$  call exit(status_int)

  call blasius_f_curves(n, f, G, H, H_prime, eta, delta, epsilon)

  call calc_r_function(eta, H, r_function, n, delta, Pr)

  te = 1.0d0 + r_function(0)*0.50d0*(gamma-1.0d0)*mach_no**2
  a0 = tw/t0 - te

  call calc_y_functions_alternate(n, eta, y0, Pr, delta, H)

  do i = 0, n

     func_to_integrate = r_function
     call trapezoidal_integration(eta, func_to_integrate, eta(0), eta(i), delta, n, integral)
     r_function_average(i) = integral

     func_to_integrate = y0
     call trapezoidal_integration(eta, func_to_integrate, eta(0), eta(i), delta, n, integral)
     y0_average(i) = integral

     eta_rhs(i) = eta(i) + (gamma-1.0d0)*0.50d0*mach_no**2*r_function_average(i) + a0*y0_average(i)

  end do

  r_function_average(0) = 0.0d0
  y0_average(0) = 0.0d0

  do i = 0, n

     if ( i .eq. 0 ) then
        r_function_average_deriv(i) = (r_function_average(i+1) - r_function_average(i))/(delta)
        y0_average_deriv(i) = (y0_average(i+1) - y0_average(i))/(delta)
     else if ( i .eq. n ) then
        r_function_average_deriv(i) = (r_function_average(i) - r_function_average(i-1))/delta
        y0_average_deriv(i) = (y0_average(i) - y0_average(i-1))/delta
     else if ( i .ge. 2 .and. i .le. n-2 ) then
        r_function_average_deriv(i) = 8.0d0*(r_function_average(i+1) - r_function_average(i-1) - r_function_average(i+2) - r_function_average(i-2))/(12.0d0*delta)
        y0_average_deriv(i) = 8.0d0*(y0_average(i+1) - y0_average(i-1) - y0_average(i+2) - y0_average(i-2))/(12.0d0*delta)
     else
        r_function_average_deriv(i) = (r_function_average(i+1) - r_function_average(i-1))/(2.0d0*delta)
        y0_average_deriv(i) = (y0_average(i+1) - y0_average(i-1))/(2.0d0*delta)
     end if

  end do

  open(unit = 50,file="eta_rhs.dat")
  write(50,*)'variables="eta", "eta_rhs"'
  do i = 0, n
     write(50,'(2E30.18)') i*delta, eta_rhs(i)
  end do
  close(50)

  open(unit=45,file="blasius_functions.dat")
  write(45,*) 'variables = "i","f","fp","fdp","r_eta","r_eta_avg","r_eta_avg_prime","y0","y0_avg","y0_avg_prime"'
  do i = 0, n
     eta(i) = i*delta
     write(45, '(10E30.18)') i*delta, f(i), G(i), H(i), r_function(i), r_function_average(i), r_function_average_deriv(i), y0(i), y0_average(i), y0_average_deriv(i)
  end do
  close(45)

!!$ On obtaining the boundary layer x and y solutions from the eta coordinates.

  do k = 0, nz
     do j = 0, ny
        do i = 0, nx
           x(i) = i*dx
           y(j) = j*dy
           z(k) = k*dz
!!$           y(i,j,k) = eta_rhs(j)*dsqrt(x(i))/transform_const
        end do
     end do
  end do

!!$ Background velocity and temperature filled throughout the domain
  u = u_inf
  T = t0   

  open(unit = 45, file = "boundary_layer_profile.dat")
  write(45,*)'variables="x","delta"'
  do i = 0, nx
     write(45,'(2E30.18)') x(i), blayer_thickness*dsqrt(x(i))
  end do
  close(45)

  call exit(status_int)

  do k = 0, nz
     do i = 1, nx
        j = 0
        do while ( y(j)/(blayer_thickness*dsqrt(x(i))  .ge. 1.0d0)
           
           T(i,j,k) = 1.0d0 + (gamma-1.0d0)*0.50d0*mach_no**2*r_function(j) + a0*y0(j)
           u(i,j,k) = 0.50d0*G(j)

        

           if ( y(j) .eq. 0.0d0 ) then
!!$              write(*,'(5E30.18)') u(i,j,k), 1.0d0*i, 1.0d0*j, 1.0d0*k, psi_deriv_y(i,j,k)
              rho(i,j,k) = rho_wall
           else
              temp = 1.0d0 + (gamma-1.0d0)*0.50d0*mach_no**2*r_function_average_deriv(j) + a0*y0_average_deriv(j)
              eta_deriv_y(i,j,k) = 0.50d0*dsqrt(u_inf/(kinematic_visc*x(i)*C_const))/temp
              psi_deriv_y(i,j,k) = dsqrt(kinematic_visc*u_inf*C_const)*dsqrt(x(i))*G(j)*eta_deriv_y(i,j,k)
              rho(i,j,k) = (rho_inf/u_inf*u(i,j,k))*psi_deriv_y(i,j,k)
           end if

        end do
     end do
  end do

!!$  write(*,*) maxval(rho), minval(rho)

!!$  call exit(status_int)

!!$ Writing profile to file	
  open(unit = 40, file="profile.dat")
  write(40,*) 'variables="x","y","z","u","T"'
  write(40,*) 'zone I=',nx+1,',J=',ny+1,',K=',nz+1,',F=POINT'
  do k = 0, nz
     do j = 0, ny
        do i = 0, nx
           write(40,'(5E30.18)') x(i), y(j), z(k), u(i,j,k)*u_inf, T(i,j,k)*t0
        end do
     end do
  end do
  close(40)

end program compressible_laminar_boundary_layer
