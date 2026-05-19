!- 1次元線形浅水方程式+コリオリ：Runge Kutta法
program main
  implicit none

  integer(8),parameter :: m = 1000

  real(8),parameter :: h_m = 1.d2       !- depth [m]
  real(8),parameter :: l_m = 1.d7        !- length [m]
  real(8),parameter :: g_mps2 = 10.d0    !- 重力加速度 [m/s2]
  real(8),parameter :: f_psec = 1.d-4    !- コリオリ・パラメータ [/s]
  real(8),parameter :: eta0_m = 1.d0   !- 初期振幅 [m]
  real(8),parameter :: pi = 3.1415926535897932384d0
  real(8),parameter :: tend_sec = 2.d5   !- 計算期間 [sec]
  real(8),parameter :: dt_sec = 2.d2     !- t方向の刻み幅 [sec]
  real(8)    :: x(0:m,0:2)               !- 変数 x(:,0) eta [m], x(:,1) u [m/s]
  real(8)    :: k1(0:m,0:2), k2(0:m,0:2), k3(0:m,0:2), k4(0:m,0:2)
  real(8)    :: t, dx_m, x_m
  integer(8) :: ii, n, nm

  nm = int( tend_sec / dt_sec )          !- 総ステップ数
  dx_m = l_m / dble(m)

  !- 初期値
  do ii = 0, m-1
    x_m = ( dble(ii) + 0.5d0 ) * dx_m
!    x(ii,0) = eta0_m * cos( pi * x_m / l_m )
    x(ii,0) = eta0_m * tanh(  ( x_m - l_m * 0.5d0 ) / ( 10.d0 * dx_m ) )
  enddo
  x(m,0) = 0.d0  !- 壁の中(使わない)
  x(:,1) = 0.d0
  x(:,2) = 0.d0
  t    = 0.d0

  open(1,file='shallow.dat',form='unformatted',access='stream')
  write(1) m, nm
  write(1) dx_m, dt_sec
  write(1) x
  print *,'                   t,                  eta(middle),                 u(middle)'
  do n = 0, nm-1
     k1 = dt_sec*f(x,t)
     k2 = dt_sec*f( x+0.5d0*k1, t + 0.5d0*dt_sec)
     k3 = dt_sec*f( x+0.5d0*k2, t + 0.5d0*dt_sec)
     k4 = dt_sec*f( x+k3,       t + dt_sec)
     x = x + ( k1 +2.d0*k2 +2.d0*k3 +k4 ) / 6.d0
     t = t + dt_sec
     if ( ( mod(n,(nm/10)) == 0 ).or. ( n == nm-1 ) ) print*,t,x(m/2,0:1)
     write(1) x
  enddo
  close(1)

contains

  function f(x,t)
    real(8) :: f(0:m,0:2)
    real(8) :: x(0:m,0:2), t
    integer(8) :: i

    f(m,0) = 0.d0  !- 壁の中
    f(0,1:2) = 0.d0  !- u 境界条件
    f(m,1:2) = 0.d0

    do i = 0, m-1
      f(i,0) =  - h_m / dx_m * ( x(i+1,1) - x(i,1) )
    enddo
    do i = 1, m-1
      f(i,1) =  f_psec * x(i,2) - g_mps2 / dx_m * ( x(i,0) - x(i-1,0) )
    enddo
    do i = 1, m-1
      f(i,2) =  - f_psec * x(i,1)
    enddo
  end function f

end program main
