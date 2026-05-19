!- 2次元線形浅水方程式：Runge Kutta法
program main
  implicit none

  integer(8),parameter :: im = 200, jm = 200

  real(8),parameter :: h_m = 1.d3       !- depth [m]
  real(8),parameter :: lx_m = 1.d7        !- length [m]
  real(8),parameter :: ly_m = 1.d7        !- length [m]
  real(8),parameter :: g_mps2 = 10.d0    !- 重力加速度 [m/s2]
  real(8),parameter :: rho_kgpm3 = 1.d3  !- 密度 [kg/m3]
  real(8),parameter :: taux_amp_npm2 = 0.d-1    !- 風応力振幅 [N/m2]
  real(8),parameter :: pi = 3.1415926535897932384d0
  real(8),parameter :: tend_sec = 2.d5   !- 計算期間 [sec]
  real(8),parameter :: dt_sec = 5.d2     !- t方向の刻み幅 [sec]
  real(8),parameter :: eta0_m = 1.d1     !- 以下、初期値。振幅 [m]
  real(8),parameter :: x0_m = lx_m * 0.5d0, y0_m = ly_m * 0.5d0    !- 山の中心 [m]
  real(8),parameter :: sigma_m = 5.d5    !- 山の水平大きさ [m]
  real(8)    :: x(0:im,0:jm,0:2)               !- 変数 x(:,0) eta [m], x(:,1) u [m/s], x(:,2) v [m/s]
  real(8)    :: k1(0:im,0:jm,0:2), k2(0:im,0:jm,0:2), k3(0:im,0:jm,0:2), k4(0:im,0:jm,0:2)
  real(8)    :: taux_npm2(0:im,0:jm)       !- 風応力 [N/m^2]
!  real(8),parameter :: f_psec = 1.d-4    !- コリオリ・パラメータ [/s]
  real(8)    :: f_psec(0:jm)              !- コリオリ・パラメータ [/s]
  real(8)    :: t, dx_m, x_m, dy_m, y_m
  integer(8) :: i, j, n, nm

  nm = int( tend_sec / dt_sec )          !- 総ステップ数
  dx_m = lx_m / dble(im)
  dy_m = ly_m / dble(jm)

  !- コリオリ・パラメータ
  do j = 0, jm
    y_m = dble(j) * dy_m
    f_psec(j) =  ( y_m - 0.5d0*ly_m ) * 0.d-11 + 1.d-4
  enddo

  !- 風応力
  do j = 0, jm
    do i = 0, im
      x_m = dble(i) * dx_m
      y_m = dble(j) * dy_m
      taux_npm2(i,j) = - taux_amp_npm2 * cos( pi * y_m / ly_m )
    enddo
  enddo

  !- 初期値
  do j = 0, jm-1
    do i = 0, im-1
      x_m = ( dble(i) + 0.5d0 ) * dx_m
      y_m = ( dble(j) + 0.5d0 ) * dy_m
      x(i,j,0) = eta0_m * exp( - ( (x_m-x0_m)**2 + (y_m-y0_m)**2 ) / ( 2.d0 * sigma_m**2 ) )
    enddo
  enddo
  x(im,:,0) = 0.d0  !- 壁の中(使わない)
  x(:,jm,0) = 0.d0
  x(:,:,1:2) = 0.d0
  t    = 0.d0

  !- 出力の用意
  open(1,file='shallow2d.dat',form='unformatted',access='stream')
  write(1) im, jm, nm
  write(1) dx_m, dy_m, dt_sec
  write(1) x
  print *,'                   t,                  eta(center),                 u,                 v'

  !- メインループ
  do n = 0, nm-1
     k1 = dt_sec*f(x,t)
     k2 = dt_sec*f( x+0.5d0*k1, t + 0.5d0*dt_sec)
     k3 = dt_sec*f( x+0.5d0*k2, t + 0.5d0*dt_sec)
     k4 = dt_sec*f( x+k3,       t + dt_sec)
     x = x + ( k1 +2.d0*k2 +2.d0*k3 +k4 ) / 6.d0
     x(im,:,0) = 0.d0  !- 壁の中
     x(:,jm,0) = 0.d0
     x(0,:,1:2) = 0.d0  !- u,v 境界条件
     x(im,:,1:2) = 0.d0
     x(:,0,1:2) = 0.d0
     x(:,jm,1:2) = 0.d0
     t = t + dt_sec
     if ( ( mod(n,(nm/10)) == 0 ).or. ( n == nm-1 ) ) print*,t,x(im/2,jm/2,0:2)
     write(1) x
  enddo
  close(1)

contains

  !- 時間変化率の計算
  function f(x,t)
    real(8) :: f(0:im,0:jm,0:2)
    real(8) :: x(0:im,0:jm,0:2), t
    integer(8) :: i, j

    f(:,:,:) = 0.d0  !- 予測しない格子では時間変化ゼロを返しておく
    do j = 0, jm-1
      do i = 0, im-1
        f(i,j,0) =  - 0.5d0 * h_m *( ( x(i+1,j,1) - x(i,j,1) + x(i+1,j+1,1) - x(i,j+1,1) ) / dx_m &
             &                      +( x(i,j+1,2) - x(i,j,2) + x(i+1,j+1,2) - x(i+1,j,2) ) / dy_m )
      enddo
    enddo
    do j = 1, jm-1
      do i = 1, im-1
        f(i,j,1) = f_psec(j) * x(i,j,2) - 0.5d0 * g_mps2 / dx_m &
                  & * ( x(i,j,0) - x(i-1,j,0) + x(i,j-1,0) - x(i-1,j-1,0) ) &
                  & + taux_npm2(i,j) / rho_kgpm3 / h_m
        f(i,j,2) = - f_psec(j) * x(i,j,1) - 0.5d0 * g_mps2 / dy_m &
                  & * ( x(i,j,0) - x(i,j-1,0) + x(i-1,j,0) - x(i-1,j-1,0) )
      enddo
    enddo

  end function f

end program main
