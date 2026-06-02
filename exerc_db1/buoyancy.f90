!- 浮力振動：オイラー法
program main
  implicit none

  real(8),parameter :: dn2 = 1.d-4       !- ブラントバイサラ振動数の2乗
  real(8),parameter :: tend_sec = 1.d3   !- 計算期間
  real(8),parameter :: dt_sec = 10.d0    !- t方向の刻み幅
  real(8)    :: x(0:1)                   !- 変数 x(0): z [m], x(1): w (= dz/dt) [m/s]
  real(8)    :: t_sec
  integer(8) :: i, nm

  nm = int( tend_sec / dt_sec )                   !- 総ステップ数

  x(0)  = 1.d0  !- 初期値
  x(1)  = 0.d0
  t_sec = 0.d0

  open(10,file='buoyancy.dat',form='unformatted',access='stream')
  write(10) nm
  write(10) dt_sec 
  write(10) x

  write(6,*) '                   t,                        z,                        w'
  do i = 0, nm-1
     x = x + dt_sec * f(x,t_sec)    !- x,fとも2要素。オイラー法
     t_sec = t_sec + dt_sec
     if ( ( mod(i,(nm/10)) == 0 ).or.( i == nm-1 ) ) write(6,*) t_sec,x
     write(1) x
  enddo
  close(1)

contains

  function f(x,t_sec)
    real(8) :: f(0:1)
    real(8) :: x(0:1), t_sec   !- t_secは使わないが、常微分方程式と合わせるために入力に含める
    f(0) = x(1)
    f(1) = -dn2 * x(0)
  end function f

end program main
