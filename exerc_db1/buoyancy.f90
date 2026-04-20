
!- 浮力振動：オイラー法
program main
  implicit none

  real(8),parameter :: dn2 = 1.d-4       !- ブラントバイサラ振動数の2乗
  real(8),parameter :: tend = 1.d3       !- 計算期間
  real(8),parameter :: h = 10.d0         !- t方向の刻み幅
  real(8)    :: x(0:1)                   !- 変数 x(0): z, x(1): w (= dz/dt)
  real(8)    :: t
  integer(8) :: i, nm

  nm = int( tend / h )                   !- 総ステップ数

  x(0)=1
  x(1)=0
  t=0

  open(1,file='buoyancy.dat',form='unformatted',access='stream')
  write(1) nm
  write(1) h 
  write(1) x
  print *,'                   t,                        z,                        w'
  do i = 0, nm-1
     x = x + h * f(x,t)    !- x,fとも2要素。オイラー法
     t = t + h
     if ( ( mod(i,(nm/10)) == 0 ).or. ( i == nm-1 ) ) print*,t,x
     write(1) x
  enddo
  close(1)

contains

  function f(x,t)
    real(8) :: f(0:1)
    real(8) :: x(0:1), t   !- tは使わないが、常微分方程式と合わせるために入力に含める
    f(0) = x(1)
    f(1) = -dn2*x(0)
  end function f

end program main
