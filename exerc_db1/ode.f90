!- 常微分方程式：オイラー法
program ode
  implicit none

  integer(8),parameter :: nm = 10  !- ステップ数
  real(8)    :: t, h, x, xa
  integer(8) :: i

  h = 1.d0 / dble(nm) ! t方向の刻み幅の設定, 1.d0 は倍精度 1. x 10^0 を明示
  x = 0.5d0
  t = 0.d0

  open(1,file='ode.dat',form='unformatted',access='stream')
  write(1) nm
  write(1) h

  print *,'                   t,                       x,                      true,                     error'
  do i = 0, nm-1
     x = x + h*f(x,t) ! オイラー法
     t = t + h
     xa = t**3 + 1.d0 / ( 1.d0 + exp(-t) )  !- 解析解(真値)
     print*,t,x,xa,x-xa
     write(1) x,xa
  enddo
  close(1)

contains          ! 以下の関数副プログラムを内部関数とする
  function f(x,t) ! 関数副プログラムの定義
    real(8)            :: f
    real(8),intent(in) :: t, x
    f = 3.d0*t**2 - t**3 - t**6 + ( 2.d0*t**3 + 1.d0 )*x - x**2
  end function f

end program ode
