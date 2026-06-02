!- 常微分方程式：オイラー法
program ode
  implicit none

  integer(8),parameter :: nm = 10  !- ステップ数
  real(8)    :: t, dt, x, xa
  integer(8) :: i

  dt = 1.d0 / dble(nm)   !- t方向の刻み幅の設定, 1.d0 は倍精度 1. x 10^0 を明示
  x = 0.5d0              !- 初期値
  t = 0.d0               !- 初期のt

  open(10,file='ode.dat',form='unformatted',access='stream')
  write(10) nm
  write(10) dt

  write(6,*) '                   t,                       x,                      true,                     error'
  do i = 0, nm-1
     x = x + dt * f(x,t) ! オイラー法
     t = t + dt
     xa = t**3 + 1.d0 / ( 1.d0 + exp(-t) )  !- 解析解(真値)
     write(6,*) t,x,xa,x-xa
     write(10) x,xa
  enddo
  close(10)

contains          ! 以下の関数副プログラムを内部関数とする
  function f(x,t) ! 関数副プログラムの定義
    real(8)            :: f
    real(8),intent(in) :: t, x
    f = 3.d0*t**2 - t**3 - t**6 + ( 2.d0*t**3 + 1.d0 )*x - x**2
  end function f

end program ode
