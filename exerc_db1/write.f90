!- 配列データの出力をテストする
program writetest

  implicit none                            ! 暗黙の型宣言の抑止
  integer(8),parameter :: im=2,jm=3,km=4   ! 8バイトの整数型の定数の定義
  integer(8) :: i,j,k                      ! 8バイトの整数型の変数の定義
  real(8)    :: a(0:im-1,0:jm-1,0:km-1)    ! 倍精度実数型の3次元配列の定義
  integer(8) :: b(0:im-1,0:jm-1,0:km-1)    ! 8バイトの整数型の3次元配列の定義

  do i = 0, im-1
     do j = 0, jm-1
        do k = 0, km-1        
           a(i,j,k) = dble( 100*i + 10*j + k )   ! 3次元配列に倍精度実数を代入
        end do
     end do
  end do

  do i = 0, im-1
     do j = 0, jm-1
        do k = 0, km-1        
           b(i,j,k) = - ( 100*i + 10*j + k )     ! 3次元配列に整数を代入
        end do
     end do
  end do

  open(10,file='data.dat',form='unformatted',access='stream') ! 書き込み用ファイルを10番に開く
  write(10) a   ! 配列aを書き込み
  write(10) b   ! 配列bを書き込み
  close(10)     ! 書き込み用ファイルを閉じる

end program writetest
