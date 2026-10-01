!# 海水のポテンシャル密度を計算する
program sample

!  use linear_equation_of_state, only: eos
  use nonlinear_equation_of_state, only: eos
  use profile,                     only: profile_data, print_profile
  implicit none  !- 変数は必ず宣言して使う、と設定する
  
  type(profile_data) :: A, B
  real(4)            :: x,y
  real(8)            :: d,t,s
  integer            :: k


  !##  書式つきファイル A.dat から入力データを読み込んで密度を計算する
  open(10,file='A.dat',form='formatted')  !- 10番ゲートにファイルを開く
  read(10,'(2f8.1)') x,y  !- 10番ゲートで8桁の数（小数点以下1桁）を2つ読む
  A%lat=x
  A%lon=y
  do k = 1, 102
    read(10,'(3f7.1)') d,t,s  !- 7桁の数（小数点以下1桁）を3つ読む
    if( t == -999. ) exit      !- 欠損値 -999. が水温に入っていれば入力を終了
    A%depth(k)=d
    A%temp (k)=t
    A%salt (k)=s
    A%dens (k)=eos(A%temp(k),A%salt(k),0.d0)  !- モジュール中の関数を使って水深0mでの密度を得る
  enddo
  A%num=k-1  !- do loopからexitで出たときの k から1を引いてデータの数
  A%station_name = 'Sta. A'


  !##  書式つきファイル B.dat も計算
  open(10,file='B.dat',form='formatted')
  read(10,'(2f8.1)') x,y
  B%lat=x
  B%lon=y
  do k = 1, 102
    read(10,'(3f7.1)') d,t,s
    if( t == -999. ) exit
    B%depth(k)=d
    B%temp (k)=t
    B%salt (k)=s
    B%dens (k)=eos(B%temp(k),B%salt(k),0.d0)
  enddo
  B%num=k-1
  B%station_name = 'Sta. B'


  !## 画面出力 (モジュール profile のサブルーチン print_profile を使う)
  call print_profile(A)
  call print_profile(B)

end program sample
