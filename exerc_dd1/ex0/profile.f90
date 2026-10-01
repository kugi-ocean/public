!# 鉛直プロファイルデータ
module profile

  implicit none
  private !- デフォルトでprivate属性とする
          ! （変数はこのモジュールでのみ参照可能）
  public profile_data, print_profile
          ! profile_dataとprint_profileのみpublic属性を与える
          ! （profile_dataとprint_profileはこのモジュールを
          !    呼び出す親プログラムで参照可能）

  type profile_data ! profile_dataという名前の型の定義
     character(8) :: station_name   !- 観測所の名前
     integer      :: num            !- 観測データ数
     real(4)      :: lat,lon
     real(8),dimension(102):: depth,temp,salt,dens
  end type profile_data

contains ! 内部プログラム

  subroutine print_profile(D)
  !- 構造体変数 D を受け取って、画面出力する。
  !-   (fortranでは大文字小文字を区別しないことに注意)

    type(profile_data):: D
    integer:: k
    
    write(6,*) 'Station : ',D%station_name
    write(6,*) 'Profile Num = ',D%num

    do k = 1, D%num
       write(6,'(4f7.1)') D%depth(k),D%temp(k),D%salt(k),D%dens(k)
    enddo

  end subroutine print_profile

end module profile
