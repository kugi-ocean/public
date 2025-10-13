!# 鉛直1次元モデル：海面エクマン層 (MKS単位)
program main
  implicit none

  character(10):: &
       expid='ex1-case1' ! 実験名

  integer,parameter:: &
       km =    ! z方向格子数

  real(8),parameter:: &
       tx_npm2  = , & ! 東西方向風応力
       ty_npm2  = , & ! 南北方向風応力
       r0_kgpm3 = 1.d3, & ! 基準海水密度
       f0_psec  = , & ! コリオリ係数
       av_m2ps  = 1.d-2 , & ! 鉛直粘性係数
       dz_m     = , & ! 格子幅
       dt_sec   = , & ! 時間刻み幅
       time_to_end_sec     = 86400.d0 * 10.d0 , & ! 実験終了時刻（１０日）
       output_interval_sec = 3600.d0 * 1.d0 , & ! 出力時間間隔（１時間）
       af       = 0.5d0 ! アセリンフィルター係数

  real(8) :: &
       ! 東西流速
       ua_mps(0:km+1), & ! n-1
       ub_mps(0:km+1), & ! n
       uc_mps(0:km+1), & ! n+1
       ! 南北流速
       va_mps(0:km+1), & ! n-1
       vb_mps(0:km+1), & ! n
       vc_mps(0:km+1)    ! n+1

  ! 作業用
  real(8):: &
       time_sec, & ! 時間
       time_to_output_sec, & ! 次期出力時刻
       gu, gv
  integer:: &
       i,j,k,n
  character(80):: &
       buff

  ! 初期条件設定
  ua_mps(:)=0.d0
  ub_mps(:)=0.d0
  va_mps(:)=0.d0
  vb_mps(:)=0.d0

  time_sec = 0.d0
  time_to_output_sec = time_sec + output_interval_sec

  do
     if( time_sec > time_to_end_sec ) exit

     ! 予報変数の計算
     ! u
     do k=1,km
        gu = 
        uc_mps(k) =
     end do

     ! v
     do k=1,km
        gv =
        vc_mps(k) =
     end do

     ! 境界条件
     uc_mps(km+1) =
     uc_mps(0   ) =

     vc_mps(km+1) =
     vc_mps(0   ) =

     ! n+1ステップが求まったので、次の時間ステップに向けて配列をシフトさせる
     ! 今回                       次回
     ! n-1 step :ua_mps,va_mps
     ! n   step :ub_mps,vb_mps => ua_mps,va_mps
     ! n+1 step :uc_mps,vc_mps => ub_mps,vb_mps

     ! ua_mps<=ub_mps - アセリンフィルター
     ua_mps(:) = ub_mps(:) + 0.5d0 * af * ( ua_mps(:) -2.d0 * ub_mps(:) + uc_mps(:) )
     va_mps(:) = vb_mps(:) + 0.5d0 * af * ( va_mps(:) -2.d0 * vb_mps(:) + vc_mps(:) )
     ! ub_mps<=uc_mps
     ub_mps(:) = uc_mps(:)
     vb_mps(:) = vc_mps(:)

     time_sec = time_sec + dt_sec

     ! time が time_to_output より前だったら次のステップへ(過ぎていたら値を保存)
     if( time_sec < time_to_output_sec ) cycle

     n = idnint( time_sec / output_interval_sec )

     ! 記録（解析・作図用）
     write(buff,'(a,i6.6)') trim(expid)//'.n',n
     open(10,file=trim(buff),form='unformatted',access='stream')
     write(10) time_sec, ub_mps, vb_mps
     close(10)

     ! 画面監視
     write(6,*) 'time [s] = ',time_sec
     write(6,*) 'umin: ',minval(ub_mps(1:km)),'(',minloc(ub_mps(1:km)),')'
     write(6,*) 'umax: ',maxval(ub_mps(1:km)),'(',maxloc(ub_mps(1:km)),')'
     write(6,*) 'vmin: ',minval(vb_mps(1:km)),'(',minloc(vb_mps(1:km)),')'
     write(6,*) 'vmax: ',maxval(vb_mps(1:km)),'(',maxloc(vb_mps(1:km)),')'

     ! 次の記録時間を設定
     time_to_output_sec = dble(n+1) * output_interval_sec
  end do

end program main
