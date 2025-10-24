!# 鉛直2次元モデル: 赤道上の風成循環 (MKS単位)
program main
  implicit none

  character(*),parameter ::  expid='ex2-case1' ! 実験名

  integer,parameter:: &
       im = 150, & ! x方向格子数
       km = 15    ! z方向格子数

  real(8),parameter:: &
       tx_npm2  = 0.1d0,  & ! 東西方向風応力 [N/m2]
       gr_mps2  = 9.8d0,  & ! 重力加速度     [m/s2]
       r0_kgpm3 = 1.d3,   & ! 基準海水密度   [kg/m3]
       f0_psec  = 0.d0,   & ! コリオリ係数   [/s]
       ah_m2ps  = 5.d4,   & ! 水平渦粘性係数 [m2/s]
       av_m2ps  = 1.d-2,  & ! 鉛直渦粘性係数 [m2/s]
       dx_m     = 10.0d3, & ! x方向格子幅 [m]
       dz_m     = 100.d0, & ! z方向格子幅 [m]
       dt_sec   = 30.d0,   & ! 時間刻み幅 [s]
       time_to_start_sec   = 0.d0, & ! 実験開始時刻
       time_to_end_sec     = 86400.d0*30.d0*12.d0, & ! 実験終了時刻（１年）
       !time_to_end_sec     = 60.d0, & ! 実験終了時刻（１年）
       output_interval_sec = 86400.d0*10.d0, & ! 出力時間間隔（１０日）
       !output_interval_sec = 30.d0*1.d0, & ! 出力時間間隔（１０日）
       asf                 = 0.5d0 ! アセリンフィルター係数

  real(8),dimension(:,:):: &
       ua_mps(0:im+1,0:km+1), & ! 東西流速 [m/s] at n-1 step
       ub_mps(0:im+1,0:km+1), & !                   n   step
       uc_mps(0:im+1,0:km+1), & !                   n+1 step
       va_mps(0:im+1,0:km+1), & ! 南北流速 [m/s] at n-1 step
       vb_mps(0:im+1,0:km+1), & !                   n   step
       vc_mps(0:im+1,0:km+1)    !                   n+1 step

  real(8),dimension(:,:):: & 
       ea_m(0:im+1), & ! 水位 [m]
       eb_m(0:im+1), &
       ec_m(0:im+1)

  real(8),dimension(:,:,:)::  &
       ww_mps(0:im+1,0:km+1), & ! 鉛直流速 [m/s]  (:,km+1)は使わない
       pp_npm2(0:im+1,0:km+1)   ! 圧力     [N/m2] (:,0),(:,km+1)は使わない

  real(8):: &
       time_sec, & ! 時間
       time_to_output_sec, & ! 次期出力時刻
       gu, gv, ge, &
       pre, adx, adz, cor, dfx, dfz, frc

  integer:: &
       i,k,n ! 作業用

  character(80):: &
       buff ! 作業用

  ! 初期化
  uc_mps(:,:)  = 0.d0
  vc_mps(:,:)  = 0.d0
  ec_m(:)      = 0.d0
  ww_mps(:,:)  = 0.d0
  pp_npm2(:,:) = 0.d0

  ! 初期設定
  if( time_to_start_sec == 0.d0 ) then ! 初期化
     ua_mps(:,:) =
     ub_mps(:,:) =
     va_mps(:,:) =
     vb_mps(:,:) =
     ea_m(:)     =
     eb_m(:)     =
     time_sec    =

  else ! 継続

     open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
     read(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m
     close(10)

  endif

  time_to_output_sec = time_sec + output_interval_sec

  !## Main loop
  do
     if( time_sec > time_to_end_sec ) exit

     ! 診断量の計算
     ! w
     do k = 1, km
        do i = 1, im
           ww_mps(i,k) =
        enddo
     enddo
     !ww_mps(:,km)=0.  上端でゼロとして

     ! p
     do i = 1, im
        do k =
           pp_npm2(i,k) =
        enddo
     enddo

     ! 予報量の計算
     ! u
     do k = 1, km
        do i = 1, im
           uc_mps(i,k) =
        enddo
     enddo


     ! v
     do k= 1, km
        do i = 1, im
           vc_mps(i,k) =
        enddo
     enddo

     ! eta
     do i= 1, im
        ec_m(i) = 
     enddo

     ! 東西境界条件
     uc_mps(0   ,:) =
     uc_mps(im  ,:) =

     vc_mps(0   ,:) =
     vc_mps(im+1,:) =

     ! 上下境界条件
     uc_mps(:,km+1) = 
     uc_mps(:,0   ) = 

     vc_mps(:,km+1) = 
     vc_mps(:,0   ) = 

     ! n+1ステップが求まったので配列をシフトさせる
     ! ua_mps=ub_mps, ub_mps=uc_mps 
     ua_mps(:,:) = ub_mps(:,:)+0.5d0*asf*(ua_mps(:,:)-2.d0*ub_mps(:,:)+uc_mps(:,:))
     va_mps(:,:) = vb_mps(:,:)+0.5d0*asf*(va_mps(:,:)-2.d0*vb_mps(:,:)+vc_mps(:,:))
     ea_m(:)     = eb_m(:)  +0.5d0*asf*(ea_m(:)  -2.d0*eb_m(:)  +ec_m(:)  )

     ! n+1 step -> ub_mps,vb_mps,eb_m,tb,sb
     ub_mps(:,:) = uc_mps(:,:)
     vb_mps(:,:) = vc_mps(:,:)
     eb_m(:)     = ec_m(:)

     time_sec = time_sec + dt_sec
     
     if( time_sec < time_to_output_sec ) cycle

     n = idnint( time_sec / output_interval_sec )

     ! 診断量の再計算（記録用）
     ! w
     do k = 1, km
        do i = 1, im
           ww_mps(i,k) = 
        enddo
     enddo
     !ww_mps(:,km)=0.

     ! p
     do i = 1, im
        do k = 
           pp_npm2(i,k) =
        enddo
     enddo

     ! 記録用（解析・作図用）
     write(buff,'(a,i6.6)') trim(expid)//'.n',n
     open(10,file=trim(buff),form='unformatted',access='stream')
     write(10) time_sec,ub_mps,vb_mps,eb_m,ww_mps,pp_npm2
     close(10)


     ! 継続用
     open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
     write(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m
     close(10)


     ! 画面で監視
     write(6,*) 'time [s] = ',time_sec
     write(6,*) 'umin: ',minval(ub_mps(1:im,1:km)),'(',minloc(ub_mps(1:im,1:km)),')'
     write(6,*) 'umax: ',maxval(ub_mps(1:im,1:km)),'(',maxloc(ub_mps(1:im,1:km)),')'
     write(6,*) 'vmin: ',minval(vb_mps(1:im,1:km)),'(',minloc(vb_mps(1:im,1:km)),')'
     write(6,*) 'vmax: ',maxval(vb_mps(1:im,1:km)),'(',maxloc(vb_mps(1:im,1:km)),')'
     write(6,*) 'wmin: ',minval(ww_mps(1:im,1:km)),'(',minloc(ww_mps(1:im,1:km)),')'
     write(6,*) 'wmax: ',maxval(ww_mps(1:im,1:km)),'(',maxloc(ww_mps(1:im,1:km)),')'
     write(6,*) 'emin: ',minval(eb_m(1:im)),'(',minloc(eb_m(1:im)),')'
     write(6,*) 'emax: ',maxval(eb_m(1:im)),'(',maxloc(eb_m(1:im)),')'

     time_to_output_sec = time_to_output_sec+output_interval_sec
  end do

end program main
