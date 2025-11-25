!- 3次元密度一様モデル
!-   Stommel (1948) の設定で実験する場合 (他はデフォルト設定)
!-   im = 100, jm = 62, km = 1
!-   dx = 100km, dy = 100km, dz = 200m
!-   fric_non = 0.02 / dz
program main
  implicit none

  character(10):: &
       expid='ex3-case1' ! 実験名

  integer,parameter:: &
       im = 100, & ! x方向格子数
       jm =  62, & ! y方向格子数
       km =  1   ! z方向格子数

  real(8),parameter:: &
       t0_npm2  = -0.1d0, & ! 東西方向風応力
       gr_mps2  = 9.80d0, & ! 重力加速度
       r0_kgpm3 = 1.d3, & ! 基準海水密度
       f0_psec  = 0.d-4, & ! コリオリ係数
       bt_pmsec = 1.d-11, & ! ベータ
       er_m     = 6370.d3,& ! 地球半径
       om_psec  = 7.29d-5,& ! 地球自転速度
       ah_m2ps  = 1.d5, & ! 水平渦粘性係数
       av_m2ps  = 1.d-3, & ! 鉛直渦粘性係数
       dx_m     = 100.d3, & ! x方向格子幅
       dy_m     = 100.d3, & ! y方向格子幅
       dz_m     = 200.d0, & ! z方向格子幅
       dt_sec   = 100.d0, & ! 時間刻み幅
       time_to_start_sec   = 0.d0, & ! 実験開始時刻
       time_to_end_sec     = 86400.d0*30.d0*12.d0, &   ! 実験終了時刻（１年）
       !time_to_end_sec     = dt_sec * 10.d0, & ! 実験終了時刻（１年）
       output_interval_sec = 86400.d0*10.d0,       & ! 出力時間間隔（１０日）
       !output_interval_sec = dt_sec , & ! 出力時間間隔（１０日）
       asf      = 0.5d0 ! アセリンフィルター係数
  real(8),parameter :: fric_non = 0.0d0 / dz_m  !- linear friction coefficient (for Stommel 1948)

  real(8),dimension(:,:,:):: &
       ua_mps(0:im+1,0:jm+1,0:km+1), & ! 東西流速
       ub_mps(0:im+1,0:jm+1,0:km+1), &
       uc_mps(0:im+1,0:jm+1,0:km+1), &
       va_mps(0:im+1,0:jm+1,0:km+1), & ! 南北流速
       vb_mps(0:im+1,0:jm+1,0:km+1), &
       vc_mps(0:im+1,0:jm+1,0:km+1)
  
  real(8),dimension(:,:):: & 
       ea_m(0:im+1,0:jm+1), & ! 水位
       eb_m(0:im+1,0:jm+1), &
       ec_m(0:im+1,0:jm+1)

  real(8),dimension(:,:,:):: &
       ww_mps(0:im+1,0:jm+1,0:km+1), & ! 鉛直流速
       pp_npm2(0:im+1,0:jm+1,0:km+1)    ! 圧力

  real(8),dimension(:,:):: &
       tx_npm2(0:im+1,0:jm+1), & ! 風応力（東西）
       ty_npm2(0:im+1,0:jm+1)    ! 風応力（南北）

  real(8),dimension(:):: &
       fs_psec(0:jm+1) ! コリオリ係数

  real(8):: &
       time_sec, & ! 時間
       time_to_output_sec, & ! 次期出力時刻
       dpi, & ! π
       gu, gv, ge, gt, gs, &
       pre, adx, ady, adz, cor, dfx, dfy, dfz, frc
  
  integer:: &
       i,j,k,n ! 作業用

  character(80):: &
       buff ! 作業用

  ! コリオリ係数
  do j = 0, jm+1
     fs_psec(j) = 
     !fs_psec(j) = 2.d0 * om_psec * sin( (dble(j)-0.5d0) * dy_m /er_m )
  end do

  ! 強制力
  dpi = 4.d0 * atan(1.d0)
  do j = 1, jm
     do i = 1, im
        tx_npm2(i,j) =
        ty_npm2(i,j) =
     end do
  end do

  ! 初期化
  uc_mps(:,:,:) = 0.d0
  vc_mps(:,:,:) = 0.d0
  ec_m(:,:)     = 0.d0

  ww_mps(:,:,:)  = 0.d0
  pp_npm2(:,:,:) = 0.d0

  ! 初期設定
  if( time_to_start_sec == 0.d0 ) then ! 初期化
     ua_mps(:,:,:) =
     ub_mps(:,:,:) =
     va_mps(:,:,:) =
     vb_mps(:,:,:) =
     ea_m(:,:)     =
     eb_m(:,:)     =
     time_sec      =

  else ! 継続

     open(10,file=trim(expid)//'.cnt',form='unformatted',access='stream')
     read(10) time_sec,ua_mps,ub_mps,va_mps,vb_mps,ea_m,eb_m
     close(10)

  endif

  time_to_output_sec = time_sec + output_interval_sec

  do
     if( time_sec > time_to_end_sec ) exit

     ! 診断量の計算
     ! w
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              ww_mps(i,j,k) = 
           end do
        end do
     end do
     ww_mps(:,:,km) = 0.d0

     ! p
     do j = 1, jm
        do i = 1, im
           pp_npm2(i,j,km) =
           do k = km-1, 1, -1
              pp_npm2(i,j,k) =
           end do
        end do
     end do

     ! 予報量の計算
     ! u
     do k = 1, km
        do j = 1, jm
           do i = 1, im - 1
              uc_mps(i,j,k) =
           end do
        end do
     end do
     
     ! v
     do k = 1, km
        do j = 1, jm - 1
           do i = 1, im
              vc_mps(i,j,k) =
           end do
        end do
     end do

     ! e
     do j = 1, jm
        do i = 1, im
           ec_m(i,j) =
        end do
     end do

     ! 東西境界条件
     uc_mps(0   ,:,:) =
     uc_mps(im  ,:,:) =
     uc_mps(im+1,:,:) = 0.d0  !- dummy

     vc_mps(0   ,:,:) =
     vc_mps(im+1,:,:) =

     ec_m(0   ,:)     =
     ec_m(im+1,:)     =

     ! 南北境界条件
     uc_mps(:,0   ,:) =
     uc_mps(:,jm+1,:) =

     vc_mps(:,0   ,:) =
     vc_mps(:,jm  ,:) =
     vc_mps(:,jm+1,:) = 0.d0

     ec_m(:,0   )     =
     ec_m(:,jm+1)     =

     ! 上下境界条件
     uc_mps(:,:,km+1) =
     uc_mps(:,:,0   ) =

     vc_mps(:,:,km+1) =
     vc_mps(:,:,0   ) =

     ! n+1ステップが求まったので配列をシフトさせる
     ! ua_mps=ub_mps, ub_mps=uc_mps 
     ua_mps(:,:,:)    = ub_mps(:,:,:) +0.5d0*asf* ( ua_mps(:,:,:) - 2.d0*ub_mps(:,:,:) + uc_mps(:,:,:) )
     va_mps(:,:,:)    = vb_mps(:,:,:) +0.5d0*asf* ( va_mps(:,:,:) - 2.d0*vb_mps(:,:,:) + vc_mps(:,:,:) )
     ea_m(:,:)        = eb_m(:,:)     +0.5d0*asf* ( ea_m(:,:)     - 2.d0*eb_m(:,:)     + ec_m(:,:)     )

     ! n+1 step -> ub_mps,vb_mps,eb_m,tb,sb
     ub_mps(:,:,:)    = uc_mps(:,:,:)
     vb_mps(:,:,:)    = vc_mps(:,:,:)
     eb_m(:,:)        = ec_m(:,:)

     time_sec = time_sec + dt_sec

     if( time_sec < time_to_output_sec ) cycle

     n = idnint( time_sec / output_interval_sec )

     ! 診断量の再計算（記録用）
     ! w
     do k = 1, km
        do j = 1, jm
           do i = 1, im
              ww_mps(i,j,k) =
           end do
        end do
     end do
     ww_mps(:,:,km) = 0.d0

     ! p
     do j = 1, jm
        do i = 1, im
           pp_npm2(i,j,km) =
           do k = km-1, 1, -1
              pp_npm2(i,j,k) =
           end do
        end do
     end do

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
     write(6,*) 'umin: ',minval(ub_mps(1:im,1:jm,1:km)),'(',minloc(ub_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'umax: ',maxval(ub_mps(1:im,1:jm,1:km)),'(',maxloc(ub_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'vmin: ',minval(vb_mps(1:im,1:jm,1:km)),'(',minloc(vb_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'vmax: ',maxval(vb_mps(1:im,1:jm,1:km)),'(',maxloc(vb_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'wmin: ',minval(ww_mps(1:im,1:jm,1:km)),'(',minloc(ww_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'wmax: ',maxval(ww_mps(1:im,1:jm,1:km)),'(',maxloc(ww_mps(1:im,1:jm,1:km)),')'
     write(6,*) 'emin: ',minval(eb_m(1:im,1:jm)),'(',minloc(eb_m(1:im,1:jm)),')'
     write(6,*) 'emax: ',maxval(eb_m(1:im,1:jm)),'(',maxloc(eb_m(1:im,1:jm)),')'

     time_to_output_sec = time_to_output_sec + output_interval_sec

  end do

end program main
