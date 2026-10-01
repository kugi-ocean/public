!## 線形近似した状態方程式 eos を提供する
module linear_equation_of_state

  implicit none

  private ! 標準でprivate属性とする
          ! （変数はこのモジュールでのみ参照可能）
  public eos ! eos のみpublic属性を与える
             ! （eosはこのモジュールを呼び出す親プログラムで参照可能）

  real(8),parameter:: &
       t2r =1.65d-4, & ! 熱膨張率 temp=10., salt=34.
       s2r =7.61d-4, & ! 塩収縮率 temp=10., salt=34.
       !s2r = 0.0, & ! 塩収縮率 temp=10., salt=34.
       t0 = 10.d0, & ! 基準温位
       s0 = 34.d0, & ! 基準塩分（34.0 psuからの偏差とする）
       r0 = 1025.d0
       
contains ! 以下は内部プログラム
  function eos(t,s,p)

    real(8):: &
         t, s, p, eos
    
    eos=r0-r0*t2r*(t-t0)+r0*s2r*(s-s0)
  end function eos
end module linear_equation_of_state



!## 近似していない非線形の状態方程式 eos を提供する
module nonlinear_equation_of_state
  !     EQUATION OF STATE FOR SEAWATER AS A FUNCTION OF
  !     POTENTIAL TEMPERATURE, SALINITY, AND PRESSURE
  !     
  !     Jackett et al. (2006) J. ATM. OCEAN. TECH.
  !     McDougall et al. (2003) J. ATM. OCEAN. TECH.
  !
  !     NOTE
  !     potential density referred at p=PR is also calculated
  !     from this subroutine with replacing p with PR

  implicit none
  private
  public eos
  
  real(8),parameter:: &
       NC0   =  9.9984085444849347d+2, &
       NT1   =  7.3471625860981584d+0, &
       NT2   = -5.3211231792841769d-2, &
       NT3   =  3.6492439109814549d-4, &
       NS1   =  2.5880571023991390d+0, &
       NS1T1 = -6.7168282786692354d-3, &
       NS2   =  1.9203202055760151d-3, &
       NP1   =  1.1798263740430364d-2, &
       NP1T2 =  9.8920219266399117d-8, &
       NP1S1 =  4.6996642771754730d-6, &
       NP2   = -2.5862187075154352d-8, &
       NP2T2 = -3.2921414007960662d-12

  real(8),parameter:: &
       DC0   =  1.0, &
       DT1   =  7.2815210113327091d-3, &
       DT2   = -4.4787265461983921d-5, &
       DT3   =  3.3851002965802430d-7, &
       DT4   =  1.3651202389758572d-10, &
       DS1   =  1.7632126669040377d-3, &
       DS1T1 = -8.8066583251206474d-6, &
       DS1T3 = -1.8832689434804897d-10, &
       DSx   =  5.7463776745432097d-6, &
       DSxT2 =  1.4716275472242334d-9, &
       DP1   =  6.7103246285651894d-6, &
       DP2T3 = -2.4461698007024582d-17, &
       DP3T1 = -9.1534417604289062d-18

  real(8):: &
       pn,pd


contains

  function eos(t,s,p)

    implicit none

    real(8):: &
         t,s,p,eos

    pn   =NC0 &
         +NT1  *t  &
         +NT2  *t**2 &
         +NT3  *t**3 &
         +NS1  *s &
         +NS1T1*s*t &
         +NS2  *s**2 &
         +NP1  *p &
         +NP1T2*p*t**2 &
         +NP1S1*p*s &
         +NP2  *p**2 &
         +NP2T2*p**2*t**2
    
    pd   =DC0 &
         +DT1  *t &
         +DT2  *t**2 &
         +DT3  *t**3 &
         +DT4  *t**4 &
         +DS1  *s &
         +DS1T1*s*t &
         +DS1T3*s*t**3 &
         +DSx  *s**1.5 &
         +DSxT2*s**1.5*t**2 &
         +DP1  *p &
         +DP2T3*p**2*t**3 &
         +DP3T1*p**3*t

    eos=pn/pd

  end function eos

end module nonlinear_equation_of_state
