!
!     water hammer code
!
!     original code : hammer.for @ hydraulic foumula program s3-2
!     水理公式集例題プログラム集平成13年版【例題3-2】を改変
!
!*****************************************************************************
    module common_variables
      implicit none
      integer, parameter :: k = kind(1.0d0)
!
      real(k), parameter :: gra = 9.8                ! 重力加速度
      real(k), parameter :: pi = 4 * atan(1.)        ! 円周率 3.141592654
!
      real(k), parameter :: h0 = 160.0               ! 貯水池水位
      real(k), parameter :: al = 400.0               ! 管路長
      real(k), parameter :: dp  =   2.0              ! 管路内径
      real(k), parameter :: ff  =   0.01             ! 管路摩擦損失係数
!
      real(k), parameter :: vs =   3.14              ! 初期流速
      real(k), parameter :: tv =   1.8               ! バルブ閉塞時間
!
      real(k), parameter :: cc   = 1000.0            ! 水中音速
!
      integer, parameter :: nxx=500                  ! 距離分割数の配列の上限
      integer, parameter :: nixx=10                  ! 途中経過記録格子点数の配列の上限
!
      integer, parameter :: ntau = 6                 ! 計算総周期数（何周期分計算するか）
      integer, parameter :: nx = 400                 ! 距離分割数
      integer, parameter :: nix = 4                  ! 途中経過を記録する最終の格子点番号
!
      real(k) he, qs, tau, hs, bb, rr
      real(k) dt, dx
      integer nt

      real(k) xx(0:nxx), hh(0:nxx), h1(0:nxx), hmax(0:nxx), hmin(0:nxx), qq(0:nxx), q1(0:nxx)
      integer ix(0:nixx)
!
    end module common_variables
!
!*****************************************************************************
    program hammer
      use common_variables
      implicit none
      integer :: i, j, j1
      real(k) :: t1
!
!                                                  出力用ファイルのオープン
      open(1,file='hammer1.dat',status='unknown')
      open(2,file='hammer2.dat',status='unknown')
      open(3,file='hammer3.dat',status='unknown')
!
      call input                      ! 計算条件および初期条件
      call output(0, 0.0)
!                                       時間について繰返し
      do 10 j = 0, ntau*nt-1
          j1 = j+1                                   ! 次の時刻
          t1 = j1*dt
          call hammer_cal(t1)                        ! 水撃圧の計算
          call output(j1,t1)                         ! ファイル出力
   10 continue
      
      write(3,100) (i,xx(i),hmax(i),hmin(i), i=0,nx)
  100 format(i6, 3f12.3)
      close(1)
      close(2)
      close(3)
!
      stop
!
    end ! program hammer

!*****************************************************************************
!
!  subroutine: input
!
!  計算条件および初期条件の設定
!
!*****************************************************************************
    subroutine input
      use common_variables
      implicit none
      integer :: i, iw
      real(k) :: ts, rho, kw, es
      real(k) :: aa
!
      ix(0) =   0     ! 途中経過を記録する格子点番号
      ix(1) = 100
      ix(2) = 200
      ix(3) = 300
      ix(4) = 400
!
!      hs    :バルブにおける初期ピエゾ水頭
!
!      x     :配列:格子点位置
!      ix    :配列:途中経過を記録する格子点の番号
!      h,h1  :配列:各格子点でのある時刻および次の時刻のピエゾ水頭
!      q,q1  :配列:各格子点でのある時刻および次の時刻の流量

      tau = 2 * al / cc                 ! 水撃波の往復時間
      dx = al / nx                     ! 距離刻み
      dt = dx / cc                      ! 時間刻み(=dx/cc=tau/2/nx)
      nt = 2 * nx                      ! １周期の時間分割数(=2*nx)
!
      aa = pi / 4 * dp**2                ! 管路断面積
      bb = cc / gra / aa                    ! 係数B
      rr = ff * cc / 2 / gra / dp / aa**2     ! 係数R
!
      qs = vs * aa                      ! 初期流量
!
!                                      格子点の設定
          xx(0) = -al
      do 10 i = 1, nx-1
          xx(i) = i * dx - al
   10 continue
          xx(nx) = 0.0
!
!                                      初期値の設定
      do 20 i = 0, nx
          qq(i) = qs
          hh(i) = h0 - rr/cc*qs*abs(qs)*(xx(i)+al)
          hmax(i) = hh(i)
          hmin(i) = hh(i)
   20 continue
!
      hs = hh(nx)
!
!                             計算条件のファイル出力
      do 30 iw = 1, 3
          write(iw,100) '貯水池水位', h0
          write(iw,100) '落差', he
          write(iw,100) '管路長', al
          write(iw,100) '管路内径', dp
          write(iw,100) '管路断面積', aa
          write(iw,100) '管路摩擦損失係数', ff
          write(iw,100) '初期流量', qs
          write(iw,100) '初期流速', vs
          write(iw,100) '伝播速度', cc
          write(iw,100) 'バルブ閉塞時間', tv
          write(iw,100) '水撃波の往復時間', tau
          write(iw,*)
          write(iw,110) '距離分割数', nx
          write(iw,100) '距離刻み', dx
          write(iw,*)
          write(iw,110) '計算総周期数', ntau
          write(iw,110) '１周期の時間分割数', nt
          write(iw,100) '時間刻み', dt
          write(iw,*)
   30 continue
!
      write(1,110) '経過記録格子点数', nix+1
      write(1,200) 'IX', 'X'
      write(1,210) (ix(i),xx(ix(i)), i=0,nix)
      write(1,*)
      write(1,200) 'J', 'T', ('H',i,'Q',i, i=0,nix)
!
      write(2,300) 'J', 'T', 'I', 'X', 'H', 'Q'
      
      write(3,300) 'I', 'X', 'HMAX'
!
  100 format(a, ':', t20, f12.3)
  110 format(a, ':', t20, i12)
  200 format(t3,a, t13,a, 100(9x,a,i2,9x,a,i2))
  210 format(i6,f12.3)
  300 format(t3,a, t12,a, t21,a, t30,a, t42,a, t54,a) 
!
      return

    end subroutine input
!
!*****************************************************************************
!
!  subroutine: hammer_cal
!
!  水撃圧を１ステップ計算
!
!*****************************************************************************
    subroutine hammer_cal(t1)
      use common_variables
      implicit none
      real(k), intent(in) :: t1
      integer :: i
      real(k) :: phy, ak1, ak2
!
!                                                  上流端条件
      h1(0) = h0
      q1(0) = qq(1) + (h1(0)- hh(1))/bb - rr/bb*dt*qq(1)*abs(qq(1))
!
      do 10 i = 1, nx-1                           ! x方向に繰返し
!
          h1(i) = (hh(i-1)+hh(i+1))/2 + bb*(qq(i-1)-qq(i+1))/2 + rr*dt*(qq(i+1)*abs(qq(i+1))-qq(i-1)*abs(qq(i-1)))/2
          q1(i) = (qq(i-1)+qq(i+1))/2 + (hh(i-1)-hh(i+1))/bb/2 - rr*dt*(qq(i+1)*abs(qq(i+1))+qq(i-1)*abs(qq(i-1)))/bb/2
!
   10 continue
!                                                  下流端条件
      if (t1.lt.tv) then      ! バルブ開度
          phy = 1 - t1/tv
      else
          phy =0.
      endif
!
      ak1=hh(nx-1)+bb*qq(nx-1)-rr*dt*qq(nx-1)*abs(qq(nx-1))
      ak2=(phy*bb*qs)**2 /2/hs
      h1(nx)=ak1+ak2-(ak2*(ak2+2*ak1))**0.5
      q1(nx) = phy*qs*h1(nx)/(abs(hs * h1(nx)))**0.5
!
      do 20 i=0,nx                                ! 値の更新
          hh(i)=h1(i)
          qq(i)=q1(i)
          hmax(i) = max(hh(i), hmax(i))
          hmin(i) = min(hh(i), hmin(i))
   20 continue
!
      return
!
    end subroutine hammer_cal

!*****************************************************************************
!
!  subroutine: output
!
!  計算結果のファイル出力
!
!*****************************************************************************
    subroutine output(j,t)
      use common_variables
      implicit none
      integer, intent(in) :: j
      real(k), intent(in) :: t
      integer :: i
!
!                                                     ixで指定された格子点について
      write(1,100) j, t, (hh(ix(i)),qq(ix(i)), i=0,nix)
!
      if (mod(j, nt/2).ne.0) return
!                                                     全格子点について
      do 10 i = 0, nx
          write(2,200) j, t, i, xx(i), hh(i), qq(i)
   10 continue
          write(2,*)
!
  100 format(i6, f12.5, 100(2f12.3))
  200 format(i6, f12.5, i6, 3f12.3) 
!
      return
!
    end subroutine output
