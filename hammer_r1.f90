!
!     water hammer code
!
!     original code : hammer.for @ hydraulic foumula program s3-2
!     ���������W���v���O�����W����13�N�Ły���3-2�z������
!
!*****************************************************************************
    module common_variables
      implicit none
      integer, parameter :: k = kind(1.0d0)
!
      real(k), parameter :: gra = 9.8                ! �d�͉����x
      real(k), parameter :: pi = 4 * atan(1.)        ! �~���� 3.141592654
!
      real(k), parameter :: h0 = 160.0               ! �����r����
      real(k), parameter :: al = 400.0               ! �ǘH��
      real(k), parameter :: dp  =   2.0              ! �ǘH���a
      real(k), parameter :: ff  =   0.01             ! �ǘH���C�����W��
!
      real(k), parameter :: vs =   3.14              ! ��������
      real(k), parameter :: tv =   1.8               ! �o���u�ǎ���
!
      real(k), parameter :: cc   = 1000.0            ! ��������
!
      integer, parameter :: ntau = 6                 ! �v�Z���������i���������v�Z���邩�j
      integer, parameter :: nx = 400                 ! ����������
      integer, parameter :: nix = 4                  ! �r���o�߂��L�^����i�q�_��
!
      real(k) :: qs, hs, bb, rr, dt, dx
      integer :: nt
!
    end module common_variables
!
!*****************************************************************************
    program hammer
      use common_variables
      implicit none
      real(k), dimension(0:nx) :: xx, hh, hmax, hmin, qq
      integer, dimension(0:nix) :: ix
!      xx     :�z��:�i�q�_�ʒu
!      ix     :�z��:�r���o�߂��L�^����i�q�_�̔ԍ�
!      hh     :�z��:�e�i�q�_�ł̂��鎞������ю��̎����̃s�G�]����
!      qq     :�z��:�e�i�q�_�ł̂��鎞������ю��̎����̗���
!
      integer :: i, j, j1
      real(k) :: t1
!
      open(1,file='hammer1.dat',status='unknown')                 ! �S�����E�r���o�ߋL�^�i�q�_�̏o��
      open(2,file='hammer2.dat',status='unknown')                 ! �������E�S�i�q�_�̏o��
      open(3,file='hammer3.dat',status='unknown')                 ! �S�i�q�_�̍ő�E�ŏ��s�G�]�����̏o��
!
      call input(xx, hh, hmax, hmin, qq, ix)                      ! �v�Z��������я�������
      call output(0, 0.0, xx, hh, qq, ix)                         ! �t�@�C���o��
!
      do j = 0, ntau*nt-1                                         ! ���Ԃɂ��ČJ�Ԃ�
          j1 = j+1                                                ! ���̎���
          t1 = j1*dt
          call hammer_cal(t1, hh, hmax, hmin, qq)                 ! �������̌v�Z
          call output(j1, t1, xx, hh, qq, ix)                     ! �t�@�C���o��
      end do
      
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
!  �v�Z��������я��������̐ݒ�
!
!*****************************************************************************
    subroutine input(xx, hh, hmax, hmin, qq, ix)
      use common_variables
      implicit none
      real(k), dimension(0:nx), intent(out) :: xx, hh, hmax, hmin, qq
      integer, dimension(0:nix), intent(out) :: ix
!
      integer :: i, iw
      real(k) :: ts, rho, kw, es
      real(k) :: aa, tau
!
!     �r���o�߂��L�^����i�q�_�ԍ�
      ix(0) =   0     
      ix(1) = 100
      ix(2) = 200
      ix(3) = 300
      ix(4) = 400
!
      tau = 2 * al / cc                     ! �����g�̉�������
      dx = al / nx                          ! ��������
      dt = dx / cc                          ! ���ԍ���(=dx/cc=tau/2/nx)
      nt = 2 * nx                           ! �P�����̎��ԕ�����(=2*nx)
!
      aa = pi / 4 * dp**2                   ! �ǘH�f�ʐ�
      bb = cc / gra / aa                    ! �W��B
      rr = ff * cc / 2 / gra / dp / aa**2   ! �W��R
!
      qs = vs * aa                          ! ��������
!
!     �i�q�_�̐ݒ�
          xx(0) = -al
      do i = 1, nx-1
          xx(i) = i * dx - al
      end do
          xx(nx) = 0.0
!
!     �����l�̐ݒ�
      do i = 0, nx
          qq(i) = qs
          hh(i) = h0 - rr/cc*qs*abs(qs)*(xx(i)+al)
          hmax(i) = hh(i)
          hmin(i) = hh(i)
      end do
!
      hs = hh(nx)                           ! �o���u�ɂ����鏉���s�G�]����
!
!     �v�Z�����̃t�@�C���o��
      do iw = 1, 3
          write(iw,100) '�����r����', h0
          write(iw,100) '�ǘH��', al
          write(iw,100) '�ǘH���a', dp
          write(iw,100) '�ǘH�f�ʐ�', aa
          write(iw,100) '�ǘH���C�����W��', ff
          write(iw,100) '��������', qs
          write(iw,100) '��������', vs
          write(iw,100) '�`�d���x', cc
          write(iw,100) '�o���u�ǎ���', tv
          write(iw,100) '�����g�̉�������', tau
          write(iw,*)
          write(iw,110) '����������', nx
          write(iw,100) '��������', dx
          write(iw,*)
          write(iw,110) '�v�Z��������', ntau
          write(iw,110) '�P�����̎��ԕ�����', nt
          write(iw,100) '���ԍ���', dt
          write(iw,*)
      end do
!
      write(1,110) '�o�ߋL�^�i�q�_��', nix+1
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
!  ���������P�X�e�b�v�v�Z
!
!*****************************************************************************
    subroutine hammer_cal(t1, hh, hmax, hmin, qq)
      use common_variables
      implicit none
      real(k), intent(in) :: t1
      real(k), dimension(0:nx), intent(inout) :: hh, hmax, hmin, qq
!
      integer :: i
      real(k), dimension(0:nx) :: h1, q1
      real(k) :: phy, ak1, ak2
!
!     �㗬�[����
      h1(0) = h0
      q1(0) = qq(1) + (h1(0)- hh(1))/bb - rr/bb*dt*qq(1)*abs(qq(1))
!
      do i = 1, nx-1                           ! x�����ɌJ�Ԃ�
          h1(i) = (hh(i-1)+hh(i+1))/2 + bb*(qq(i-1)-qq(i+1))/2 + rr*dt*(qq(i+1)*abs(qq(i+1))-qq(i-1)*abs(qq(i-1)))/2
          q1(i) = (qq(i-1)+qq(i+1))/2 + (hh(i-1)-hh(i+1))/bb/2 - rr*dt*(qq(i+1)*abs(qq(i+1))+qq(i-1)*abs(qq(i-1)))/bb/2
      end do
!
!     �����[����
      if (t1.lt.tv) then                       ! �o���u�J�x
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
      do i=0,nx                                ! �l�̍X�V
          hh(i)=h1(i)
          qq(i)=q1(i)
          hmax(i) = max(hh(i), hmax(i))
          hmin(i) = min(hh(i), hmin(i))
      end do
!
      return
!
    end subroutine hammer_cal

!*****************************************************************************
!
!  subroutine: output
!
!  �v�Z���ʂ̃t�@�C���o��
!
!*****************************************************************************
    subroutine output(j, t, xx, hh, qq, ix)
      use common_variables
      implicit none
      real(k), dimension(0:nx), intent(in) :: xx, hh, qq
      integer, dimension(0:nix), intent(in) :: ix
!
      integer, intent(in) :: j
      real(k), intent(in) :: t
      integer :: i
!
      write(1,100) j, t, (hh(ix(i)),qq(ix(i)), i=0,nix)
!
      if (mod(j, nt/2).ne.0) return
!
      do i = 0, nx
          write(2,200) j, t, i, xx(i), hh(i), qq(i)
      end do
          write(2,*)
!
  100 format(i6, f12.5, 100(2f12.3))
  200 format(i6, f12.5, i6, 3f12.3) 
!
      return
!
    end subroutine output
