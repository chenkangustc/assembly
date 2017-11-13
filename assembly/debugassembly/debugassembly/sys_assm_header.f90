module sys_assm_header
    implicit none
    type,public::AssmGeom
      real rFuel          !Ԫ���뾶
      real GasGap         !Ԫ����϶���
      real ShellThick     !Ԫ����Ǻ��
      real AssmShellThick !�����Ǻ��
      real AcrossFlat     !�����Ա߾ࣨ�������Ǻ�ȣ�
      real Height         !����߶ȣ���������
      integer n_pin       !Ԫ���ĸ���
    contains
      procedure,public::set=>set_assmgeom
    end type AssmGeom

    type,public::Assmmesh
        integer nf
        integer ng
        integer ns
        integer ny
      contains
      procedure,public::set=>set_assmmesh
    end type Assmmesh
    
    type,public::Assmboundary
       real Tin
       real Tout
       real uin
       real uout
       real pin
       real pout
       contains
       procedure,public::set=>set_Assmboundary
       !procedure,public::init=>init_Assmboundary !���ó��ڵı߽�����
    end type Assmboundary
    
    type,public::AssmMaterial!�����Ժ�ˮ��ѧ����
        real,allocatable::rho(:,:)!������
        real,allocatable::shc(:,:)
        real,allocatable::ctc(:,:)
        real,allocatable::dvs(:,:)
        real,allocatable::htc(:)
     contains
       procedure,public::init=>init_AssmMaterial
    end type AssmMaterial
     
    type,public::AssmThermal!                
        real,allocatable::Temperature(:,:) !pvt
        real,allocatable::pressure(:)
        real,allocatable::velocity(:)
      contains
       procedure,public::init=>init_AssmThermal
    end type AssmThermal
     private::set_assmgeom
     private::set_assmmesh
     private::set_Assmboundary
     private::init_AssmMaterial
     private::init_AssmThermal
    contains
     subroutine set_assmgeom(this,rFuel,GasGap,ShellThick,AssmShellThick,AcrossFlat,Height,n_pin)
       implicit none
       class(assmgeom),intent(in out)::this
       real,intent(in)::rFuel          !Ԫ���뾶
       real,intent(in)::GasGap         !Ԫ����϶���
       real,intent(in)::ShellThick     !Ԫ����Ǻ��
       real,intent(in)::AssmShellThick !�����Ǻ��
       real,intent(in)::AcrossFlat     !�����Ա߾ࣨ�������Ǻ�ȣ�
       real,intent(in)::Height         !����߶ȣ���������
       integer,intent(in)::n_pin       !Ԫ���ĸ���
       this%rFuel=rFuel
       this%GasGap=GasGap
       this%ShellThick=ShellThick
       this%AssmShellThick=AssmShellThick
       this%AcrossFlat=AcrossFlat
       this%Height=Height
       this%n_pin=n_pin
     end subroutine set_assmgeom
     
     subroutine set_assmmesh(this,nf,ng,ns,ny)
        implicit none
        class(assmmesh),intent(in out)::this
        integer,intent(in)::nf
        integer,intent(in)::ng
        integer,intent(in)::ns
        integer,intent(in)::ny
        this%nf=nf
        this%ng=ng
        this%ns=ns
        this%ny=ny     
     end subroutine set_assmmesh
     
     subroutine set_assmboundary(this,Tin,uin,pin)
       implicit none
       class(assmboundary),intent(in out)::this
       real,intent(in)::Tin
       real,intent(in)::uin
       real,intent(in)::pin
       this%Tin=Tin
       this%uin=uin
       this%pin=pin
     end subroutine set_assmboundary
     
     !subroutine init_AssmMaterial(this,LBE,he,T91)  
     subroutine init_AssmMaterial(this,Nf,Ng,Ns,Ny)
       implicit none
       class(AssmMaterial),intent(in out)::this
       integer,intent(in)::Nf
       integer,intent(in)::Ng
       integer,intent(in)::Ns
       integer,intent(in)::Ny
       integer M,N!local
       integer i,j!local
       M=Ny+1
       N=Nf+Ng+Ns+1
          Do i=0,M,1
             Do j=0,N,1
                 if (j>=0.and.j<=Nf)then !о�������� UO2
                  !RHOI(i,j)=10980
                  this%RHO(i,j)=10980
                  this%SHC(i,j)=300.0
                  this%CTC(i,j)=4.33
                  this%DVS(i,j)=0.0
                 elseif (j>Nf.and.j<=Nf+Ng) then!��϶������ He
                  !RHOI(i,j)=1.785
                  this%RHO(i,j)=1.785
                  this%SHC(i,j)=1.260
                  this%CTC(i,j)=0.124
                  this%DVS(i,j)=0.0
                 elseif (j>Nf+Ng.and.j<=Nf+Ng+Ns)then!���������� Ti
                  !RHOI(i,j)=7900
                  this%RHO(i,j)=7900
                  this%SHC(i,j)=502.42
                  this%CTC(i,j)=18.84
                  this%DVS(i,j)=0.0
                 else!��������
                  !RHOI(i,j)=10470  
                  this%RHO(i,j)=10470
                  this%SHC(i,j)=159
                  this%CTC(i,j)=3.61
                  this%DVS(i,j)=5.0 !����ճ�ȣ������˶�ճ��
                 endif                 
             enddo
         enddo    
     end subroutine init_AssmMaterial
     
     subroutine init_AssmThermal(this,Temperature,Pressure,Velocity)
        implicit none
        class(AssmThermal),intent(in out)::this
        real,intent(in)::Temperature
        real,intent(in)::Pressure
        real,intent(in)::Velocity
        this%Temperature=Temperature
        this%Pressure=Pressure
        this%Velocity=Velocity
     end subroutine init_AssmThermal
end module sys_assm_header