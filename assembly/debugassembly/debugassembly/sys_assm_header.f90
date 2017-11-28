module sys_assm_header
    use mathkerel
    use sys_property
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
        real,allocatable::r(:,:)
        real,allocatable::z(:,:)
      contains
      procedure,public::set=>set_assmmesh
      end type Assmmesh
    
    type,public::boundary
        real inlet
        real outlet
    end type boundary
    
    type,public::iteration
        real,allocatable::Temperature(:,:) !pvt
        real,allocatable::pressure(:)
        real,allocatable::velocity(:)
    end type iteration
      
    type,public::th_boundary
       type(boundary)::p
       type(boundary)::u
       type(boundary)::T
       !type(boundary)::rho
       contains
       procedure,public::init=>init_th_boundary
       !procedure,public::init=>init_th_boundary !���ó��ڵı߽�����
    end type th_boundary
    
    type,public::hydraulic
        real fric
        real aflow
        real wet
        real de
    contains
        procedure,public::set=>set_hydraulic
        procedure,public::cal=>cal_hydraulic
    end type hydraulic
    
        
    type,public::material!�����Ժ�ˮ��ѧ����
        real,allocatable::rho(:,:)!������
        real,allocatable::shc(:,:)
        real,allocatable::ctc(:,:)
        real,allocatable::dvs(:,:)
        real,allocatable::htc(:)
     contains
       procedure,public::init=>init_material
    end type material
     
    type,public::thermal!                
        real,allocatable::Temperature(:,:) !pvt
        real,allocatable::pressure(:)
        real,allocatable::velocity(:)
      contains
       procedure,public::init=>init_thermal
      end type thermal
      
    type,public::AssmInit
      real Ti!��ʼ�¶�
      real Pi!��ʼѹ��
      real Ui!��ʼ�ٶ�
      real Tin
      real Pin
      real Uin
    contains
      procedure,public::set=>set_assminit
    end type AssmInit
    
    type,public::confactor
      real alpha!��ʼ�¶�
      real sigma!��ʼѹ��
    contains
      procedure,public::set=>set_confactor
    end type confactor
    
    type,public::assmpow
        real,allocatable::power(:)
        real,allocatable::fq_core(:)
    end type
        
    type,public::sys_time
        type(th_boundary)::th_boundary
        type(thermal)::thermal
        type(material)::material
        type(assmpow)::power
    end type sys_time
    
     private::set_assmgeom
     private::set_assmmesh
     private::set_assminit
     private::set_confactor
     private::set_hydraulic
     private::init_th_boundary!����ʱ��仯������init
     private::init_material
     private::init_thermal
     private::cal_hydraulic
     !private::cal_grid
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
     
     subroutine init_th_boundary(this,Tin,uin,pin)
       implicit none
       class(th_boundary),intent(in out)::this
       real,intent(in)::Tin
       real,intent(in)::uin
       real,intent(in)::pin
       this%T%inlet=Tin
       this%u%inlet=uin
       this%p%inlet=pin
       !this%rho%inlet=get_density(Tin)
     end subroutine init_th_boundary
     
     !subroutine init_material(this,LBE,he,T91)  
     subroutine init_material(this,Nf,Ng,Ns,Ny)
       implicit none
       class(material),intent(in out)::this
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
     end subroutine init_material
     
     subroutine init_thermal(this,Temperature,Pressure,Velocity)
        implicit none
        class(thermal),intent(in out)::this
        real,intent(in)::Temperature
        real,intent(in)::Pressure
        real,intent(in)::Velocity
        this%Temperature=Temperature
        this%Pressure=Pressure
        this%Velocity=Velocity
     end subroutine init_thermal
     
     subroutine set_AssmInit(this,Ti,Pi,Ui,Tin,Pin,Uin)
        implicit none
        class(AssmInit),intent(in out)::this
        real,intent(in)::Ti
        real,intent(in)::Pi
        real,intent(in)::Ui
        real,intent(in)::Tin
        real,intent(in)::Pin
        real,intent(in)::Uin
        this%Ti=Ti
        this%Pi=Pi
        this%Ui=Ui
        this%Tin=Tin
        this%Pin=Pin
        this%Uin=Uin
     end subroutine set_AssmInit
     
     subroutine set_confactor(this,alpha,sigma)
        implicit none
        class(confactor),intent(in out)::this
        real,intent(in)::alpha
        real,intent(in)::sigma
        this%alpha=alpha
        this%sigma=sigma
     end subroutine set_confactor
 
     subroutine set_hydraulic(this,fric)
        implicit none
        class(hydraulic),intent(in out)::this
        real,intent(in)::fric
        this%fric=fric
     end subroutine set_hydraulic
     
     subroutine cal_hydraulic(this,rc,pd)
        implicit none
        class(hydraulic),intent(in out)::this
        real,intent(in)::rc,pd
        call get_hyconstant(rc,pd,this%Aflow,this%wet,this%de)
     end subroutine cal_hydraulic
     
end module sys_assm_header