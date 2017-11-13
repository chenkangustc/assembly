module sys_assm_header
    implicit none
    type,public::AssmGeom
      real rFuel          !元件半径
      real GasGap         !元件气隙厚度
      real ShellThick     !元件外壳厚度
      real AssmShellThick !组件外壳厚度
      real AcrossFlat     !组件外对边距（包含包壳厚度）
      real Height         !组件高度（活性区）
      integer n_pin       !元件的个数
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
       !procedure,public::init=>init_Assmboundary !设置出口的边界条件
    end type Assmboundary
    
    type,public::AssmMaterial!热物性和水力学参数
        real,allocatable::rho(:,:)!热物性
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
       real,intent(in)::rFuel          !元件半径
       real,intent(in)::GasGap         !元件气隙厚度
       real,intent(in)::ShellThick     !元件外壳厚度
       real,intent(in)::AssmShellThick !组件外壳厚度
       real,intent(in)::AcrossFlat     !组件外对边距（包含包壳厚度）
       real,intent(in)::Height         !组件高度（活性区）
       integer,intent(in)::n_pin       !元件的个数
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
                 if (j>=0.and.j<=Nf)then !芯块热物性 UO2
                  !RHOI(i,j)=10980
                  this%RHO(i,j)=10980
                  this%SHC(i,j)=300.0
                  this%CTC(i,j)=4.33
                  this%DVS(i,j)=0.0
                 elseif (j>Nf.and.j<=Nf+Ng) then!气隙热物性 He
                  !RHOI(i,j)=1.785
                  this%RHO(i,j)=1.785
                  this%SHC(i,j)=1.260
                  this%CTC(i,j)=0.124
                  this%DVS(i,j)=0.0
                 elseif (j>Nf+Ng.and.j<=Nf+Ng+Ns)then!包壳热物性 Ti
                  !RHOI(i,j)=7900
                  this%RHO(i,j)=7900
                  this%SHC(i,j)=502.42
                  this%CTC(i,j)=18.84
                  this%DVS(i,j)=0.0
                 else!流体物性
                  !RHOI(i,j)=10470  
                  this%RHO(i,j)=10470
                  this%SHC(i,j)=159
                  this%CTC(i,j)=3.61
                  this%DVS(i,j)=5.0 !动力粘度，而非运动粘度
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