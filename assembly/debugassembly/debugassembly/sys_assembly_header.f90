module sys_assembly_header
    use sys_assm_header
    use sys_re_input_global
    implicit none
    type,public::sys_assembly!描述一个组件的特征和使用方法
      private
      real::fric  !摩擦因子
      type(AssmGeom)::geom !assm_geom
      type(AssmMesh)::mesh !Assm_mesh
      type(AssmMaterial)::property !Assm_material 热物性和水力学参数
      type(AssmBoundary)::boundary !Assm_boundary
      real::power(:,:) !Assm_power(zone,layer)
      type(AssmThermal)::Thermal  !pvt
    contains
      procedure,public::alloc=>alloc_assembly
      procedure,public::clean=>free_assembly
      procedure,public::set=>set_assembly!输入卡输入
      procedure,public::init=>init_assembly
      !procedure,public::Steady=>cal_Assembly_Steady
      !procedure,public::Transient=>cal_Assembly_Transient
    end type sys_assembly
     private::alloc_assembly
     private::clean_assembly
     private::set_assembly
     private::init_assembly
!     private::cal_Assembly_Steady
!     private::cal_Assembly_Transient
    contains
     subroutine init_assembly(this)
      implicit none
      class(sys_assembly),intent(in out)::this
      !热物性初始化
      call this%property%init(this%mesh%Nf,this%mesh%Ng,this%mesh%Ns,this%mesh%Ny)
      !热工参数初始化
      call this%Thermal%init(reInputdata%Ti,reInputdata%Pi,reInputdata%ui)
     endsubroutine init_assembly
    
     subroutine alloc_assembly(this)
      implicit none
      class(sys_assembly),intent(in out)::this
      !local
      integer::i_allocate
      integer::M,N
      M=this%mesh%Ny+1
      N=this%mesh%Nf+this%mesh%Ng+this%mesh%Ns+1
      !check allocated first
      call this%clean()
      
      allocate(this%property%rho(0:M,0:N))
      allocate(this%property%shc(0:M,0:N))
      allocate(this%property%ctc(0:M,0:N))
      allocate(this%property%htc(0:M))
      
      allocate(this%thermal%Temperature(0:M,0:N))
      allocate(this%thermal%Pressure(1:this%mesh%Ny))
      allocate(this%thermal%Velocity(1:this%mesh%Ny-1))
     end subroutine alloc_assembly
     
     subroutine Free_assembly(this)
      implicit none
      class(sys_assembly),intent(in out)::this
      if(allocated(this%property%rho))  deallocate(this%property%rho)
      if(allocated(this%property%shc))  deallocate(this%property%shc)
      if(allocated(this%property%ctc))  deallocate(this%property%ctc)
      if(allocated(this%property%htc))  deallocate(this%property%htc)
      
      if(allocated(this%thermal%temperature))  deallocate(this%thermal%temperature)
      if(allocated(this%thermal%pressure))  deallocate(this%thermal%pressure)
      if(allocated(this%thermal%Velocity))  deallocate(this%thermal%Velocity)
     end subroutine Free_assembly
     
     subroutine set_assembly(this)
      implicit none
      class(sys_assembly),intent(in out)::this
      call this%geom%set(reInputdata%xf,reInputdata%xg,reInputdata%xs,reInputdata%xos,reInputdata%acf,reInputdata%Height,reInputdata%npin)
      call this%mesh%set(reInputdata%nf,reInputdata%ng,reInputdata%ns,reInputdata%ny)
      this%fric=reInputdata%f
     end subroutine set_assembly
end module sys_assembly_header