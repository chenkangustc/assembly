module sys_assembly_header
    use sys_assm_header
    use sys_re_input_global
    implicit none
    type,public::sys_assembly!����һ�������������ʹ�÷���
      !private
      real::fric  !Ħ������
      type(AssmGeom)::geom !assm_geom
      type(AssmMesh)::mesh !Assm_mesh
      type(material)::property !Assm_material �����Ժ�ˮ��ѧ����
      type(boundary)::boundary !Assm_boundary
      type(AssmInit)::initdata
      type(confactor)::confactor_
      type(assmpow)::pow
      !real,allocatable::power(:,:) !Assm_power(zone,layer)
      !real,allocatable::fq_core(:,:)
      type(thermal)::Thermal  !pvt
    contains
      procedure,public::alloc=>alloc_assembly
      procedure,public::clean=>free_assembly
      procedure,public::set=>set_assembly!���뿨����
      procedure,public::init=>init_assembly
      !procedure,public::Steady=>cal_Assembly_Steady
      !procedure,public::Transient=>cal_Assembly_Transient
    end type sys_assembly
     private::alloc_assembly
     private::free_assembly
     private::set_assembly
     private::init_assembly
!     private::cal_Assembly_Steady
!     private::cal_Assembly_Transient
    contains
     subroutine init_assembly(this)
      implicit none
      class(sys_assembly),intent(in out)::this
      !�����Գ�ʼ��
      call this%property%init(this%mesh%Nf,this%mesh%Ng,this%mesh%Ns,this%mesh%Ny)
      !�ȹ�������ʼ��
      call this%Thermal%init(this%initdata%Ti,this%initdata%Pi,this%initdata%ui)
      !�߽�������ʼ��
      call this%boundary%init(this%initdata%Tin,this%initdata%uin,this%initdata%pin) 
      !��Դ��ʼ��
      this%pow%power=0.0
      this%pow%fq_core=0.0
     endsubroutine init_assembly
    
     subroutine alloc_assembly(this)
      implicit none
      class(sys_assembly),intent(in out)::this
      !local
      integer::i_allocate
      integer::M,N,Ny
      Ny=this%mesh%Ny
      M=Ny+1
      N=this%mesh%Nf+this%mesh%Ng+this%mesh%Ns+1
      !check allocated first
      call this%clean()
      
      allocate(this%property%rho(0:M,0:N))
      allocate(this%property%shc(0:M,0:N))
      allocate(this%property%ctc(0:M,0:N))
	  allocate(this%property%dvs(0:M,0:N))
      allocate(this%property%htc(0:M))
      
      allocate(this%thermal%Temperature(0:M,0:N))
      allocate(this%thermal%Pressure(1:Ny))
      allocate(this%thermal%Velocity(1:Ny-1))
      
      allocate(this%pow%power(1:Ny))
      allocate(this%pow%fq_core(1:Ny))
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
      
      if(allocated(this%pow%power))  deallocate(this%pow%power)
      if(allocated(this%pow%fq_core))  deallocate(this%pow%fq_core)
     end subroutine Free_assembly
     
     subroutine set_assembly(this,reInputdata)
      implicit none
      class(sys_assembly),intent(in out)::this
      type(sys_re_input),intent(in)::reInputdata
      !���ü��β���
      call this%geom%set(reInputdata%xf,reInputdata%xg,reInputdata%xs,reInputdata%xos,reInputdata%acf,reInputdata%Height,reInputdata%npin)
      !�����������
      call this%mesh%set(reInputdata%nf,reInputdata%ng,reInputdata%ns,reInputdata%ny)
      !���ó�ʼֵ
      call this%initdata%set(reInputdata%Ti,reInputdata%Pi,reInputdata%Ui,reInputdata%Tin,reInputdata%Pin,reInputdata%Uin)
      !������������
      call this%confactor_%set(reInputdata%alpha,reInputdata%sigma)
      this%fric=reInputdata%f
     end subroutine set_assembly
     

end module sys_assembly_header
