module sys_assembly_header
    use sys_assm_header
    use sys_re_input_global
    
    implicit none
    type,public::sys_assembly!描述一个组件的特征和使用方法
      !private
      !real::fric  !摩擦因子
      type(hydraulic)::hydrau
      type(AssmGeom)::geom !assm_geom
      type(AssmMesh)::mesh !Assm_mesh
      type(material)::property !Assm_material 热物性和水力学参数
      type(th_boundary)::th_boundary !Assm_th_boundary
      type(AssmInit)::initdata
      type(confactor)::confactor_
      type(assmpow)::pow
      !real,allocatable::power(:,:) !Assm_power(zone,layer)
      !real,allocatable::fq_core(:,:)
      type(thermal)::Thermal  !pvt
    contains
      procedure,public::alloc=>alloc_assembly
      procedure,public::clean=>free_assembly
      procedure,public::set=>set_assembly!输入卡输入
      procedure,public::init=>init_assembly
      procedure,public::grid=>cal_grid
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
      !热物性初始化
      call this%property%init(this%mesh%Nf,this%mesh%Ng,this%mesh%Ns,this%mesh%Ny)
      !热工参数初始化
      call this%Thermal%init(this%initdata%Ti,this%initdata%Pi,this%initdata%ui)
      !边界条件初始化
      call this%th_boundary%init(this%initdata%Tin,this%initdata%uin,this%initdata%pin) 
      !热源初始化
      this%pow%power=0.0
      this%pow%fq_core=0.0
      !网格
      call this%grid()
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
      
      allocate(this%property%rho(0:M,0:N))!(1:M,1:N)
      allocate(this%property%shc(0:M,0:N))
      allocate(this%property%ctc(0:M,0:N))
	  allocate(this%property%dvs(0:M,0:N))
      allocate(this%property%htc(0:M))
      
      allocate(this%thermal%Temperature(M,N))
      allocate(this%thermal%Pressure(M))
      allocate(this%thermal%Velocity(Ny-1))
      
      allocate(this%mesh%r(0:M,0:N))
      allocate(this%mesh%z(0:M,0:N))
      
      allocate(this%pow%power(M))
      allocate(this%pow%fq_core(M))
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
      
      if(allocated(this%mesh%r))  deallocate(this%mesh%r)
      if(allocated(this%mesh%z))  deallocate(this%mesh%z)
      
      if(allocated(this%pow%power))  deallocate(this%pow%power)
      if(allocated(this%pow%fq_core))  deallocate(this%pow%fq_core)
     end subroutine Free_assembly
     
     subroutine set_assembly(this,reInputdata)
      implicit none
      class(sys_assembly),intent(in out)::this
      type(sys_re_input),intent(in)::reInputdata
      !设置几何参数
      call this%geom%set(reInputdata%xf,reInputdata%xg,reInputdata%xs,reInputdata%xos,reInputdata%acf,reInputdata%Height,reInputdata%npin)
      !设置网格参数
      call this%mesh%set(reInputdata%nf,reInputdata%ng,reInputdata%ns,reInputdata%ny)
      !设置初始值
      call this%initdata%set(reInputdata%Ti,reInputdata%Pi,reInputdata%Ui,reInputdata%Tin,reInputdata%Pin,reInputdata%Uin)
      !设置收敛因子
      call this%confactor_%set(reInputdata%alpha,reInputdata%sigma)
      call this%hydrau%set(reInputdata%f)
     end subroutine set_assembly
     
     subroutine cal_grid(this)
       implicit none
       class(sys_assembly),intent(in out)::this
       !local
       real Df,Dg,Ds,Dy 
       integer  M,N,i,j
     
     Df=this%geom%rFuel/this%mesh%Nf
     Dg=this%geom%GasGap/this%mesh%Ng
     Ds=this%geom%ShellThick/this%mesh%Ns
     Dy=this%geom%Height/this%mesh%Ny
     M=this%mesh%Ny+1
     N=this%mesh%Nf+this%mesh%Ng+this%mesh%Ns+1
     
     Do i=0,M,1
         do j=0,N,1
            if (j==0)then
               this%mesh%r(i,j)=0.0
            elseif(j==1)then
               this%mesh%r(i,j)=Df/2.0
            elseif(j>1.and.j<=this%mesh%Nf) then
               this%mesh%r(i,j)=this%mesh%r(i,j-1)+Df
            elseif(j==this%mesh%Nf+1)then
               this%mesh%r(i,j)=this%mesh%r(i,j-1)+Df/2.0+Dg/2.0
            elseif (j>this%mesh%Nf+1.and.j<=this%mesh%Nf+this%mesh%Ng) then
               this%mesh%r(i,j)=this%mesh%r(i,j-1)+Dg
            elseif (j==this%mesh%Nf+this%mesh%Ng+1)then
               this%mesh%r(i,j)=this%mesh%r(i,j-1)+ Dg/2.0+Ds/2.0
            elseif (j>this%mesh%Nf+this%mesh%Ng+1.and.j<=this%mesh%Nf+this%mesh%Ng+this%mesh%Ns)then
               this%mesh%r(i,j)=this%mesh%r(i,j-1)+Ds
            else!流体的径向坐标，没有实际意义
               this%mesh%r(i,j)=this%mesh%r(i,j-1)+Ds
            endif
            
            if(i==0)then
              this%mesh%z(i,j)=0.0
            elseif (i==1)then
              this%mesh%z(i,j)=Dy/2.0
            elseif (i>1.and.i<M)then
              this%mesh%z(i,j)=this%mesh%z(i-1,j)+Dy
            elseif (i==M)then
              this%mesh%z(i,j)=this%mesh%z(i-1,j)+Dy/2.0
            endif
         enddo
      enddo   
     end subroutine cal_grid
     
end module sys_assembly_header
