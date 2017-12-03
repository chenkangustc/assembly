!计算一个时间步长
subroutine cal_Assembly_Transient(this,dt)
 class(sys_assembly),intent(in out)::this
 real,intent(in)::dt
 !real,intent(in)::ltime
 !real,intent(out)::ctime
 !local
 integer Ny,i
 real,allocatable::pmodify(:)
 !real,allocatable::pguess(:)
 real,allocatable::ui(:),Ti(:,:)
 real,allocatable::rhoi(:,:),rho(:,:),rhofi(:)
 real,allocatable::ap(:)
 real btotal,drho
 Ny=this%mesh%Ny
 allocate(pmodify,rhoi,rhol,ui,ap,Ti(:,:))
 
 ui=this%thermal%u
 rhoi=this%property%density
 rhol=this%property%density
 rho=this%property%density
 
 drho=1.0!
 do while(drho>sigma)
     btotal=1.0
  do while(btotal>sigmab)
  !pguess(:)=this%thermal%p(:)

  call solve_momentum(this,rhofi,ui,dt)!需要输入当前迭代步的rho,uz 
  call solve_pressureCorrection(this,ap,rhofi,dt,pmodify)
  call modify_PV(this,pmodify)
  !计算btotal
    btotal=0.0
    do i=1,Ny,1
      btotal=btotal+abs(b(i))
    enddo
  end do
  call solve_temperature(this,Ti,rhoi,dt)
  call 计算drho
 end do
  call thermal set
 enddo

!********************
 !solve_momentum   假设压力，求解一维瞬态动量方程，求解一个时间步长
 !输入:几何尺寸、网格N、边界条件uin/pin，物性，初始PV，时间步长dt
 !输出:一个时间步长以后的PVT
!********************
subroutine cal_Assembly_Transient(this,dt)
 class(sys_assembly),intent(in out)::this
 real,intent(in)::dt
 !local rhoi/ui/Ti/dt/ap/pmodify/
 integer M,N,Ny,i
 real btotal,drho!判断因子
 real,allocatable::pmodify(:)
 real,allocatable::ui(:),Ti(:,:)
 real,allocatable::rhoi(:,:),rhofi(:)
 real,allocatable::ap(:)
 
 Ny=this%mesh%Ny
 M=Ny+1
 N=this%mesh%Nf+this%mesh%Ng+this%mesh%Ns+1
 allocate(pmodify(Ny),ui(Ny-1),Ti(M-1,N),rhoi(0:M,0:N),rhofi(0:M),ap(Ny-1))
 
 pmodify=0.0
 ui=this%thermal%velocity
 rhoi=this%property%density
 rhofi=this%property%density(:,N)
 ap=0.0
 drho=1.0
 do while(drho>assm%confactor_%alpha)
     btotal=1.0
  do while(btotal>assm%confactor_%sigma)
   call solve_momentum(this,rhofi,ui,dt,ap)!需要输入当前迭代步的rho,uz 
   call solve_pressureCorrection(this,ap,rhofi,dt,pmodify,btotal)
   call modify_PV(this,ap,pmodify)
  end do
  call solve_temperature(this,Ti,rhoi,dt)
  call update_property(this,drho)!物性更新
 end do
end subroutine cal_Assembly_Transient
    
    
    
