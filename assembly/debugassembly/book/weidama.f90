!����һ��ʱ�䲽��
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

  call solve_momentum(this,rhofi,ui,dt)!��Ҫ���뵱ǰ��������rho,uz 
  call solve_pressureCorrection(this,ap,rhofi,dt,pmodify)
  call modify_PV(this,pmodify)
  !����btotal
    btotal=0.0
    do i=1,Ny,1
      btotal=btotal+abs(b(i))
    enddo
  end do
  call solve_temperature(this,Ti,rhoi,dt)
  call ����drho
 end do
  call thermal set
 enddo

!********************
 !solve_momentum   ����ѹ�������һά˲̬�������̣����һ��ʱ�䲽��
 !����:���γߴ硢����N���߽�����uin/pin�����ԣ���ʼPV��ʱ�䲽��dt
 !���:һ��ʱ�䲽���Ժ��PVT
!********************
subroutine cal_Assembly_Transient(this,dt)
 class(sys_assembly),intent(in out)::this
 real,intent(in)::dt
 !local rhoi/ui/Ti/dt/ap/pmodify/
 integer M,N,Ny,i
 real btotal,drho!�ж�����
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
   call solve_momentum(this,rhofi,ui,dt,ap)!��Ҫ���뵱ǰ��������rho,uz 
   call solve_pressureCorrection(this,ap,rhofi,dt,pmodify,btotal)
   call modify_PV(this,ap,pmodify)
  end do
  call solve_temperature(this,Ti,rhoi,dt)
  call update_property(this,drho)!���Ը���
 end do
end subroutine cal_Assembly_Transient
    
    
    
