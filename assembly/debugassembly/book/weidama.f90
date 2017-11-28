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

subroutine solve_temperature(assm,Ti,rhoi,dt)
 implicit none
 class(sys_assembly),intent(in out)::assm
 real,intent(in)::Ti(:,:)
 real,intent(in)::rhoi(:,:)
 real dt
 !local
 call cal_th_convection(assm)
 call cal_th_temperature(assm,Ti,rhoi,dt)
 
end subroutine solve_temperature
    
    
    
    
    
    
    
subroutine cal_th_temperature(assm,Ti,rhoi,dt)
 implicit none
 class(sys_assembly),intent(in out)::assm
 real,intent(in)::Ti(:,:),rhoi(:,:)
 real,intent(in)::dt
 !local
    real Area,Xt,xf,xg,xs,Df,Dg,Ds,Dy,uin,Tin
    integer  M,N,i,j,k,Nf,Ng,Ns,Ny
    real,allocatable::Tj(:,:)
    real,allocatable::RHO(:,:),SHC(:,:),CTC(:,:),DVS(:,:)
    real,allocatable::aw(:,:),ae(:,:),ap(:,:),as(:,:),an(:,:),api(:,:),bs(:,:),q(:,:)
     Area=this%hydrau%aflow
     uin=this%boundary%v%inlet
     Tin=this%boundary%T%inlet
     xf=this%geom%rFuel
     xg=this%geom%GasGap
     xs=this%geom%ShellThick
     Nf=this%mesh%Nf
     Ng=this%mesh%Ng
     Ns=this%mesh%Ns
     M=this%mesh%Ny+1
     N=this%mesh%Nf+this%mesh%Ng+this%mesh%Ns+1
     Df=Xf/Nf
     Dg=Xg/Ng
     Ds=Xs/Ns
     Dy=Length/Ny
     
    allocate(Tj(1:M-1,1:N))
    allocate(RHO(0:M,0:N),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N))
    allocate(aw(1:M-1,1:N),ae(1:M-1,1:N),ap(1:M-1,1:N),as(1:M-1,1:N),an(1:M-1,1:N),api(1:M-1,1:N),bs(1:M-1,1:N),q(1:M-1,1:N))
    rho=this%property%rho
    shc=this%property%shc
    ctc=this%property%ctc
    dvs=this%property%dvs
    
    Do i=1,M-1,1
        Do j=1,N,1
         if (j==1)then!轴对称边界的控制体
          aw(i,j)=0.0
          ae(i,j)=(this%mesh%r(i,j)+Df/2.0)*CTC(i,j)/Df
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Df/dt
          ap(i,j)=ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Df*q(i,j)
         elseif (j>1.and.j<Nf)then!fuel内部控制体
          aw(i,j)=(this%mesh%r(i,j)-Df/2.0)*CTC(i,j)/Df
          ae(i,j)=(this%mesh%r(i,j)+Df/2.0)*CTC(i,j)/Df
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Df/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Df*q(i,j)
         elseif(j==Nf)then!f-g边界左侧控制体
          aw(i,j)=(this%mesh%r(i,j)-Df/2.0)*CTC(i,j)/Df
          ae(i,j)=2*(this%mesh%r(i,j)+Df/2.0)*(Df+Dg)/(Df/CTC(i,j)+Dg/CTC(i,j+1))/(Df+Dg)
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Df/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Df*q(i,j)  
         elseif (j==Nf+1)then!f-g边界右侧控制体
          aw(i,j)=2*(this%mesh%r(i,j)-Dg/2.0)*(df+dg)/(df/CTC(i,j-1)+dg/CTC(i,j))/(df+dg)
          ae(i,j)=(this%mesh%r(i,j)+Dg/2.0)*CTC(i,j)/Dg
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Dg*q(i,j)
         elseif(j>Nf+1.and.j<Nf+Ng)then!g气隙内部控制体
          aw(i,j)=(this%mesh%r(i,j)-Dg/2.0)*CTC(i,j)/Dg
          ae(i,j)=(this%mesh%r(i,j)+Dg/2.0)*CTC(i,j)/Dg
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Dg*q(i,j)
         elseif(j==Nf+Ng)then!g-c边界左侧控制体
          aw(i,j)=(this%mesh%r(i,j)-Dg/2.0)*CTC(i,j)/Dg
          ae(i,j)=2*(this%mesh%r(i,j)+Dg/2.0)*(Dg+Ds)/(Dg/CTC(i,j)+Ds/CTC(i,j+1))/(Dg+Ds)
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Dg*q(i,j)
         elseif(j==Nf+Ng+1)then!g-c边界右侧控制体
          aw(i,j)=2*(this%mesh%r(i,j)-Ds/2.0)*(dg+ds)/(dg/CTC(i,j-1)+ds/CTC(i,j))/(dg+ds)
          ae(i,j)=(this%mesh%r(i,j)+Ds/2.0)*CTC(i,j)/Ds
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Ds*q(i,j)
         elseif(j>Nf+Ng+1.and.j<Nf+Ng+Ns)then!c包壳内部控制体
          aw(i,j)=(this%mesh%r(i,j)-Ds/2.0)*CTC(i,j)/Ds
          ae(i,j)=(this%mesh%r(i,j)+Ds/2.0)*CTC(i,j)/Ds
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Ds*q(i,j)
         elseif(j==Nf+Ng+Ns)then!s-fluid边界左侧控制体
          aw(i,j)=(this%mesh%r(i,j)-Ds/2.0)*CTC(i,j)/Ds
          ae(i,j)=(this%mesh%r(i,j)+Ds/2.0)/(1.0/htc(i)+ds/(2.0*CTC(i,j)))
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*this%mesh%r(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=this%mesh%r(i,j)*Ds*q(i,j)
         elseif(j==Nf+Ng+Ns+1)then!fluid控制体
          Xt=Xf+Xg+Xs !包壳外径
          if(i==1)then!流体入口的控制体
           aw(i,j)=Dy/SHC(i,j)*2.0*3.14*Xt/Area*1.0/(1.0/HTC(i)+Ds/(2*CTC(i,j-1)))
           ae(i,j)=0.0
           as(i,j)=0.0
           an(i,j)=0.0
           api(i,j)=RHOI(i,j)*Dy/dt
           ap(i,j)=aw(i,j)+api(i,j)+RHO(i-1,j)*uin
           bs(i,j)=Dy/SHC(i,j)*q(i,j)+RHO(i-1,j)*uin*T(i-1,j)          
          else!流体内部以及出口控制体
           aw(i,j)=Dy/SHC(i,j)*2.0*3.14*Xt/Area*1.0/(1.0/HTC(i)+Ds/(2*CTC(i,j-1)))
           ae(i,j)=0.0
           as(i,j)=0.0
           an(i,j)=0.5*(RHO(i,j)+RHO(i-1,j))*u(i-1)
           api(i,j)=RHOI(i,j)*Dy/dt
           ap(i,j)=an(i,j)+aw(i,j)+api(i,j)
           bs(i,j)=Dy/SHC(i,j)*q(i,j)        
          endif
         endif              
      enddo
    enddo

     do k=1,100,1
       do i=1,M-1,1
           do j=1,N,1
               if(k==1) then
                Tj(i,j)=Ti(i,j)
               else
                if(j==1)then
                 Tj(i,j)=(ae(i,j)*Tj(i,j+1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                elseif(j>1.and.j<N)then
                 Tj(i,j)=(aw(i,j)*Tj(i,j-1)+ae(i,j)*Tj(i,j+1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                elseif(j==N) then
                  if(i==1)then
                    Tj(i,j)=(aw(i,j)*Tj(i,j-1)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)
                  else
                    Tj(i,j)=(aw(i,j)*Tj(i,j-1)+an(i,j)*Tj(i-1,j)+bs(i,j)+api(i,j)*Ti(i,j))/ap(i,j)  
                 endif
                endif
              endif
           enddo
       enddo
     enddo
     
     do i=0,M,1
         do j=0,N,1
            if(i>0.and.i<M)then
                if(j==0)then
                  this%thermal%Temperature(i,j)=Tj(i,j+1)
                else
                  this%thermal%Temperature(i,j)=Tj(i,j)
                endif
            elseif(i==0)then
                if(j==0)then
                  this%thermal%Temperature(i,j)=Tj(i+1,j+1)
                elseif(j>0.and.j<N)then
                  this%thermal%Temperature(i,j)=Tj(i+1,j)
                elseif(j==N)then!入口
                  this%thermal%Temperature(i,j)=Tin
                endif
           elseif(i==M)then
                if(j==0)then
                  this%thermal%Temperature(i,j)=Tj(i-1,j+1)
                else
                  this%thermal%Temperature(i,j)=Tj(i-1,j)
                endif
            endif           
         enddo
     enddo
 
end subroutine cal_th_temperature