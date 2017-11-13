module testTransient
    use sys_re_input_global
    use assm_global!assm1
    implicit none
    private
    public::cal_Assembly_Transient 
    public::cal_Assembly_Steady
 contains
    !单独调试时用assm%，等到耦合时用this%
    !subroutine函数只负责计算，输入卡输入，边界条件更新，数据初始化（物性、热工参数）由另外的函数负责
    subroutine cal_Assembly_Transient(assm1,power, fq_core, tidx, ltime, ctime)!计算一个时间步长
      implicit none
      !class(sys_assembly)::assm1
      real(KREAL), intent(in)  :: power  !(nth%na, nth%nr)
      real(KREAL), intent(in)  :: fq_core!(nth%na, nth%nr)                    ! power peak from core calculation
      integer, intent(in)      :: tidx
      real(KREAL), intent(in)  :: ltime
      real(KREAL), intent(in)  :: ctime!must
      !local
      !参数输入
      call input(xf,xg,xs,length,pd,f,Tin,pout,Tic,uic,tmax,sigma,sigmab,alpha,nf,ng,ns,ny,nt)
      !中间参数计算
      call predata(xf,xg,xs,length,pd,Tin,pout,tmax,nf,ng,ns,ny,nt,df,dg,ds,dy,rc,dt,M,N,area,wet,de)
      !数组空间分配
        allocate(XX(0:M,0:N),YY(0:M,0:N))
        allocate(T(0:M,0:N),Ti(0:M,0:N))
        allocate(RHO(0:M,0:N),RHOI(0:M,0:N),RHOF(0:M),RHOFI(0:M),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N),htc(0:M))
        allocate(u(1:Ny-1),ui(1:Ny-1),ulast(1:Ny-1))
        allocate(ap(1:Ny-1),pm(1:Ny),p(1:Ny),q(1:M-1,1:N),b(1:Ny))
        allocate(Post_T_fuel_outlet(1:Nt),Post_T_fluid_outlet(1:Nt),Post_rho_fluid_outlet(1:Nt),Tfg(0:M),Tgs(0:M),Tsf(0:M),uint(0:Nt))!后处理相关参数,后续可更新设置
       !初始化
       call initial(Nf,Ng,Ns,Ny,M,N,Nt,RHO,RHOI,RHOF,RHOFI,SHC,CTC,DVS,Ti,Tic,Tin,ui,uic,uin,uint,ulast,dt,p,pin,pout,q,T)!热物性赋值，边界条件赋值，Tmin时刻参数赋值,收敛因子赋初值,中间参数赋值(ulast)
       !求坐标
       call grid(Xf,Xg,Xs,Length,Nf,Ng,Ns,Ny,M,N,XX,YY)
       !计算
       call single_channal_transient(f,De,Area,wet,xf,xg,xs,length,dy,nf,ng,ns,ny,M,N,RHO,RHOI,RHOF,RHOFI,SHC,CTC,DVS,htc,dt,Tin,Ti,uin,p,pout,ui,ulast,u,uout,q,T,sigma,sigmab,alpha)
    end subroutine cal_Assembly_Transient
end module testTransient