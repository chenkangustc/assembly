module driving_cal_transient
    use assm_global
    use mathkerel
    implicit none
    private
    public::cal_Assembly_transient
    public::solve_momentum
contains
    subroutine cal_Assembly_transient(assm)
      implicit none
      !class(sys_assembly),intent(in out)::this 
      type(sys_assembly),intent(in out)::assm
      !real,intent(in)::power(:)
      !real,intent(in)::fq_core(:)
      integer Ny
      !Ny=SIZE(power)

    end subroutine cal_Assembly_transient
    
    subroutine solve_momentum(Ny,f,De,dx,dt,p,pout,uin,ui,ulast,RHO,RHOI,ap,u,uout)!app为aP数组，主要输出为压力修正方程所用。
      integer N,Ny,i
      real f,De,dx,dt,uin,uout,pout,pin
      !real,allocatable::p(:),ui(:),ulast(:),u(:),RHO(:),RHOI(:)
      real p(1:Ny),ui(1:Ny-1),ulast(1:Ny-1),u(1:Ny-1),RHO(0:Ny+1),RHOI(0:Ny+1),ap(1:Ny-1)
      real,allocatable::aw(:),ae(:),api(:),bs(:),A(:,:),b(:),uc(:)
     
      N=Ny-1!矢量控制体个数
      !allocate(p(0:Ny+1),ui(0:Ny),ulast(0:Ny),u(0:Ny),RHO(0:Ny+1),RHOI(0:Ny+1))
      allocate(aw(1:N),ae(1:N),api(1:N),bs(1:N),A(1:N,1:N),b(1:N),uc(1:N))
      pin=(3*p(1)-p(2))/2.0
     !计算各个控制体的常系数和源项
     do i=1,N,1
         if(i==1)then
             aw(i)=0.0
             ae(i)=0.0
             api(i)=1.50*dx/dt*(2*RHOI(i)+RHOI(i+1))/3.0
             ap(i)=api(i)+f/(2*De)*1.50*dx*ulast(i)*(2*RHO(i)+RHO(i+1))/3.0+RHO(0)*uin
             bs(i)=pin-p(i+1)-(2*RHO(i)+RHO(i+1))/3.0*9.8*1.5*dx+RHO(0)*uin*uin             
         elseif(i>1.and.i<N)then
             aw(i)=ulast(i-1)*RHO(i)
             ae(i)=0.0
             api(i)=dx/dt*(RHOI(i)+RHOI(i+1))/2.0
             ap(i)=aw(i)+api(i)+f/(2*De)*dx*ulast(i)*(RHO(i)+RHO(i+1))/2.0
             bs(i)=p(i)-p(i+1)-(RHO(i)+RHO(i+1))/2.0*9.8*dx
         elseif(i==N)then
             aw(i)=RHO(N)*ulast(N-1)
             ae(i)=0.0
             api(i)=1.50*dx/dt*(RHOI(N)+2*RHOI(N+1))/3.0
             ap(i)=aw(i)+api(i)+f/(2*De)*1.50*dx*ulast(N)*(RHO(N)+2*RHO(N+1))/3.0
             bs(i)=p(i)-pout-(RHO(i)+2*RHO(i+1))/3.0*9.8*1.5*dx            
         endif
        enddo
      !用直接发求矩阵，第一步，写出系数矩阵
      A=0.0
       do i=1,N,1
           if(i==1)then
               A(i,i)=ap(i)
           elseif(i>1.and.i<N)then
               A(i,i)=ap(i)
               A(i,i-1)=-aw(i)
           elseif(i==N)then
               A(i,i)=ap(i)
               A(i,i-1)=-aw(i)
           endif
           b(i)=bs(i)+api(i)*ui(i)
       enddo
    
      !第二步，调用TDMA方法求解
      call tdma(N,A,b,uc)
       !print *,RHO(0)*ui(0)*ui(0)
      !do i=1,N,1
      !  if(i==1)then
      !      uc(i)=b(i)/ap(i)
      !  else
      !      uc(i)=(b(i)+aw(i)*uc(i-1))/ap(i)
            !print *,aw(i)*uc(i-1)
       ! endif
      !enddo

      do i=1,N,1
        if(uc(i)>0.0)then
        u(i)=uc(i)
        else
        u(i)=0.0
        endif
      enddo
        uout=uc(N)
     ! print *,ap
    end subroutine solve_momentum
    
    
end module driving_cal_transient    