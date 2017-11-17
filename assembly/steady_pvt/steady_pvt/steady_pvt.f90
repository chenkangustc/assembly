!****************************************************************************
!
!  PROGRAM: Console84
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Console84

    implicit none
    !��������
    real Xf,Xg,Xs,rc,Length,Area,Df,Dg,Ds,Dy,pd!���β���
    real f,De,wet!ˮ��ѧ����
    real tmax,dt,dx,sigma,sigmab    !ʱ�����
    real Tin,uin,uout,pin,pout,uini,Tic,uic  !�߽�����
    real btotal,alpha !�����о�,alpha��ѹ���������̵��ɳ�����
    real,allocatable:: Ti(:,:),ui(:),RHOI(:,:)!��ʼ����
    integer  Nf,Ng,Ns,Ny,Nt,M,N,i,k,j!�������
    real,allocatable::XX(:,:),YY(:,:)!��������
    real,allocatable::RHO(:,:),SHC(:,:),CTC(:,:),DVS(:,:),q(:,:),htc(:)!���Բ�������Դ
    real,allocatable::u(:),T(:,:),p(:)!�����
    real,allocatable::RHOF(:),RHOFI(:)
    real,allocatable::ulast(:),ap(:),pm(:),b(:)!�м����,b��ѹ���������̵�Դ������ж�P-V����
    real,allocatable::Post_T_fuel_outlet(:),Post_T_fluid_outlet(:),Post_rho_fluid_outlet(:),Tfg(:),Tgs(:),Tsf(:),uint(:)!������ز���,�����ɸ�������,��ű߽��¶�
    
    !��������
    call input(xf,xg,xs,length,pd,f,Tin,pout,Tic,uic,tmax,sigma,sigmab,alpha,nf,ng,ns,ny,nt)
    !�м��������
    call predata(xf,xg,xs,length,pd,Tin,pout,tmax,nf,ng,ns,ny,nt,df,dg,ds,dy,rc,dt,M,N,area,wet,de)
    !����ռ����
      allocate(XX(0:M,0:N),YY(0:M,0:N))
      allocate(T(0:M,0:N),Ti(0:M,0:N))
      allocate(RHO(0:M,0:N),RHOI(0:M,0:N),RHOF(0:M),RHOFI(0:M),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N),htc(0:M))
      allocate(u(1:Ny-1),ui(1:Ny-1),ulast(1:Ny-1))
      allocate(ap(1:Ny-1),pm(1:Ny),p(1:Ny),q(1:M-1,1:N),b(1:Ny))
      allocate(Post_T_fuel_outlet(1:Nt),Post_T_fluid_outlet(1:Nt),Post_rho_fluid_outlet(1:Nt),Tfg(0:M),Tgs(0:M),Tsf(0:M),uint(0:Nt))!������ز���,�����ɸ�������
    !��ʼ��
    call initial(Nf,Ng,Ns,Ny,M,N,Nt,RHO,RHOI,RHOF,RHOFI,SHC,CTC,DVS,Ti,Tic,Tin,ui,uic,uin,uint,ulast,dt,p,pin,pout,q,T)!�����Ը�ֵ���߽�������ֵ��Tminʱ�̲�����ֵ,�������Ӹ���ֵ,�м������ֵ(ulast)
    !������
    call grid(Xf,Xg,Xs,Length,Nf,Ng,Ns,Ny,M,N,XX,YY)
    !����
    call single_channal(f,De,Area,wet,xf,xg,xs,length,dy,nf,ng,ns,ny,M,N,RHO,RHOI,RHOF,RHOFI,SHC,CTC,DVS,htc,dt,Tin,Ti,uin,p,pout,ui,ulast,u,uout,q,T,sigma,sigmab,alpha)
    
    print *,'T=',Tmax,'��������¶�T=',T(M,N)
    print *,'T=',Tmax,'о����������¶�T=',maxval(T(0,:))
    print *,'T=',Tmax,'о���������¶�T=',maxval(Tfg(:))
    print *,'T=',Tmax,'�����ڱ�������¶�T=',maxval(Tgs(:))
    print *,'T=',Tmax,'�������������¶�T=',maxval(Tsf(:))
    print *,'T=',Tmax,'��������V=',uout
    print *,'T=',Tmax,'���ѹ��P=',pin
    print *,'T=',Tmax,'����ѹ��P=',pin,p,pout
    call CFD_post(Nf,Ng,Ns,Ny,M,N,Nt,T,uin,u,uout,pin,p,pout,RHO,Post_T_fluid_outlet,Post_rho_fluid_outlet,Post_T_fuel_outlet,Tfg,Tgs,Tsf)
    
    contains
    subroutine single_channal(f,De,Area,wet,xf,xg,xs,length,dy,nf,ng,ns,ny,M,N,RHO,RHOI,RHOF,RHOFI,SHC,CTC,DVS,htc,dt,Tin,Ti,uin,p,pout,ui,ulast,u,uout,q,T,sigma,sigmab,alpha)
      real f,De,Area,wet,xf,xg,xs,length,dy,dt,Tin,uin,pout,uout,sigma,sigmab,alpha
      real drho,btotal
      integer nf,ng,ns,ny,M,N
      real T(0:M,0:N),Ti(0:M,0:N) 
      real RHO(0:M,0:N),RHOI(0:M,0:N),RHOF(0:M),RHOFI(0:M),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N),htc(0:M)
      real u(1:Ny-1),ui(1:Ny-1),ulast(1:Ny-1)
      real ap(1:Ny-1),pm(1:Ny),p(1:Ny),q(1:M-1,1:N),b(1:Ny)
    
    do j=1,Nt,1!ʱ�䲽����!                                     
        drho=1.0
       do while(drho>sigma)!rho����
           btotal=1.0
           do while(btotal>sigmab)!p-v����
             !���ٶ� 
             call solve_momentum(Ny,f,De,dy,dt,p,pout,uin,ui,ulast,RHOF,RHOFI,ap,u,uout)  
             !if(1.eq.2)then
             call solve_pressureCorrection(Ny,ap,RHOF,RHOFI,u,uin,pout,dy,dt,pm,b)!�˴�����һ������õ���u
             !print *,pm
             !����ѹ�����ٶ�
             do i=1,Ny-1,1
               p(i)=p(i)+alpha*pm(i)
               if(i==1)then
                 u(i)=u(i)+1.5*(pm(i)-pm(i+1))/ap(i)
               elseif(i>1.and.i<Ny-1)then
                 u(i)=u(i)+(pm(i)-pm(i+1))/ap(i)
               elseif(i==Ny-1)then
               !u(i)=u(i)+(pm(i)-pout)/ap(i)!pmout=0.0
                 u(i)=u(i)+pm(i)/ap(i)
               endif
               ulast(i)=u(i)
            enddo
            !ѹ�����������оݼ���
            btotal=0.0
            do i=1,Ny,1
              btotal=btotal+abs(b(i))
            enddo
            !�߽����
            uout=u(Ny-1)
            pin=1.50*p(1)-0.5*p(2)
            p(Ny)=(2*pout+p(Ny-1))/3.0
        enddo!k,p-v����
    !endif
    !if(1.eq.2) then
    !�����¶ȳ�
    !print *,uin,u
      call solve_temperature(De,Area,wet,Xf,Xg,Xs,Length,Nf,Ng,Ns,Ny,M,N,RHO,RHOI,SHC,CTC,DVS,htc,dt,Tin,Ti,uin,u,q,T)
    !print *,T(:,N)
    !solve_temperature(De,Area,wet,Xf,Xg,Xs,Length,Nf,Ng,Ns,Ny,M,N,IA,RHO,RHOI,SHC,CTC,DVS,dt,Ti,u,q,T)
    !print *,T(:,N)
    !�����µ����Բ���
      do i=1,M-1,1
         RHO(i,N)=11096.0-1.3326*T(i,N)
      enddo
      RHO(0,N)=RHO(1,N)
      RHO(M,N)=RHO(M-1,N)
    
      drho=0.0
      do i=0,M,1
        drho=drho+abs((RHO(i,N)-RHOF(i))/RHOF(i))
      enddo
      !print *,drho
      do i=0,M,1
        RHOF(i)=RHO(i,N)
        RHOFI(i)=RHOI(i,N)
      enddo
     !endif
     !���
     enddo !do while(drho<sigma)
    !��ʼ�ٶȡ���ʼѹ������ʼ�¶����¸�ֵ
     ui=u
     ulast=ui
     Ti=T
    !�߽��������¸�ֵ 
     uin=uint(j)
       ! call get_uin(j,dt,uini,uin)
       ! if (uin>0.3*uini)then
       !     uin=uini/(1+j*dt/5.5)
       ! else
       !     uin=0.3*uini
       ! endif

    !����ֵ
    do i=0,M,1
        Tfg(i)=(CTC(i,Nf)*(Xg/Ng)*T(i,Nf)+CTC(i,Nf+1)*(Xf/Nf)*T(i,Nf+1))/(CTC(i,Nf)*(Xg/Ng)+CTC(i,Nf+1)*(Xf/Nf))!о����߽�
        Tgs(i)=(CTC(i,Nf+Ng)*(Xs/Ns)*T(i,Nf+Ng)+CTC(i,Nf+Ng+1)*(Xg/Ng)*T(i,Nf+Ng+1))/(CTC(i,Nf+Ng)*(Xs/Ns)+CTC(i,Nf+Ng+1)*(Xg/Ng))!�����ڱ߽�
        Tsf(i)=(htc(i)*T(i,N)+2*CTC(i,N-1)/(Xs/Ns)*T(i,N-1))/(htc(i)+2*CTC(i,N-1)/(Xs/Ns))!������߽�
    enddo
     Post_T_fluid_outlet(j)=T(M,N)!��������¶�
     !Post_rho_fluid_outlet(j)=RHO(M,N)
     Post_rho_fluid_outlet(j)=maxval(T(:,0))!������������¶�
     !Post_rho_fluid_outlet(j)=CTC(2,N-1)*2/(Xs/Ns)*(T(2,N-1)-Tsf(2))!���ڿ�����߽�����
     Post_T_fuel_outlet(j)=maxval(Tgs(:))!�����ڱ�������¶�
      !Post_T_fuel_outlet(j)=SHC(M-1,N)*RHO(M-1,N)*u(M-2)*area*(T(M,N)-T(0,N))!3.14*(xf+xg+xs)*(xf+xg+xs)*length*q(1,2)!���ߵĹ�����ʱ��仯
      !Post_T_fuel_outlet(j)=T(M,1)!��������¶�
      !print *,'t=',j*dt
      !print *,uin,u,uout
      !print *,pin,p,pout
      !print *,T(:,1)
     enddo !ʱ�䲽������
    endsubroutine single_channal
    
    subroutine input(xf,xg,xs,length,pd,f,Tin,pout,Tic,uic,tmax,sigma,sigmab,alpha,nf,ng,ns,ny,nt)
      real xf,xg,xs,length,pd,f,Tin,pout,Tic,uic,tmax,sigma,sigmab,alpha
      integer nf,ng,ns,ny,nt
      
      open(unit=1,file='E:\documents\doctors degree\software\tansistant\parallel_channal\input.txt')
      read(1,*) xf,xg,xs,length,pd,nf,ng,ns,ny,f,Tin,pout,Tic,uic,tmax,nt,sigma,sigmab,alpha
      close(1)
      !open(unit=2,file='E:\documents\doctors degree\software\tansistant\parallel channal\pow.txt')
      !read(2,*) q
      !open(unit=3,file='E:\documents\doctors degree\software\tansistant\parallel channal\uin.txt')
      !read(3,*) uin
    endsubroutine input 
    
   subroutine predata(xf,xg,xs,length,pd,Tin,pout,tmax,nf,ng,ns,ny,nt,df,dg,ds,dy,rc,dt,M,N,area,wet,de)!�м����Ԥ����M,N,ˮ��ѧ����,(��������ռ�)
     real xf,xg,xs,length,pd,Tin,pout,tmax,area,wet,de,df,dg,ds,dy,rc,dt
     integer nf,ng,ns,ny,nt,M,N
     !�м��������
     Df=Xf/Nf
     Dg=Xg/Ng
     Ds=Xs/Ns
     Dy=Length/Ny
     rc=xf+xg+xs
     dt=Tmax/Nt
     M=Ny+1
     N=Nf+Ng+Ns+1
     !ˮ���������� wet,de
     call get_hyconstant(rc,pd,area,wet,de)
   endsubroutine predata
    
   subroutine initial(Nf,Ng,Ns,Ny,M,N,Nt,RHO,RHOI,RHOF,RHOFI,SHC,CTC,DVS,Ti,Tic,Tin,ui,uic,uin,uint,ulast,dt,p,pin,pout,q,T)
         integer  Nf,Ng,Ns,Ny,Nt,M,N,i,j,k
         real RHO(0:M,0:N),RHOI(0:M,0:N),RHOF(0:M),RHOFI(0:M),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N),T(0:M,0:N),Ti(0:M,0:N),q(1:M-1,1:N)
         real Tin,pin,pout,uin,dt,Tic,uic
         real ui(1:Ny-1),p(1:Ny),uint(0:Nt),ulast(1:Ny-1)
         !k=0!����Tminʱ�̱߽���������
         !���ó�ʼѹ����
         do i=1,Ny,1
             if (i<Ny)then
                 p(i)=60000.0-5000.0*(i-1)
             elseif(i==Ny)then
                 p(i)=(2*pout+p(i-1))/3.0
             endif
         enddo
         !��������ٶȺ���
          open(unit=1,file='E:\documents\doctors degree\software\tansistant\parallel_channal\uin.txt')
          read(1,*) uint
          close(1)
         !���ó�ʼ����ٶ�
         uin=uint(0)
         !call get_uin(k,dt,uini,uin)
         !�������ѹ��
         pin=1.50*p(1)-0.5*p(2)
         !���ó�ʼ�¶�
         Ti(0:M,0:N)=Tic
         !��ʼ������¶�
         T(0,N)=Tin
         !������Դ
         q=0.0
         Do i=1,M-1,1
             Do j=1,N,1
                if(j<=Nf)  q(i,j)=2.827*1e7
             enddo
         enddo
         !print *,q
         !���ٶȸ�ֵ tminʱ�̵��ٶ�
         ui=uic
         ulast=ui
         !��ʼ�����Բ������¶ȸ�ֵ
         Do i=0,M,1
             Do j=0,N,1
                 if (j>=0.and.j<=Nf)then !о�������� UO2
                  RHOI(i,j)=10980
                  RHO(i,j)=10980
                  SHC(i,j)=300.0
                  CTC(i,j)=4.33
                  DVS(i,j)=0.0
                 elseif (j>Nf.and.j<=Nf+Ng) then!��϶������ He
                  RHOI(i,j)=1.785
                  RHO(i,j)=1.785
                  SHC(i,j)=1.260
                  CTC(i,j)=0.124
                  DVS(i,j)=0.0
                 elseif (j>Nf+Ng.and.j<=Nf+Ng+Ns)then!���������� Ti
                  RHOI(i,j)=7900
                  RHO(i,j)=7900
                  SHC(i,j)=502.42
                  CTC(i,j)=18.84
                  DVS(i,j)=0.0
                 else!��������
                  RHOI(i,j)=10470  
                  RHO(i,j)=10470
                  SHC(i,j)=159
                  CTC(i,j)=3.61
                  DVS(i,j)=5.0 !����ճ�ȣ������˶�ճ��
                 endif                 
             enddo
         enddo    
         do i=0,M,1
           RHOF(i)=RHO(i,N)
           RHOFI(i)=RHOI(i,N)
         enddo
    end subroutine initial
    
    subroutine grid(Xf,Xg,Xs,Length,Nf,Ng,Ns,Ny,M,N,XX,YY)
     real Xf,Xg,Xs,Length,Df,Dg,Ds,Dy 
     integer  Nf,Ng,Ns,Ny,M,N,i,j
     real XX(0:M,0:N)
     real YY(0:M,0:N)
     
     Df=Xf/Nf
     Dg=Xg/Ng
     Ds=Xs/Ns
     Dy=Length/Ny

     
     Do i=0,M,1
         do j=0,N,1
            if (j==0)then
               XX(i,j)=0.0
            elseif(j==1)then
               XX(i,j)=Df/2.0
            elseif(j>1.and.j<=Nf) then
               XX(i,j)=XX(i,j-1)+Df
            elseif(j==Nf+1)then
               XX(i,j)=XX(i,j-1)+Df/2.0+Dg/2.0
            elseif (j>Nf+1.and.j<=Nf+Ng) then
               XX(i,j)=XX(i,j-1)+Dg
            elseif (j==Nf+Ng+1)then
               XX(i,j)=XX(i,j-1)+ Dg/2.0+Ds/2.0
            elseif (j>Nf+Ng+1.and.j<=Nf+Ng+Ns)then
               XX(i,j)=XX(i,j-1)+Ds
            else!����ľ������꣬û��ʵ������
               XX(i,j)=XX(i,j-1)+Ds
            endif
            
            if(i==0)then
              YY(i,j)=0.0
            elseif (i==1)then
              YY(i,j)=Dy/2.0
            elseif (i>1.and.i<M)then
              YY(i,j)=YY(i-1,j)+Dy
            elseif (i==M)then
              YY(i,j)=YY(i-1,j)+Dy/2.0
            endif
         enddo
      enddo   
    endsubroutine grid
    
     subroutine solve_momentum(Ny,f,De,dx,dt,p,pout,uin,ui,ulast,RHO,RHOI,ap,u,uout)!appΪaP���飬��Ҫ���Ϊѹ�������������á�
      integer N,Ny,i
      real f,De,dx,dt,uin,uout,pout,pin
      !real,allocatable::p(:),ui(:),ulast(:),u(:),RHO(:),RHOI(:)
      real p(1:Ny),ui(1:Ny-1),ulast(1:Ny-1),u(1:Ny-1),RHO(0:Ny+1),RHOI(0:Ny+1),ap(1:Ny-1)
      real,allocatable::aw(:),ae(:),api(:),bs(:),A(:,:),b(:),uc(:)
     
      N=Ny-1!ʸ�����������
      !allocate(p(0:Ny+1),ui(0:Ny),ulast(0:Ny),u(0:Ny),RHO(0:Ny+1),RHOI(0:Ny+1))
      allocate(aw(1:N),ae(1:N),api(1:N),bs(1:N),A(1:N,1:N),b(1:N),uc(1:N))
      pin=(3*p(1)-p(2))/2.0
     !�������������ĳ�ϵ����Դ��
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
      !��ֱ�ӷ�����󣬵�һ����д��ϵ������
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
    
      !�ڶ���������TDMA�������
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
    
    subroutine solve_pressureCorrection(Ny,ap,RHO,RHOI,ulast,uin,pout,dx,dt,pm,b)
       integer Ny,i
       real uin,pout,dx,dt
       real ap(1:Ny-1),RHO(0:Ny+1),RHOI(0:Ny+1),ulast(1:Ny-1),pm(1:Ny),bp(1:Ny),be(1:Ny),bw(1:Ny),b(1:Ny),A(1:Ny-1,1:Ny-1),bb(1:Ny-1),pmm(1:Ny-1)
       
       do i=1,Ny,1
           if(i==1)then
             be(i)=1.5*(RHO(i)+RHO(i+1))/2.0/ap(i)
             bw(i)=0.0
             bp(i)=be(i)
             b(i)=RHO(0)*uin-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt
           elseif(i==2)then
             be(i)=(RHO(i)+RHO(i+1))/2.0/ap(i)
             bw(i)=1.5*(RHO(i)+RHO(i-1))/2.0/ap(i-1)
             bp(i)=be(i)+bw(i)
             b(i)=(RHO(i)+RHO(i-1))/2.0*ulast(i-1)-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt               
           elseif(i>2.and.i<Ny-1)then
             be(i)=(RHO(i)+RHO(i+1))/2.0/ap(i)
             bw(i)=(RHO(i)+RHO(i-1))/2.0/ap(i-1)
             bp(i)=be(i)+bw(i)
             b(i)=(RHO(i)+RHO(i-1))/2.0*ulast(i-1)-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt    
           elseif(i==Ny-1)then
             be(i)=0.0
             bw(i)=(RHO(i)+RHO(i-1))/2.0/ap(i-1)
             bp(i)=bw(i)+(RHO(i)+RHO(i+1))/2.0/ap(i)
             b(i)=(RHO(i)+RHO(i-1))/2.0*ulast(i-1)-(RHO(i)+RHO(i+1))/2.0*ulast(i)+(RHOI(i)-RHO(i))*dx/dt  
           elseif(i==Ny)then
             be(i)=0.0
             bw(i)=0.0
             bp(i)=0.0
             b(i)=0.0
           endif
       enddo
       
       A=0.0
       do i=1,Ny-1,1
           if(i==1)then
               A(i,i)=bp(i)
               A(i,i+1)=-be(i)
               bb(i)=b(i)
           elseif(i>1.and.i<Ny-1)then
               A(i,i)=bp(i)
               A(i,i+1)=-be(i)
               A(i,i-1)=-bw(i)
               bb(i)=b(i)
           elseif(i==Ny-1)then
               A(i,i)=bp(i)
               A(i,i-1)=-bw(i)
               bb(i)=b(i)
           endif
       enddo
        call tdma(Ny-1,A,bb,pmm)
        do i=1,Ny,1
            if(i<Ny)then
                pm(i)=pmm(i)
            else
                pm(i)=0.0
            endif
        enddo
       !print *,pmm
    end subroutine solve_pressureCorrection
    
  subroutine solve_temperature(De,Area,wet,Xf,Xg,Xs,Length,Nf,Ng,Ns,Ny,M,N,RHO,RHOI,SHC,CTC,DVS,htc,dt,Tin,Ti,uin,u,q,T)
    real De,Area,Xf,Xg,Xs,Xt,Length,Df,Dg,Ds,Dy,dt,wet,uin,velocity,Tin
    integer  Nf,Ng,Ns,Ny,M,N,i,j,k
   
    real T(0:M,0:N),Ti(0:M,0:N),Tj(1:M-1,1:N)
    real u(1:Ny-1),htc(0:M)
    real RHO(0:M,0:N),RHOI(0:M,0:N),SHC(0:M,0:N),CTC(0:M,0:N),DVS(0:M,0:N)
    real aw(1:M-1,1:N),ae(1:M-1,1:N),ap(1:M-1,1:N),as(1:M-1,1:N),an(1:M-1,1:N),api(1:M-1,1:N),bs(1:M-1,1:N),q(1:M-1,1:N)
    !print *,uin,u
    !�����������ϵ��
    do i=1,M-1,1
      if (i==1)then
          velocity=u(1)
      elseif(i>1.and.i<Ny)then
          velocity=(u(i-1)+u(i))/2.0
      else
          velocity=u(i-1)
      endif
      call get_convection(De,Area,wet,RHO(i,N),velocity,DVS(i,N),SHC(i,N),CTC(i,N),htc(i))!DVS(i,N)����ճ�� Pa*s
    enddo
    htc(0)=htc(1)!�߽��ϵĶ�������ϵ��
    htc(M)=htc(M-1)
    !������
     Df=Xf/Nf
     Dg=Xg/Ng
     Ds=Xs/Ns
     Dy=Length/Ny
    
    Do i=1,M-1,1
        Do j=1,N,1
         if (j==1)then!��ԳƱ߽�Ŀ�����
          aw(i,j)=0.0
          ae(i,j)=(XX(i,j)+Df/2.0)*CTC(i,j)/Df
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Df/dt
          ap(i,j)=ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Df*q(i,j)
         elseif (j>1.and.j<Nf)then!fuel�ڲ�������
          aw(i,j)=(XX(i,j)-Df/2.0)*CTC(i,j)/Df
          ae(i,j)=(XX(i,j)+Df/2.0)*CTC(i,j)/Df
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Df/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Df*q(i,j)
         elseif(j==Nf)then!f-g�߽���������
          aw(i,j)=(XX(i,j)-Df/2.0)*CTC(i,j)/Df
          ae(i,j)=2*(XX(i,j)+Df/2.0)*(Df+Dg)/(Df/CTC(i,j)+Dg/CTC(i,j+1))/(Df+Dg)
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Df/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Df*q(i,j)  
         elseif (j==Nf+1)then!f-g�߽��Ҳ������
          aw(i,j)=2*(XX(i,j)-Dg/2.0)*(df+dg)/(df/CTC(i,j-1)+dg/CTC(i,j))/(df+dg)
          ae(i,j)=(XX(i,j)+Dg/2.0)*CTC(i,j)/Dg
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Dg*q(i,j)
         elseif(j>Nf+1.and.j<Nf+Ng)then!g��϶�ڲ�������
          aw(i,j)=(XX(i,j)-Dg/2.0)*CTC(i,j)/Dg
          ae(i,j)=(XX(i,j)+Dg/2.0)*CTC(i,j)/Dg
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Dg*q(i,j)
         elseif(j==Nf+Ng)then!g-c�߽���������
          aw(i,j)=(XX(i,j)-Dg/2.0)*CTC(i,j)/Dg
          ae(i,j)=2*(XX(i,j)+Dg/2.0)*(Dg+Ds)/(Dg/CTC(i,j)+Ds/CTC(i,j+1))/(Dg+Ds)
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Dg/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Dg*q(i,j)
         elseif(j==Nf+Ng+1)then!g-c�߽��Ҳ������
          aw(i,j)=2*(XX(i,j)-Ds/2.0)*(dg+ds)/(dg/CTC(i,j-1)+ds/CTC(i,j))/(dg+ds)
          ae(i,j)=(XX(i,j)+Ds/2.0)*CTC(i,j)/Ds
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Ds*q(i,j)
         elseif(j>Nf+Ng+1.and.j<Nf+Ng+Ns)then!c�����ڲ�������
          aw(i,j)=(XX(i,j)-Ds/2.0)*CTC(i,j)/Ds
          ae(i,j)=(XX(i,j)+Ds/2.0)*CTC(i,j)/Ds
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Ds*q(i,j)
         elseif(j==Nf+Ng+Ns)then!s-fluid�߽���������
          aw(i,j)=(XX(i,j)-Ds/2.0)*CTC(i,j)/Ds
          ae(i,j)=(XX(i,j)+Ds/2.0)/(1.0/htc(i)+ds/(2.0*CTC(i,j)))
          as(i,j)=0.0
          an(i,j)=0.0
          api(i,j)=RHO(i,j)*SHC(i,j)*XX(i,j)*Ds/dt
          ap(i,j)=aw(i,j)+ae(i,j)+api(i,j)
          bs(i,j)=XX(i,j)*Ds*q(i,j)
         elseif(j==Nf+Ng+Ns+1)then!fluid������
          Xt=Xf+Xg+Xs !�����⾶
          if(i==1)then!������ڵĿ�����
           aw(i,j)=Dy/SHC(i,j)*2.0*3.14*Xt/Area*1.0/(1.0/HTC(i)+Ds/(2*CTC(i,j-1)))
           ae(i,j)=0.0
           as(i,j)=0.0
           an(i,j)=0.0
           api(i,j)=RHOI(i,j)*Dy/dt
           ap(i,j)=aw(i,j)+api(i,j)+RHO(i-1,j)*uin
           bs(i,j)=Dy/SHC(i,j)*q(i,j)+RHO(i-1,j)*uin*T(i-1,j)          
          else!�����ڲ��Լ����ڿ�����
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

!����ע�Ͳ���Ϊ�������Դ��ľ�����⣬������ֱ�ӷ���⣬�ڵ������в���Ҫ��    
!    A=0.0
!    b=0.0
!    Do i=1,M-1,1
!     Do j=1,N,1
!      k=(i-1)*N+j
!      if(j==1)then
!        A(k,k)=ap(i,j)
!        A(k,k+1)=-ae(i,j)
!     elseif(j>1.and.j<N)then
!        A(k,k)=ap(i,j)
!        A(k,k-1)=-aw(i,j)
!        A(k,k+1)=-ae(i,j)
!     elseif(j==N)then
!       if(i==1)then
!        A(k,k)=ap(i,j)
!        A(k,k-1)=-aw(i,j)
!        A(k,k+N)=-as(i,j)
!       elseif(i>1.and.i<M-1)then
!        A(k,k)=ap(i,j)
!        A(k,k-1)=-aw(i,j)
!        A(k,k+N)=-as(i,j)
!        A(k,k-N)=-an(i,j)
!       elseif(i==M-1)then
!        A(k,k)=ap(i,j)
!        A(k,k-1)=-aw(i,j)
!        A(k,k-N)=-an(i,j)
!       endif
!      endif 
!        b(k)=bs(i,j)+api(i,j)*Ti(i,j)
!     enddo
!    enddo

!�ſɱȵ���������
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
                  T(i,j)=Tj(i,j+1)
                else
                  T(i,j)=Tj(i,j)
                endif
            elseif(i==0)then
                if(j==0)then
                  T(i,j)=Tj(i+1,j+1)
                elseif(j>0.and.j<N)then
                  T(i,j)=Tj(i+1,j)
                elseif(j==N)then!���
                  T(i,j)=Tin
                endif
           elseif(i==M)then
                if(j==0)then
                  T(i,j)=Tj(i-1,j+1)
                else
                  T(i,j)=Tj(i-1,j)
                endif
            endif           
         enddo
     enddo
        !print *,Tj!����������
        !print *,dvs(:,N)

    end subroutine 
    
    subroutine get_hyconstant(rc,pd,Aflow,wet,de)
       real rc,p,pd !r�ǰ�����뾶 p�ǶԱ߾�
       real Aflow,Ashell,Atotal,wet,de
       
           p=pd*2*rc
           Ashell=3.14*rc*rc
           Atotal=0.5*sqrt(3.0)*p*p
           Aflow=Atotal-Ashell
           wet=2*3.14*rc+p/sqrt(3.0)*6
           de=4*Aflow/wet
    endsubroutine get_hyconstant
    
    subroutine get_Nusselt(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,Nu)!Һ̬�ؽ�����Ŭ������
     real flow_area
     real wetted_perimeter    
     real density
     real velocity
     real De
     real Re 
     real viscosity
     real capacity
     real conductivity
     real Pr
     real Pe
     real Nu
     
     De=4*flow_area/wetted_perimeter
     Re=4*density*velocity*De/viscosity
     Pr=viscosity*capacity/conductivity
     Pe=Re*Pr
     Nu=5.0+2.5D-2*Pe**0.8
    end subroutine get_Nusselt
    
    subroutine get_convection(lenth,flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,convection) !lenth����������,ȡˮ��ѧֱ��  Nu=(h*x)/k
     real lenth
     real flow_area
     real wetted_perimeter
     real density
     real velocity
     real De
     real viscosity,capacity,convection,conductivity
     real Pr,Re,Pe,Nu
     
     call get_Nusselt(flow_area,wetted_perimeter,density,velocity,viscosity,capacity,conductivity,Nu)  
     convection=Nu*conductivity/lenth
    end subroutine get_convection
    
    subroutine CFD_post(Nf,Ng,Ns,Ny,M,N,Nt,T,uin,u,uout,pin,p,pout,RHO,Ttout,RHOtout,Ttout_fuel,Tfg,Tgs,Tsf)!Post_rho_fluid_outlet,Post_T_fuel_outlet
      integer Nf,Ng,Ns,Ny,M,N,Nt
      real uin,uout,pin,pout
      real T(0:M,0:N),Ttout(1:Nt),RHOtout(1:Nt),Ttout_fuel(1:Nt),u(1:Ny-1),p(1:Ny),RHO(0:M,0:N),Tfg(0:M),Tgs(0:M),Tsf(0:M)
     !print *,
     !print *,T(:,0)
      open(1,file='E:\documents\doctors degree\software\tansistant\single\output\temperature.txt')
      write(1,*)T
      close(1)    
      open(2,file='E:\documents\doctors degree\software\tansistant\single\output\N.txt')
      write(2,*)Nf,Ng,Ns,Ny,M,N,Nt
      close(2)      
      !print *,T(:,0)
      open(3,file='E:\documents\doctors degree\software\tansistant\single\output\tt.txt')
      write(3,*)Ttout,RHOtout,Ttout_fuel
      close(3)
      
      open(4,file='E:\documents\doctors degree\software\tansistant\single\output\pressure.txt')
      write(4,*) pin,p,pout
      close(4)
      
      open(5,file='E:\documents\doctors degree\software\tansistant\single\output\density.txt')
      write(5,*) RHO
      close(5)
      
      open(6,file='E:\documents\doctors degree\software\tansistant\single\output\velocity.txt')
      write(6,*) uin,u,uout
      close(6)
      
      open(7,file='E:\documents\doctors degree\software\tansistant\single\output\Tboundary.txt')
      write(7,*) Tfg,Tgs,Tsf
      close(7)
    end subroutine CFD_post
     
     subroutine tdma(N,A,B,u)
      real A(N,N)
      real B(N)
      real u(N)
      integer N
      real,dimension(N)::aa,bb,cc,dd,x,y

      do i=1,N,1
          if(i==1)then
          aa(1)=0.0
          bb(1)=A(1,1)
          cc(1)=A(1,2)
          dd(1)=B(1)
          elseif (i==N) then
           aa(N)=A(N,N-1)
           bb(N)=A(N,N)
           cc(N)=0.0
           dd(N)=B(N)
          else
           aa(i)=A(i,i-1)
           bb(i)=A(i,i)
           cc(i)=A(i,i+1)
           dd(i)=B(i)
          endif
      enddo
      
      x(1)=cc(1)/bb(1)
      y(1)=dd(1)/bb(1)
      
      do i=2,N,1
          x(i)=cc(i)/(bb(i)-x(i-1)*aa(i))
          y(i)=(dd(i)-y(i-1)*aa(i))/(bb(i)-x(i-1)*aa(i))
      enddo
      
      u(N)=y(N)
      
      do i=N-1,1,-1
          u(i)=y(i)-x(i)*u(i+1)
      enddo
     endsubroutine tdma    
     
!     subroutine get_uin(k,dt,uini,uin)
!        integer k
!        real dt,uini,uin
!        if (k<=80)then
!         uin=uini
!         else
!          if (uin>0.3*uini)then 
!            uin=uini/(1+(k-80)*dt/5.5)
!          else
!           uin=0.3*uini
!          endif
!         endif
!     end subroutine get_uin
    end program Console84

