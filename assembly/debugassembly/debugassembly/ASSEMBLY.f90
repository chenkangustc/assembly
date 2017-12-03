!***************************************************************************************
!
!  PROGRAM: debugassembly
!
!  PURPOSE:  Entry point for the console application.
!
!  pow(na,nr),fq_core(na,nr)     平均功率密度，功率峰因子
!  nr = SIZE(assembly, dim=1)    径向的组件数目
!  na = SIZE(assembly, dim=2)    轴向的节块数目，原输入变量不需要操作赋值的另外用局部变量表达
!
!***************************************************************************************
    program debugassembly
     !两个实例：reInputdata/assm1
     use sys_re_input_global
     use assm_global
     use driving_pre_process
     use driving_output
     use sys_power_header
    implicit none
     real ltime
     real ctime
     real ttotal
     real dt
     real,allocatable::power(:,:),fq_core(:,:)
     real,allocatable::time_array(:),tout(:),tpow(:),tuin(:)
     integer M,N,Nt,i,j
     !pre_process
     call sys_pre_process()
     !*******************************************
     ttotal=150.0
     dt=1.0
     M=size(assm1%thermal%temperature,dim=1)
     N=size(assm1%thermal%temperature,dim=2)
     Nt=ttotal/dt
     allocate(power(M,N),fq_core(M,N))
     allocate(time_array(Nt),tout(Nt),tpow(Nt),tuin(Nt)) 
     fq_core=0.0
     power=0.0
     do i=1,M,1
         do j=1,N,1
          if(j<=assm1%mesh%Nf) power(i,j)=2.827*1e7
         enddo
     enddo
     !*********************************************
     call assm1%steady(power,fq_core)!power come from other data,so it should be an interface in place with the data
     ltime=0.0
     ctime=0.0
     i=1
     do while(ctime<ttotal)
     !call assm1%transient(power, fq_core, tidx, ltime, ctime)
      time_array(i)=ctime
      tout(i)=assm1%th_boundary%T%outlet
      tpow(i)=power(1,1)
      tuin(i)=assm1%th_boundary%u%inlet
      i=i+1
      ctime=ctime+dt
      call update_power(power,fq_core,ltime,ctime)
      call assm1%transient(power, fq_core,ltime,ctime)
      print*,'ctime=',ctime
      !print*,assm1%th_boundary%T%inlet,assm1%Thermal%temperature(:,1),assm1%th_boundary%T%outlet
      ltime=ctime
     !output
     enddo
     call Run_output() 
      
      open(1,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\ctime.txt')
      write(1,*) time_array
      close(1)
      open(2,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tout.txt')
      write(2,*) tout
      close(2) 
      open(3,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tpow.txt')
      write(3,*) tpow
      close(3) 
      open(4,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\tuin.txt')
      write(4,*) tuin
      close(4)

    end program debugassembly

