
!****************************************************************************
!
!  PROGRAM: debugassembly
!
!  PURPOSE:  Entry point for the console application.
!
!  pow(na,nr),fq_core(na,nr)     平均功率密度，功率峰因子
!  nr = SIZE(assembly, dim=1)    径向的组件数目
!  na = SIZE(assembly, dim=2)    轴向的节块数目，原输入变量不需要操作赋值的另外用局部变量表达
!
!****************************************************************************

    program debugassembly
     !两个实例：reInputdata/assm1
     use sys_re_input_global
     use assm_global
     use driving_pre_process
     use driving_output
    implicit none
     !预处理
     call sys_pre_process()
     !计算
     !call assm1%steady(power,fq_core)
     !call assm1%transient(power, fq_core, tidx, ltime, ctime)
     
     !输出
     call Run_output()     

    end program debugassembly

