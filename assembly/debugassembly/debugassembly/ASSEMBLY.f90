
!****************************************************************************
!
!  PROGRAM: debugassembly
!
!  PURPOSE:  Entry point for the console application.
!
!  pow(na,nr),fq_core(na,nr)     ƽ�������ܶȣ����ʷ�����
!  nr = SIZE(assembly, dim=1)    ����������Ŀ
!  na = SIZE(assembly, dim=2)    ����Ľڿ���Ŀ��ԭ�����������Ҫ������ֵ�������þֲ��������
!
!****************************************************************************

    program debugassembly
     !����ʵ����reInputdata/assm1
     use sys_re_input_global
     use assm_global
     use driving_pre_process
     use driving_output
    implicit none
     !Ԥ����
     call sys_pre_process()
     !����
     !call assm1%steady(power,fq_core)
     !call assm1%transient(power, fq_core, tidx, ltime, ctime)
     
     !���
     call Run_output()     

    end program debugassembly

