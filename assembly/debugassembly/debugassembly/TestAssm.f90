
!****************************************************************************
!
!  PROGRAM: debugassembly
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program debugassembly
     use sys_re_input_global
    implicit none
     type(sys_re_input)::input1
     call reInputdata=>set()
     ! Body of debugassembly
     call cal_Assembly_steady(this,power,fq_core)
     call cal_Assembly_Transient(this,power, fq_core, tidx, ltime, ctime)!pow(ia,ir),fq_core(ia,ir)
    
     print *, 'Hello World'

    end program debugassembly

