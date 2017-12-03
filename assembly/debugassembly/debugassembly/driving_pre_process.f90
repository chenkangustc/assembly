module driving_pre_process
    use assm_global
    implicit none
    private
    public::Sys_pre_process
contains
    subroutine Sys_pre_process()
     implicit none
     write(*,*)'start the sys pre process:'
     !读取参数
     call reInputdata%set()
     !call reInputdata%publish()
     !参数赋值
     call assm1%set(reInputdata)
     !call assm1%geom%print()
     !分配空间
     call assm1%alloc()
     !初始化
     call assm1%init()
    end subroutine Sys_pre_process
end module driving_pre_process