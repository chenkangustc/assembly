module driving_pre_process
    use assm_global
    implicit none
    private
    public::Sys_pre_process
contains
    subroutine Sys_pre_process()
     !读取参数
     call reInputdata%set()
     !参数赋值
     call assm1%set(reInputdata)
     !分配空间
     call assm1%alloc()
     !初始化
     call assm1%init()
    end subroutine Sys_pre_process
end module driving_pre_process