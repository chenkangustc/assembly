module driving_pre_process
    use assm_global
    implicit none
    private
    public::Sys_pre_process
contains
    subroutine Sys_pre_process()
     implicit none
     write(*,*)'start the sys pre process:'
     !��ȡ����
     call reInputdata%set()
     !call reInputdata%publish()
     !������ֵ
     call assm1%set(reInputdata)
     !call assm1%geom%print()
     !����ռ�
     call assm1%alloc()
     !��ʼ��
     call assm1%init()
    end subroutine Sys_pre_process
end module driving_pre_process