module driving_pre_process
    use assm_global
    implicit none
    private
    public::Sys_pre_process
contains
    subroutine Sys_pre_process()
     !��ȡ����
     call reInputdata%set()
     !������ֵ
     call assm1%set(reInputdata)
     !����ռ�
     call assm1%alloc()
     !��ʼ��
     call assm1%init()
    end subroutine Sys_pre_process
end module driving_pre_process