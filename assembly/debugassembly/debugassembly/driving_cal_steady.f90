module driving_cal_steady
    use assm_global
    use mathkerel
    implicit none
    private
    public::cal_Assembly_Steady
contains
    subroutine cal_Assembly_Steady(assm,power,fq_core)
      implicit none
      !class(sys_assembly),intent(in out)::this 
      type(sys_assembly),intent(in out)::assm
      real,intent(in)::power(:)
      real,intent(in)::fq_core(:)
      integer Ny
      Ny=SIZE(power,dim=1)
      
      
      
    end subroutine cal_Assembly_Steady
end odule driving_cal_steady    