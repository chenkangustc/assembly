module driving_cal_transient
    use assm_global
    use mathkerel
    implicit none
    private
    public::cal_Assembly_transient
    !public::solve_momentum
    public::cal_momentumA
    !public::solve_momentumA
contains
    subroutine cal_Assembly_transient(assm)
      implicit none
      !class(sys_assembly),intent(in out)::this 
      type(sys_assembly),intent(in out)::assm
      !real,intent(in)::power(:)
      !real,intent(in)::fq_core(:)
      integer Ny
      !Ny=SIZE(power)

    end subroutine cal_Assembly_transient
    
    subroutine cal_momentumA(assm,pguess,uinit,uboundary,dt,A,b)
     implicit none
     class(sys_assembly),intent(in out)::assm
     real,intent(in)::pguess(:)
     real,intent(in)::uinit(:)
     type(boundary),intent(in)::uboundary
     real,intent(in)::dt
     real,intent(in out)::A(:,:)
     real,intent(in out)::b(:)
     
     
    end subroutine cal_momentumA

    
    
end module driving_cal_transient    