module driving_output
    use assm_global
    implicit none
    private
    public::Run_output
contains
    subroutine Run_output()
      open(1,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\temperature.txt')
      write(1,*) assm1%Thermal%Temperature
      close(1) 
    end subroutine Run_output
end module driving_output
