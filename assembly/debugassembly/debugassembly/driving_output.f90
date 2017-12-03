module driving_output
    use assm_global
    implicit none
    private
    public::Run_output
contains
    subroutine Run_output()
      open(1,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\mesh.txt')
      write(1,*) assm1%mesh%Nf,assm1%mesh%Ng,assm1%mesh%Ns,assm1%mesh%Ny
      close(1) 
      open(2,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\output\temperature.txt')
      write(2,*) assm1%thermal%temperature
      !write(2,*) assm1%pow%power
      close(2) 

      !print*,assm1%th_boundary%p%outlet
      !liquid PVT¡¢rho distribution
      !solid temperature distribution
      print*,assm1%th_boundary%u%inlet,assm1%Thermal%velocity,assm1%th_boundary%u%outlet
      print*,assm1%th_boundary%p%inlet,assm1%Thermal%pressure,assm1%th_boundary%p%outlet
      print*,assm1%th_boundary%T%inlet,assm1%Thermal%temperature(:,assm1%mesh%Nf+assm1%mesh%Ng+assm1%mesh%Ns+1),assm1%th_boundary%T%outlet
    end subroutine Run_output
end module driving_output
