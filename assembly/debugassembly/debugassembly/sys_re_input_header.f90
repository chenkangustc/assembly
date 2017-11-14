module sys_re_input_header
    implicit none
    type,public::sys_re_input
       !public
       integer nf,ng,ns,ny,npin
       real xf,xg,xs,xos,acf,height,f,pout,flowin,sigma,sigmab,alpha
       real Tin,uin,pin
       real Ti,ui,pi
    contains
     procedure,public::set=>set_inputdata
    end type sys_re_input
     private::set_inputdata
    contains
     subroutine set_inputdata(this)
      implicit none
      class(sys_re_input)::this
       open(unit=1,file='E:\documents\doctors degree\software\tansistant\parallel\moduledebug\assembly\debugassembly\re_input.txt')
       !read(1,*) xf,xg,xs,height,nf,ng,ns,ny,f,Tin,pout,Tic,uic,tmax,nt,sigma,sigmab,alpha
       read(1,*) this%xf,this%xg,this%xs,this%xos,this%acf,this%height,this%npin,this%nf,this%ng,this%ns,this%ny,this%f,this%Tin,this%pout,this%flowin,this%sigma,this%sigmab,this%alpha,this%Ti,this%ui,this%pi,this%uin,this%pin
       close(1)
     end subroutine set_inputdata
end module sys_re_input_header