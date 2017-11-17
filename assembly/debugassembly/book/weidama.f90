subroutine cal_Assembly_Transient(this,dt)
 class(sys_assembly),intent(in out)::this
 real dt
 real pmodify
 !local
 type(sys_time)::last  !上一时刻的热工
 type(thermal):: iteration  !当前迭代步的热工
 type(boundary)::biter 
 real btotal,drho
 
 iteration%p=last%p !给出压力猜测值
 drho=1.0!
 do while(drho>sigma)
     btotal=1.0
  do while(btotal>sigmab)   
  call solve_momentum(this,dt,last,iteration)
  call solve_pressureCorrection(iteration,pmodify)
  call modify_PV(this,pmodify,iteration)
  !计算btotal
  end do
  call solve_temperature(this,last,iteration)
  call 计算drho
 end do
  call thermal set
    enddo

!********************
 solve_momentum   P已知的情况下，求解一维瞬态动量方程
!********************
    subroutine solve_momentum(this,pguess,init,boundary,dt,itertion)!(几何网格水力学)（p）（初始条件）（边界条件）（瞬态计算选项）
        class(sys_assembly),intent(in out)::this
        real pguess,init,boundary,dt,itertion
        !local
        real A(1:N,1:N),b(1:N)
        call cal_momentumA（this,pguess,init,boundary,dt,A，b）!this可以进一步抽象成计算类。
        call solve_momentumA（N,A,b,u_iteration）
    end solve_momentum


