subroutine cal_Assembly_Transient(this,dt)
 class(sys_assembly),intent(in out)::this
 real dt
 real pmodify
 !local
 type(sys_time)::last  !��һʱ�̵��ȹ�
 type(thermal):: iteration  !��ǰ���������ȹ�
 type(boundary)::biter 
 real btotal,drho
 
 iteration%p=last%p !����ѹ���²�ֵ
 drho=1.0!
 do while(drho>sigma)
     btotal=1.0
  do while(btotal>sigmab)   
  call solve_momentum(this,dt,last,iteration)
  call solve_pressureCorrection(iteration,pmodify)
  call modify_PV(this,pmodify,iteration)
  !����btotal
  end do
  call solve_temperature(this,last,iteration)
  call ����drho
 end do
  call thermal set
    enddo

!********************
 solve_momentum   P��֪������£����һά˲̬��������
!********************
    subroutine solve_momentum(this,pguess,init,boundary,dt,itertion)!(��������ˮ��ѧ)��p������ʼ���������߽���������˲̬����ѡ�
        class(sys_assembly),intent(in out)::this
        real pguess,init,boundary,dt,itertion
        !local
        real A(1:N,1:N),b(1:N)
        call cal_momentumA��this,pguess,init,boundary,dt,A��b��!this���Խ�һ������ɼ����ࡣ
        call solve_momentumA��N,A,b,u_iteration��
    end solve_momentum


