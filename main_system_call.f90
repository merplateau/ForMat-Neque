program main_system_call
    implicit none
    
    integer :: status
    character(len=256) :: cmd
    
    ! print *, "使用系统调用方式调用MATLAB..."
    
    ! 方法1：直接调用MATLAB运行脚本
    cmd = 'cmd /c matlab -batch k_solver'
    ! print *, "执行命令: ", trim(cmd)
    
    call system(cmd)
    
   
end program main_system_call 