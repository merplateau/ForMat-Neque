program main_system_call
    implicit none
    
    integer :: status
    character(len=256) :: cmd
    
    print *, "使用系统调用方式调用MATLAB..."
    
    ! 方法1：直接调用MATLAB运行脚本
    cmd = "matlab -batch 'k_solver'"
    print *, "执行命令: ", trim(cmd)
    
    call system(cmd, status)
    
    if (status == 0) then
        print *, "MATLAB脚本执行成功"
    else
        print *, "MATLAB脚本执行失败，状态码: ", status
    end if
    
    ! 方法2：创建临时MATLAB脚本并执行
    print *, "创建临时MATLAB脚本..."
    
    ! 创建临时脚本文件
    open(unit=10, file='temp_script.m', status='replace')
    write(10, *) "cd('" // trim(getcwd()) // "')"
    write(10, *) "k_solver"
    write(10, *) "exit"
    close(10)
    
    ! 执行临时脚本
    cmd = "matlab -batch 'temp_script'"
    print *, "执行临时脚本: ", trim(cmd)
    
    call system(cmd, status)
    
    if (status == 0) then
        print *, "临时脚本执行成功"
    else
        print *, "临时脚本执行失败，状态码: ", status
    end if
    
    ! 清理临时文件
    call system("rm -f temp_script.m")
    
    print *, "程序执行完成"
    
contains
    function getcwd() result(path)
        character(len=256) :: path
        call getcwd(path)
    end function getcwd
end program main_system_call 