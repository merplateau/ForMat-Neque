program main_with_matlab
    use, intrinsic :: iso_c_binding
    implicit none
    
    ! MATLAB Engine接口
    interface
        function engOpen(startcmd) bind(c, name='engOpen')
            use, intrinsic :: iso_c_binding
            character(kind=c_char), dimension(*) :: startcmd
            type(c_ptr) :: engOpen
        end function engOpen
        
        function engClose(ep) bind(c, name='engClose')
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: ep
            integer(c_int) :: engClose
        end function engClose
        
        function engEvalString(ep, string) bind(c, name='engEvalString')
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: ep
            character(kind=c_char), dimension(*) :: string
            integer(c_int) :: engEvalString
        end function engEvalString
        
        function engGetVariable(ep, name) bind(c, name='engGetVariable')
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: ep
            character(kind=c_char), dimension(*) :: name
            type(c_ptr) :: engGetVariable
        end function engGetVariable
    end interface
    
    ! 变量声明
    type(c_ptr) :: ep
    integer(c_int) :: status
    character(len=256) :: cmd
    real(c_double) :: result_real, result_imag
    
    print *, "启动MATLAB引擎..."
    
    ! 启动MATLAB引擎
    ep = engOpen("" // c_null_char)
    if (.not. c_associated(ep)) then
        print *, "错误：无法启动MATLAB引擎"
        stop
    end if
    
    print *, "MATLAB引擎启动成功"
    
    ! 设置工作目录到当前目录
    cmd = "cd('" // trim(getcwd()) // "')" // c_null_char
    status = engEvalString(ep, cmd)
    
    ! 运行k_solver.m
    print *, "运行k_solver.m..."
    cmd = "k_solver" // c_null_char
    status = engEvalString(ep, cmd)
    
    if (status /= 0) then
        print *, "错误：无法执行k_solver.m"
    else
        print *, "k_solver.m执行成功"
    end if
    
    ! 获取结果（如果需要的话）
    ! 这里可以添加代码来获取MATLAB变量
    
    ! 关闭MATLAB引擎
    status = engClose(ep)
    if (status == 0) then
        print *, "MATLAB引擎已关闭"
    else
        print *, "警告：关闭MATLAB引擎时出现问题"
    end if
    
    print *, "程序执行完成"
    
contains
    function getcwd() result(path)
        character(len=256) :: path
        call getcwd(path)
    end function getcwd
end program main_with_matlab 