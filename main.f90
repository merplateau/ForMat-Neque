program main
    use, intrinsic :: iso_c_binding
    implicit none
    
    ! 声明外部函数接口
    interface
        function solvek(M, ni, ne, B0, f, Te, Ti, Tperp_e, Tpara_e, Tperp_i, Tpara_i, V_e, V_i) &
                     bind(c, name='solvek')
            use, intrinsic :: iso_c_binding
            real(c_double), value :: M, ni, ne, B0, f, Te, Ti, Tperp_e, Tpara_e, Tperp_i, Tpara_i, V_e, V_i
            complex(c_double) :: solvek
        end function solvek
    end interface
    
    ! 变量声明
    real(c_double) :: M, ni, ne, B0, f, Te, Ti, Tperp_e, Tpara_e, Tperp_i, Tpara_i, V_e, V_i
    complex(c_double) :: kpar_sol
    real(c_double) :: kpar_real, kpar_imag
    
    ! 设置参数值
    M = 40.0_c_double        ! 离子质量数
    ni = 1.0e18_c_double     ! 离子密度 (m^-3)
    ne = 1.0e18_c_double     ! 电子密度 (m^-3)
    B0 = 10.0_c_double       ! 磁场 (T)
    f = 13.56e9_c_double     ! 驱动频率 (Hz)
    Te = 30.0_c_double       ! 电子温度 (K)
    Ti = 3.0_c_double        ! 离子温度 (K)
    Tperp_e = 30.0_c_double  ! 电子垂直温度 (K)
    Tpara_e = 30.0_c_double  ! 电子平行温度 (K)
    Tperp_i = 3.0_c_double   ! 离子垂直温度 (K)
    Tpara_i = 3.0_c_double   ! 离子平行温度 (K)
    V_e = 0.0_c_double       ! 电子漂移速度 (m/s)
    V_i = 0.0_c_double       ! 离子漂移速度 (m/s)
    
    print *, "=== Fortran调用MATLAB函数测试 ==="
    print *, "参数设置："
    print *, "M = ", M
    print *, "ni = ne = ", ni, " m^-3"
    print *, "B0 = ", B0, " T"
    print *, "f = ", f, " Hz"
    print *, "Te = ", Te, " K"
    print *, "Ti = ", Ti, " K"
    print *, "Tperp_e = Tpara_e = ", Tperp_e, " K"
    print *, "Tperp_i = Tpara_i = ", Tperp_i, " K"
    print *, "V_e = V_i = ", V_e, " m/s"
    print *
    
    ! 调用MATLAB函数
    print *, "调用MATLAB函数 solvek..."
    kpar_sol = solvek(M, ni, ne, B0, f, Te, Ti, Tperp_e, Tpara_e, Tperp_i, Tpara_i, V_e, V_i)
    
    ! 提取实部和虚部
    kpar_real = real(kpar_sol)
    kpar_imag = aimag(kpar_sol)
    
    ! 输出结果
    print *, "计算结果："
    print *, "k_parallel = ", kpar_real, " + ", kpar_imag, "i"
    print *, "|k_parallel| = ", abs(kpar_sol)
    print *, "相位 = ", atan2(kpar_imag, kpar_real) * 180.0 / 3.14159265359, " 度"
    
    print *, "测试完成！"
    
end program main