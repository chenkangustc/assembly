module assm_global
    use sys_re_input_header
    implicit none
    type(sys_re_input)::indata
    call indata%set()
end module assm_global