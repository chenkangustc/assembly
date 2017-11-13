module assm_global
    use sys_assembly_header
    implicit none
    type(sys_assmebly)::assm1
    call assm1%set()
end module assm_global