#---------------------
# ModFEM operating system dependencies
# for all of targets, regardless from.

add_custom_target( ModFEM_OS_target ALL )

SET(MF_OS_INFO)

if (WIN32)
    message("OS: Windows")
    #do something

elseif(UNIX)

        if(EXISTS "/etc/debian_version")
            message("OS: debian-based")
            execute_process(COMMAND lsb_release -a
                            OUTPUT_VARIABLE MF_OS_INFO)
#            message("${MF_OS_INFO}")
            if("${MF_OS_INFO}" MATCHES ".*Ubuntu 14.*")
                message("OS: linux/ubutnu 14")
                target_link_libraries(ModFEM_OS_target "-Wl,--no-as-needed -lm -ldl")
            endif()

    elseif(EXIST "/etc/redhat-release")

        message("OS: redhat-based")
    endif()

elseif(UNIX)

    message(WARNING "OS: Operating system not detected, not applying OS-depandent stuff.")

endif()
