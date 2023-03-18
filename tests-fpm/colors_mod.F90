!*****************************************************************************************
module mod_colors
    implicit none
    character(len=*), parameter:: red_failed=achar(27)//'[31mFAILED:'//achar(27)//'[0m'
    character(len=*), parameter:: green_passed=achar(27)//'[32mPASSED:'//achar(27)//'[0m'
end module mod_colors
!*****************************************************************************************
