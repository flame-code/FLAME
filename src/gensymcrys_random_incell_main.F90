subroutine random_incell_main(LATSGP,CRYSSYS,NGUESS,RED_POS)
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS
real(8):: RED_POS(3,NGUESS)
if(NGUESS.lt.30) write(*,*) "WARNING: to have a good sampling of the unit cell increase NGUESS"
if(LATSGP<1 .or. LATSGP>230) stop 'ERROR: LATSGP not within 1 and 230 in random_incell_main'
if(LATSGP<50) then
    call random_incell_p1(LATSGP,CRYSSYS,NGUESS,RED_POS)
elseif(LATSGP<85) then
    call random_incell_p2(LATSGP,CRYSSYS,NGUESS,RED_POS)
elseif(LATSGP<124) then
    call random_incell_p3(LATSGP,CRYSSYS,NGUESS,RED_POS)
elseif(LATSGP<141) then
    call random_incell_p4(LATSGP,CRYSSYS,NGUESS,RED_POS)
elseif(LATSGP<190) then
    call random_incell_p5(LATSGP,CRYSSYS,NGUESS,RED_POS)
else
    call random_incell_p6(LATSGP,CRYSSYS,NGUESS,RED_POS)
endif
end subroutine random_incell_main
