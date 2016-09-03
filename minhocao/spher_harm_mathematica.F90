subroutine ylm_mathematica(l,m,theta,phi,ylm_r,ylm_i)
implicit none
integer:: l,m
real(8):: theta,phi,ylm_r,ylm_i
complex(8):: Ylm


select case(l)
case(0)
   Ylm=2.820947917738781d-1
case(1)
 select case(m)
   case(-1)
   Ylm=(3.454941494713355d-1*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=4.886025119029199d-1*Cos(theta)
   case(1)
   Ylm=-3.454941494713355d-1*2.7182818284590446d0**((0.d0,1.d0)*phi)*Sin(theta)
end select
case(2)
 select case(m)
   case(-2)
   Ylm=(3.8627420202318956d-1*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(7.725484040463791d-1*Cos(theta)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=3.1539156525252006d-1*(-1.d0 + 3.d0*Cos(theta)**2)
   case(1)
   Ylm=-7.725484040463791d-1*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*Sin(theta)
   case(2)
   Ylm=3.8627420202318956d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*Sin(theta)**2
end select
case(3)
 select case(m)
   case(-3)
   Ylm=(4.172238236327841d-1*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.0219854764332823d0*Cos(theta)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(3.2318018411415066d-1*(-1.d0 + 5.d0*Cos(theta)**2)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=3.731763325901154d-1*(-3.d0*Cos(theta) + 5.d0*Cos(theta)**3)
   case(1)
   Ylm=-3.2318018411415066d-1*2.7182818284590446d0**((0.d0,1.d0)*phi)*(-1.d0 + 5.d0*Cos(theta)**2)*Sin(theta)
   case(2)
   Ylm=1.0219854764332823d0*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*Sin(theta)**2
   case(3)
   Ylm=-4.172238236327841d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*Sin(theta)**3
end select
case(4)
 select case(m)
   case(-4)
   Ylm=(4.425326924449826d-1*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(1.2516714708983523d0*Cos(theta)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(3.345232717786446d-1*(-1.d0 + 7.d0*Cos(theta)**2)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(4.7308734787878d-1*Cos(theta)*(-3.d0 + 7.d0*Cos(theta)**2)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=1.057855469152043d-1*(3.d0 - 3.d1*Cos(theta)**2 + 3.5d1*Cos(theta)**4)
   case(1)
   Ylm=-4.7308734787878d-1*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(-3.d0 + 7.d0*Cos(theta)**2)*Sin(theta)
   case(2)
   Ylm=3.345232717786446d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*(-1.d0 + 7.d0*Cos(theta)**2)*Sin(theta)**2
   case(3)
   Ylm=-1.2516714708983523d0*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*Sin(theta)**3
   case(4)
   Ylm=4.425326924449826d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*Sin(theta)**4
end select
case(5)
 select case(m)
   case(-5)
   Ylm=(4.6413220344085815d-1*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(1.467714898305751d0*Cos(theta)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(3.4594371914684023d-1*(-1.d0 + 9.d0*Cos(theta)**2)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.6947711832608994d0*Cos(theta)*(-1.d0 + 3.d0*Cos(theta)**2)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(3.2028164857621517d-1*(1.d0 - 1.4000000000000001d1*Cos(theta)**2 + 2.1d1*Cos(theta)**4)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=1.169503224534236d-1*(1.5d1*Cos(theta) - 7.d1*Cos(theta)**3 + 6.3d1*Cos(theta)**5)
   case(1)
   Ylm=-3.2028164857621517d-1*2.7182818284590446d0**((0.d0,1.d0)*phi)*(1.d0 - 1.4000000000000001d1*Cos(theta)**2 + 2.1d1*Cos(theta)**4)*Sin(theta)
   case(2)
   Ylm=1.6947711832608994d0*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(-1.d0 + 3.d0*Cos(theta)**2)*Sin(theta)**2
   case(3)
   Ylm=-3.4594371914684023d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*(-1.d0 + 9.d0*Cos(theta)**2)*Sin(theta)**3
   case(4)
   Ylm=1.467714898305751d0*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*Sin(theta)**4
   case(5)
   Ylm=-4.6413220344085815d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*Sin(theta)**5
end select
case(6)
 select case(m)
   case(-6)
   Ylm=(4.830841135800663d-1*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(1.673452458100098d0*Cos(theta)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(3.56781262853998d-1*(-1.d0 + 1.1d1*Cos(theta)**2)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(6.513904858677158d-1*Cos(theta)*(-3.d0 + 1.1d1*Cos(theta)**2)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(3.256952429338579d-1*(1.d0 - 1.7999999999999998d1*Cos(theta)**2 + 3.3000000000000003d1*Cos(theta)**4)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(4.119755163011408d-1*Cos(theta)*(5.d0 - 3.d1*Cos(theta)**2 + 3.3000000000000003d1*Cos(theta)**4)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=6.356920226762842d-2*(-5.d0 + 1.05d2*Cos(theta)**2 - 3.15d2*Cos(theta)**4 + 2.31d2*Cos(theta)**6)
   case(1)
   Ylm=-4.119755163011408d-1*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(5.d0 - 3.d1*Cos(theta)**2 + 3.3000000000000003d1*Cos(theta)**4)*Sin(theta)
   case(2)
   Ylm=3.256952429338579d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*(1.d0 - 1.7999999999999998d1*Cos(theta)**2 + 3.3000000000000003d1*Cos(theta)**4)*Sin(theta)**2
   case(3)
   Ylm=-6.513904858677158d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(-3.d0 + 1.1d1*Cos(theta)**2)*Sin(theta)**3
   case(4)
   Ylm=3.56781262853998d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*(-1.d0 + 1.1d1*Cos(theta)**2)*Sin(theta)**4
   case(5)
   Ylm=-1.673452458100098d0*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*Sin(theta)**5
   case(6)
   Ylm=4.830841135800663d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*Sin(theta)**6
end select
case(7)
 select case(m)
   case(-7)
   Ylm=(5.000395635705507d-1*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(1.870976726712969d0*Cos(theta)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(3.6692872457643775d-1*(-1.d0 + 1.3d1*Cos(theta)**2)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(7.338574491528755d-1*Cos(theta)*(-3.d0 + 1.3d1*Cos(theta)**2)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(1.1063317311124565d-1*(3.d0 - 6.6000000000000005d1*Cos(theta)**2 + 1.43d2*Cos(theta)**4)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.5645893386229404d-1*Cos(theta)*(1.5d1 - 1.1d2*Cos(theta)**2 + 1.43d2*Cos(theta)**4)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(6.387409227708014d-2*(-5.d0 + 1.35d2*Cos(theta)**2 - 4.95d2*Cos(theta)**4 + 4.29d2*Cos(theta)**6)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=6.828427691200495d-2*(-3.5d1*Cos(theta) + 3.15d2*Cos(theta)**3 - 6.93d2*Cos(theta)**5 + 4.29d2*Cos(theta)**7)
   case(1)
   Ylm=-6.387409227708014d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*(-5.d0 + 1.35d2*Cos(theta)**2 - 4.95d2*Cos(theta)**4 + 4.29d2*Cos(theta)**6)*Sin(theta)
   case(2)
   Ylm=1.5645893386229404d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(1.5d1 - 1.1d2*Cos(theta)**2 + 1.43d2*Cos(theta)**4)*Sin(theta)**2
   case(3)
   Ylm=-1.1063317311124565d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*(3.d0 - 6.6000000000000005d1*Cos(theta)**2 + 1.43d2*Cos(theta)**4)*Sin(theta)**3
   case(4)
   Ylm=7.338574491528755d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(-3.d0 + 1.3d1*Cos(theta)**2)*Sin(theta)**4
   case(5)
   Ylm=-3.6692872457643775d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*(-1.d0 + 1.3d1*Cos(theta)**2)*Sin(theta)**5
   case(6)
   Ylm=1.870976726712969d0*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*Sin(theta)**6
   case(7)
   Ylm=-5.000395635705507d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Sin(theta)**7
end select
case(8)
 select case(m)
   case(-8)
   Ylm=(5.154289843972844d-1*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(2.0617159375891374d0*Cos(theta)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(3.7641610872849456d-1*(-1.d0 + 1.5d1*Cos(theta)**2)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(2.439455195373073d0*Cos(theta)*(-1.d0 + 5.d0*Cos(theta)**2)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(3.382915688890246d-1*(1.d0 - 2.6d1*Cos(theta)**2 + 6.5d1*Cos(theta)**4)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(8.734650749797142d-1*Cos(theta)*(3.d0 - 2.6d1*Cos(theta)**2 + 3.9000000000000004d1*Cos(theta)**4)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(3.2254835519288303d-1*(-1.d0 + 3.3000000000000003d1*Cos(theta)**2 - 1.43d2*Cos(theta)**4 + 1.43d2*Cos(theta)**6)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(7.710380440405712d-2*Cos(theta)*(-3.5d1 + 3.85d2*Cos(theta)**2 - 1.001d3*Cos(theta)**4 + 7.1499999999999995d2*Cos(theta)**6)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=9.086770491564996d-3*(3.5d1 - 1.26d3*Cos(theta)**2 + 6.93d3*Cos(theta)**4 - 1.2012d4*Cos(theta)**6 + 6.435d3*Cos(theta)**8)
   case(1)
   Ylm=-7.710380440405712d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(-3.5d1 + 3.85d2*Cos(theta)**2 - 1.001d3*Cos(theta)**4 + 7.1499999999999995d2*Cos(theta)**6)*Sin(theta)
   case(2)
   Ylm=3.2254835519288303d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*(-1.d0 + 3.3000000000000003d1*Cos(theta)**2 - 1.43d2*Cos(theta)**4 + 1.43d2*Cos(theta)**6)*Sin(theta)**2
   case(3)
   Ylm=-8.734650749797142d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(3.d0 - 2.6d1*Cos(theta)**2 + 3.9000000000000004d1*Cos(theta)**4)*Sin(theta)**3
   case(4)
   Ylm=3.382915688890246d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*(1.d0 - 2.6d1*Cos(theta)**2 + 6.5d1*Cos(theta)**4)*Sin(theta)**4
   case(5)
   Ylm=-2.439455195373073d0*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(-1.d0 + 5.d0*Cos(theta)**2)*Sin(theta)**5
   case(6)
   Ylm=3.7641610872849456d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*(-1.d0 + 1.5d1*Cos(theta)**2)*Sin(theta)**6
   case(7)
   Ylm=-2.0617159375891374d0*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*Sin(theta)**7
   case(8)
   Ylm=5.154289843972844d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*Sin(theta)**8
end select
case(9)
 select case(m)
   case(-9)
   Ylm=(5.295529414924496d-1*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(2.246702855559565d0*Cos(theta)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(3.8530636096409983d-1*(-1.d0 + 1.7000000000000002d1*Cos(theta)**2)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(8.898269248923926d-1*Cos(theta)*(-3.d0 + 1.7000000000000002d1*Cos(theta)**2)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(3.446284861115194d-1*(1.d0 - 3.d1*Cos(theta)**2 + 8.5d1*Cos(theta)**4)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(2.883368783344621d0*Cos(theta)*(1.d0 - 1.d1*Cos(theta)**2 + 1.7000000000000002d1*Cos(theta)**4)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(3.264772254350559d-1*(-1.d0 + 3.9000000000000004d1*Cos(theta)**2 - 1.9500000000000002d2*Cos(theta)**4 + 2.21d2*Cos(theta)**6)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(4.274590280672302d-1*Cos(theta)*(-7.d0 + 9.1d1*Cos(theta)**2 - 2.7300000000000004d2*Cos(theta)**4 + 2.21d2*Cos(theta)**6)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(4.556728549830324d-2*(7.d0 - 3.08d2*Cos(theta)**2 + 2.002d3*Cos(theta)**4 - 4.004d3*Cos(theta)**6 + 2.431d3*Cos(theta)**8)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=9.606427264386591d-3*(3.15d2*Cos(theta) - 4.62d3*Cos(theta)**3 + 1.8018d4*Cos(theta)**5 - 2.5740000000000003d4*Cos(theta)**7 + 1.2155d4*Cos(theta)**9)
   case(1)
   Ylm=-4.556728549830324d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*(7.d0 - 3.08d2*Cos(theta)**2 + 2.002d3*Cos(theta)**4 - 4.004d3*Cos(theta)**6 + 2.431d3*Cos(theta)**8)*Sin(theta)
   case(2)
   Ylm=4.274590280672302d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(-7.d0 + 9.1d1*Cos(theta)**2 - 2.7300000000000004d2*Cos(theta)**4 + 2.21d2*Cos(theta)**6)*Sin(theta)**2
   case(3)
   Ylm=-3.264772254350559d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*(-1.d0 + 3.9000000000000004d1*Cos(theta)**2 - 1.9500000000000002d2*Cos(theta)**4 + 2.21d2*Cos(theta)**6)*Sin(theta)**3
   case(4)
   Ylm=2.883368783344621d0*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(1.d0 - 1.d1*Cos(theta)**2 + 1.7000000000000002d1*Cos(theta)**4)*Sin(theta)**4
   case(5)
   Ylm=-3.446284861115194d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*(1.d0 - 3.d1*Cos(theta)**2 + 8.5d1*Cos(theta)**4)*Sin(theta)**5
   case(6)
   Ylm=8.898269248923926d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*(-3.d0 + 1.7000000000000002d1*Cos(theta)**2)*Sin(theta)**6
   case(7)
   Ylm=-3.8530636096409983d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*(-1.d0 + 1.7000000000000002d1*Cos(theta)**2)*Sin(theta)**7
   case(8)
   Ylm=2.246702855559565d0*2.7182818284590446d0**((0.d0,8.d0)*phi)*Cos(theta)*Sin(theta)**8
   case(9)
   Ylm=-5.295529414924496d-1*2.7182818284590446d0**((0.d0,9.d0)*phi)*Sin(theta)**9
end select
case(10)
 select case(m)
   case(-10)
   Ylm=(5.426302919442215d-1*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(2.4267164388756717d0*Cos(theta)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(3.936653893957947d-1*(-1.d0 + 1.9d1*Cos(theta)**2)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(9.642793334137448d-1*Cos(theta)*(-3.d0 + 1.9d1*Cos(theta)**2)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(1.1693604541956055d-1*(3.d0 - 1.02d2*Cos(theta)**2 + 3.23d2*Cos(theta)**4)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(2.0918155726251224d-1*Cos(theta)*(1.5d1 - 1.7000000000000002d2*Cos(theta)**2 + 3.23d2*Cos(theta)**4)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(3.3074508272523757d-1*(-1.d0 + 4.5d1*Cos(theta)**2 - 2.55d2*Cos(theta)**4 + 3.23d2*Cos(theta)**6)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(4.677441816782422d-1*Cos(theta)*(-7.d0 + 1.05d2*Cos(theta)**2 - 3.57d2*Cos(theta)**4 + 3.23d2*Cos(theta)**6)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(4.586609057205472d-2*(7.d0 - 3.6399999999999997d2*Cos(theta)**2 + 2.7300000000000004d3*Cos(theta)**4 - 6.188000000000001d3*Cos(theta)**6 + 4.199d3*Cos(theta)**8)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(5.296159947690311d-2*Cos(theta)*(6.3d1 - 1.092d3*Cos(theta)**2 + 4.914d3*Cos(theta)**4 - 7.9559999999999995d3*Cos(theta)**6 + 4.199d3*Cos(theta)**8)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=5.049690376783604d-3*(-6.3d1 + 3.465d3*Cos(theta)**2 - 3.003d4*Cos(theta)**4 + 9.009d4*Cos(theta)**6 - 1.09395d5*Cos(theta)**8 + 4.6189d4*Cos(theta)**10)
   case(1)
   Ylm=-5.296159947690311d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(6.3d1 - 1.092d3*Cos(theta)**2 + 4.914d3*Cos(theta)**4 - 7.9559999999999995d3*Cos(theta)**6 + 4.199d3*Cos(theta)**8)*Sin(theta)
   case(2)
   Ylm=4.586609057205472d-2*2.7182818284590446d0**((0.d0,2.d0)*phi)*(7.d0 - 3.6399999999999997d2*Cos(theta)**2 + 2.7300000000000004d3*Cos(theta)**4 - 6.188000000000001d3*Cos(theta)**6 + 4.199d3*Cos(theta)**8)*Sin(theta)**2
   case(3)
   Ylm=-4.677441816782422d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(-7.d0 + 1.05d2*Cos(theta)**2 - 3.57d2*Cos(theta)**4 + 3.23d2*Cos(theta)**6)*Sin(theta)**3
   case(4)
   Ylm=3.3074508272523757d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*(-1.d0 + 4.5d1*Cos(theta)**2 - 2.55d2*Cos(theta)**4 + 3.23d2*Cos(theta)**6)*Sin(theta)**4
   case(5)
   Ylm=-2.0918155726251224d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(1.5d1 - 1.7000000000000002d2*Cos(theta)**2 + 3.23d2*Cos(theta)**4)*Sin(theta)**5
   case(6)
   Ylm=1.1693604541956055d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*(3.d0 - 1.02d2*Cos(theta)**2 + 3.23d2*Cos(theta)**4)*Sin(theta)**6
   case(7)
   Ylm=-9.642793334137448d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*(-3.d0 + 1.9d1*Cos(theta)**2)*Sin(theta)**7
   case(8)
   Ylm=3.936653893957947d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*(-1.d0 + 1.9d1*Cos(theta)**2)*Sin(theta)**8
   case(9)
   Ylm=-2.4267164388756717d0*2.7182818284590446d0**((0.d0,9.d0)*phi)*Cos(theta)*Sin(theta)**9
   case(10)
   Ylm=5.426302919442215d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*Sin(theta)**10
end select
case(11)
 select case(m)
   case(-11)
   Ylm=(5.548257538066194d-1*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(2.6023634596104817d0*Cos(theta)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(4.015533996368363d-1*(-1.d0 + 2.1d1*Cos(theta)**2)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(3.110419258812877d0*Cos(theta)*(-1.d0 + 7.d0*Cos(theta)**2)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(3.567895584528425d-1*(1.d0 - 3.8d1*Cos(theta)**2 + 1.33d2*Cos(theta)**4)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(2.2565353001535278d-1*Cos(theta)*(1.5d1 - 1.9d2*Cos(theta)**2 + 3.99d2*Cos(theta)**4)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(6.702908649261444d-2*(-5.d0 + 2.55d2*Cos(theta)**2 - 1.615d3*Cos(theta)**4 + 2.261d3*Cos(theta)**6)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(7.093691738691861d-1*Cos(theta)*(-5.d0 + 8.5d1*Cos(theta)**2 - 3.23d2*Cos(theta)**4 + 3.23d2*Cos(theta)**6)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(3.2378124843913114d-1*(1.d0 - 6.d1*Cos(theta)**2 + 5.1d2*Cos(theta)**4 - 1.292d3*Cos(theta)**6 + 9.69d2*Cos(theta)**8)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.7306835713159483d-1*Cos(theta)*(2.1d1 - 4.2d2*Cos(theta)**2 + 2.142d3*Cos(theta)**4 - 3.876d3*Cos(theta)**6 + 2.261d3*Cos(theta)**8)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(1.517909905105581d-2*(-2.1d1 + 1.3650000000000002d3*Cos(theta)**2 - 1.3650000000000002d4*Cos(theta)**4 + 4.641d4*Cos(theta)**6 - 6.298500000000001d4*Cos(theta)**8 + 2.9393000000000002d4*Cos(theta)**10)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=5.284683964654306d-3*(-6.93d2*Cos(theta) + 1.5015d4*Cos(theta)**3 - 9.009d4*Cos(theta)**5 + 2.1879d5*Cos(theta)**7 - 2.30945d5*Cos(theta)**9 + 8.8179d4*Cos(theta)**11)
   case(1)
   Ylm=-1.517909905105581d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*(-2.1d1 + 1.3650000000000002d3*Cos(theta)**2 - 1.3650000000000002d4*Cos(theta)**4 + 4.641d4*Cos(theta)**6 - 6.298500000000001d4*Cos(theta)**8 + 2.9393000000000002d4*Cos(theta)**10)*Sin(theta)
   case(2)
   Ylm=1.7306835713159483d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(2.1d1 - 4.2d2*Cos(theta)**2 + 2.142d3*Cos(theta)**4 - 3.876d3*Cos(theta)**6 + 2.261d3*Cos(theta)**8)*Sin(theta)**2
   case(3)
   Ylm=-3.2378124843913114d-1*2.7182818284590446d0**((0.d0,3.d0)*phi)*(1.d0 - 6.d1*Cos(theta)**2 + 5.1d2*Cos(theta)**4 - 1.292d3*Cos(theta)**6 + 9.69d2*Cos(theta)**8)*Sin(theta)**3
   case(4)
   Ylm=7.093691738691861d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(-5.d0 + 8.5d1*Cos(theta)**2 - 3.23d2*Cos(theta)**4 + 3.23d2*Cos(theta)**6)*Sin(theta)**4
   case(5)
   Ylm=-6.702908649261444d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*(-5.d0 + 2.55d2*Cos(theta)**2 - 1.615d3*Cos(theta)**4 + 2.261d3*Cos(theta)**6)*Sin(theta)**5
   case(6)
   Ylm=2.2565353001535278d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*(1.5d1 - 1.9d2*Cos(theta)**2 + 3.99d2*Cos(theta)**4)*Sin(theta)**6
   case(7)
   Ylm=-3.567895584528425d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*(1.d0 - 3.8d1*Cos(theta)**2 + 1.33d2*Cos(theta)**4)*Sin(theta)**7
   case(8)
   Ylm=3.110419258812877d0*2.7182818284590446d0**((0.d0,8.d0)*phi)*Cos(theta)*(-1.d0 + 7.d0*Cos(theta)**2)*Sin(theta)**8
   case(9)
   Ylm=-4.015533996368363d-1*2.7182818284590446d0**((0.d0,9.d0)*phi)*(-1.d0 + 2.1d1*Cos(theta)**2)*Sin(theta)**9
   case(10)
   Ylm=2.6023634596104817d0*2.7182818284590446d0**((0.d0,1.d1)*phi)*Cos(theta)*Sin(theta)**10
   case(11)
   Ylm=-5.548257538066194d-1*2.7182818284590446d0**((0.d0,1.1d1)*phi)*Sin(theta)**11
end select
case(12)
 select case(m)
   case(-12)
   Ylm=(5.6626666374219115d-1*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(2.774128769033097d0*Cos(theta)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(4.090229723318166d-1*(-1.d0 + 2.3000000000000003d1*Cos(theta)**2)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(1.1076394452006766d0*Cos(theta)*(-3.d0 + 2.3000000000000003d1*Cos(theta)**2)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(3.625601143107851d-1*(1.d0 - 4.2d1*Cos(theta)**2 + 1.61d2*Cos(theta)**4)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(7.251202286215702d-1*Cos(theta)*(5.d0 - 7.d1*Cos(theta)**2 + 1.61d2*Cos(theta)**4)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(6.791373178178367d-2*(-5.d0 + 2.8499999999999996d2*Cos(theta)**2 - 1.995d3*Cos(theta)**4 + 3.059d3*Cos(theta)**6)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(7.623297485540853d-1*Cos(theta)*(-5.d0 + 9.5d1*Cos(theta)**2 - 3.99d2*Cos(theta)**4 + 4.37d2*Cos(theta)**6)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(6.536923664453508d-2*(5.d0 - 3.4000000000000004d2*Cos(theta)**2 + 3.23d3*Cos(theta)**4 - 9.044d3*Cos(theta)**6 + 7.429d3*Cos(theta)**8)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(8.715898219271343d-2*Cos(theta)*(4.5d1 - 1.02d3*Cos(theta)**2 + 5.814d3*Cos(theta)**4 - 1.1627999999999998d4*Cos(theta)**6 + 7.429d3*Cos(theta)**8)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.0674751643623661d-1*(-3.d0 + 2.25d2*Cos(theta)**2 - 2.55d3*Cos(theta)**4 + 9.69d3*Cos(theta)**6 - 1.4535d4*Cos(theta)**8 + 7.429d3*Cos(theta)**10)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(1.7203920019399237d-2*Cos(theta)*(-2.31d2 + 5.775d3*Cos(theta)**2 - 3.927d4*Cos(theta)**4 + 1.0659d5*Cos(theta)**6 - 1.24355d5*Cos(theta)**8 + 5.2003d4*Cos(theta)**10)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=1.3774159754583892d-3*(2.31d2 - 1.8018d4*Cos(theta)**2 + 2.25225d5*Cos(theta)**4 - 1.02102d6*Cos(theta)**6 + 2.078505d6*Cos(theta)**8 - 1.939938d6*Cos(theta)**10 + 6.760389999999999d5*Cos(theta)**12)
   case(1)
   Ylm=-1.7203920019399237d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(-2.31d2 + 5.775d3*Cos(theta)**2 - 3.927d4*Cos(theta)**4 + 1.0659d5*Cos(theta)**6 - 1.24355d5*Cos(theta)**8 + 5.2003d4*Cos(theta)**10)*Sin(theta)
   case(2)
   Ylm=1.0674751643623661d-1*2.7182818284590446d0**((0.d0,2.d0)*phi)*(-3.d0 + 2.25d2*Cos(theta)**2 - 2.55d3*Cos(theta)**4 + 9.69d3*Cos(theta)**6 - 1.4535d4*Cos(theta)**8 + 7.429d3*Cos(theta)**10)*Sin(theta)**2
   case(3)
   Ylm=-8.715898219271343d-2*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(4.5d1 - 1.02d3*Cos(theta)**2 + 5.814d3*Cos(theta)**4 - 1.1627999999999998d4*Cos(theta)**6 + 7.429d3*Cos(theta)**8)*Sin(theta)**3
   case(4)
   Ylm=6.536923664453508d-2*2.7182818284590446d0**((0.d0,4.d0)*phi)*(5.d0 - 3.4000000000000004d2*Cos(theta)**2 + 3.23d3*Cos(theta)**4 - 9.044d3*Cos(theta)**6 + 7.429d3*Cos(theta)**8)*Sin(theta)**4
   case(5)
   Ylm=-7.623297485540853d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(-5.d0 + 9.5d1*Cos(theta)**2 - 3.99d2*Cos(theta)**4 + 4.37d2*Cos(theta)**6)*Sin(theta)**5
   case(6)
   Ylm=6.791373178178367d-2*2.7182818284590446d0**((0.d0,6.d0)*phi)*(-5.d0 + 2.8499999999999996d2*Cos(theta)**2 - 1.995d3*Cos(theta)**4 + 3.059d3*Cos(theta)**6)*Sin(theta)**6
   case(7)
   Ylm=-7.251202286215702d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*(5.d0 - 7.d1*Cos(theta)**2 + 1.61d2*Cos(theta)**4)*Sin(theta)**7
   case(8)
   Ylm=3.625601143107851d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*(1.d0 - 4.2d1*Cos(theta)**2 + 1.61d2*Cos(theta)**4)*Sin(theta)**8
   case(9)
   Ylm=-1.1076394452006766d0*2.7182818284590446d0**((0.d0,9.d0)*phi)*Cos(theta)*(-3.d0 + 2.3000000000000003d1*Cos(theta)**2)*Sin(theta)**9
   case(10)
   Ylm=4.090229723318166d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*(-1.d0 + 2.3000000000000003d1*Cos(theta)**2)*Sin(theta)**10
   case(11)
   Ylm=-2.774128769033097d0*2.7182818284590446d0**((0.d0,1.1d1)*phi)*Cos(theta)*Sin(theta)**11
   case(12)
   Ylm=5.6626666374219115d-1*2.7182818284590446d0**((0.d0,1.2d1)*phi)*Sin(theta)**12
end select
case(13)
 select case(m)
   case(-13)
   Ylm=(5.770536647012669d-1*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(2.942407896701988d0*Cos(theta)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(4.161193153549645d-1*(-1.d0 + 2.5d1*Cos(theta)**2)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(1.1769631586807954d0*Cos(theta)*(-3.d0 + 2.5d1*Cos(theta)**2)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(1.2270689169954498d-1*(3.d0 - 1.3800000000000001d2*Cos(theta)**2 + 5.75d2*Cos(theta)**4)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(1.286960737459393d0*Cos(theta)*(3.d0 - 4.6000000000000005d1*Cos(theta)**2 + 1.1500000000000001d2*Cos(theta)**4)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(3.439547249859269d-1*(-1.d0 + 6.3d1*Cos(theta)**2 - 4.83d2*Cos(theta)**4 + 8.05d2*Cos(theta)**6)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(8.139454379163322d-1*Cos(theta)*(-5.d0 + 1.05d2*Cos(theta)**2 - 4.83d2*Cos(theta)**4 + 5.75d2*Cos(theta)**6)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(6.601969283084413d-2*(5.d0 - 3.8d2*Cos(theta)**2 + 3.99d3*Cos(theta)**4 - 1.2236d4*Cos(theta)**6 + 1.0925d4*Cos(theta)**8)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(9.336594498508555d-2*Cos(theta)*(4.5d1 - 1.1400000000000001d3*Cos(theta)**2 + 7.1819999999999995d3*Cos(theta)**4 - 1.5732d4*Cos(theta)**6 + 1.0925d4*Cos(theta)**8)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(3.580420547710517d-2*(-9.d0 + 7.65d2*Cos(theta)**2 - 9.69d3*Cos(theta)**4 + 4.0698d4*Cos(theta)**6 - 6.686100000000001d4*Cos(theta)**8 + 3.7145d4*Cos(theta)**10)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(4.318149653976205d-2*Cos(theta)*(-9.9d1 + 2.805d3*Cos(theta)**2 - 2.1318d4*Cos(theta)**4 + 6.3954d4*Cos(theta)**6 - 8.171899999999999d4*Cos(theta)**8 + 3.7145d4*Cos(theta)**10)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(9.655676163307987d-3*(3.3000000000000003d1 - 2.9699999999999998d3*Cos(theta)**2 + 4.2075000000000005d4*Cos(theta)**4 - 2.1318d5*Cos(theta)**6 + 4.79655d5*Cos(theta)**8 - 4.9031400000000005d5*Cos(theta)**10 + 1.85725d5*Cos(theta)**12)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=1.431452671590586d-3*(3.003d3*Cos(theta) - 9.009d4*Cos(theta)**3 + 7.65765d5*Cos(theta)**5 - 2.77134d6*Cos(theta)**7 + 4.849845d6*Cos(theta)**9 - 4.056234d6*Cos(theta)**11 + 1.300075d6*Cos(theta)**13)
   case(1)
   Ylm=-9.655676163307987d-3*2.7182818284590446d0**((0.d0,1.d0)*phi)*(3.3000000000000003d1 - 2.9699999999999998d3*Cos(theta)**2 + 4.2075000000000005d4*Cos(theta)**4 - 2.1318d5*Cos(theta)**6 + 4.79655d5*Cos(theta)**8 - 4.9031400000000005d5*Cos(theta)**10 + 1.85725d5*Cos(theta)**12)*Sin(theta)
   case(2)
   Ylm=4.318149653976205d-2*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(-9.9d1 + 2.805d3*Cos(theta)**2 - 2.1318d4*Cos(theta)**4 + 6.3954d4*Cos(theta)**6 - 8.171899999999999d4*Cos(theta)**8 + 3.7145d4*Cos(theta)**10)*Sin(theta)**2
   case(3)
   Ylm=-3.580420547710517d-2*2.7182818284590446d0**((0.d0,3.d0)*phi)*(-9.d0 + 7.65d2*Cos(theta)**2 - 9.69d3*Cos(theta)**4 + 4.0698d4*Cos(theta)**6 - 6.686100000000001d4*Cos(theta)**8 + 3.7145d4*Cos(theta)**10)*Sin(theta)**3
   case(4)
   Ylm=9.336594498508555d-2*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(4.5d1 - 1.1400000000000001d3*Cos(theta)**2 + 7.1819999999999995d3*Cos(theta)**4 - 1.5732d4*Cos(theta)**6 + 1.0925d4*Cos(theta)**8)*Sin(theta)**4
   case(5)
   Ylm=-6.601969283084413d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*(5.d0 - 3.8d2*Cos(theta)**2 + 3.99d3*Cos(theta)**4 - 1.2236d4*Cos(theta)**6 + 1.0925d4*Cos(theta)**8)*Sin(theta)**5
   case(6)
   Ylm=8.139454379163322d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*(-5.d0 + 1.05d2*Cos(theta)**2 - 4.83d2*Cos(theta)**4 + 5.75d2*Cos(theta)**6)*Sin(theta)**6
   case(7)
   Ylm=-3.439547249859269d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*(-1.d0 + 6.3d1*Cos(theta)**2 - 4.83d2*Cos(theta)**4 + 8.05d2*Cos(theta)**6)*Sin(theta)**7
   case(8)
   Ylm=1.286960737459393d0*2.7182818284590446d0**((0.d0,8.d0)*phi)*Cos(theta)*(3.d0 - 4.6000000000000005d1*Cos(theta)**2 + 1.1500000000000001d2*Cos(theta)**4)*Sin(theta)**8
   case(9)
   Ylm=-1.2270689169954498d-1*2.7182818284590446d0**((0.d0,9.d0)*phi)*(3.d0 - 1.3800000000000001d2*Cos(theta)**2 + 5.75d2*Cos(theta)**4)*Sin(theta)**9
   case(10)
   Ylm=1.1769631586807954d0*2.7182818284590446d0**((0.d0,1.d1)*phi)*Cos(theta)*(-3.d0 + 2.5d1*Cos(theta)**2)*Sin(theta)**10
   case(11)
   Ylm=-4.161193153549645d-1*2.7182818284590446d0**((0.d0,1.1d1)*phi)*(-1.d0 + 2.5d1*Cos(theta)**2)*Sin(theta)**11
   case(12)
   Ylm=2.942407896701988d0*2.7182818284590446d0**((0.d0,1.2d1)*phi)*Cos(theta)*Sin(theta)**12
   case(13)
   Ylm=-5.770536647012669d-1*2.7182818284590446d0**((0.d0,1.3d1)*phi)*Sin(theta)**13
end select
case(14)
 select case(m)
   case(-14)
   Ylm=(5.87267796860102d-1*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(3.107529086977258d0*Cos(theta)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(4.2288114577506475d-1*(-1.d0 + 2.7d1*Cos(theta)**2)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(3.734785154364099d0*Cos(theta)*(-1.d0 + 9.d0*Cos(theta)**2)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(3.734785154364099d-1*(1.d0 - 5.d1*Cos(theta)**2 + 2.25d2*Cos(theta)**4)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(1.3637507176537538d0*Cos(theta)*(3.d0 - 5.d1*Cos(theta)**2 + 1.35d2*Cos(theta)**4)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(3.4827051141890646d-1*(-1.d0 + 6.8999999999999995d1*Cos(theta)**2 - 5.75d2*Cos(theta)**4 + 1.035d3*Cos(theta)**6)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(6.174176267472802d-1*Cos(theta)*(-7.d0 + 1.61d2*Cos(theta)**2 - 8.05d2*Cos(theta)**4 + 1.035d3*Cos(theta)**6)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(3.334436284646243d-1*(1.d0 - 8.4d1*Cos(theta)**2 + 9.66d2*Cos(theta)**4 - 3.22d3*Cos(theta)**6 + 3.105d3*Cos(theta)**8)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(8.947231438933006d-1*Cos(theta)*(5.d0 - 1.4000000000000001d2*Cos(theta)**2 + 9.66d2*Cos(theta)**4 - 2.3000000000000003d3*Cos(theta)**6 + 1.7249999999999999d3*Cos(theta)**8)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(3.2455019565917604d-1*(-1.d0 + 9.5d1*Cos(theta)**2 - 1.33d3*Cos(theta)**4 + 6.118d3*Cos(theta)**6 - 1.0925d4*Cos(theta)**8 + 6.555d3*Cos(theta)**10)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(4.612955613859324d-2*Cos(theta)*(-9.9d1 + 3.135d3*Cos(theta)**2 - 2.6334d4*Cos(theta)**4 + 8.6526d4*Cos(theta)**6 - 1.20175d5*Cos(theta)**8 + 5.8995d4*Cos(theta)**10)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(9.689144811888621d-3*(3.3000000000000003d1 - 3.366d3*Cos(theta)**2 + 5.3295d4*Cos(theta)**4 - 2.98452d5*Cos(theta)**6 + 7.35471d5*Cos(theta)**8 - 8.171899999999999d5*Cos(theta)**10 + 3.3430500000000003d5*Cos(theta)**12)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(1.0749141056818556d-2*Cos(theta)*(4.29d2 - 1.4586d4*Cos(theta)**2 + 1.38567d5*Cos(theta)**4 - 5.54268d5*Cos(theta)**6 + 1.062347d6*Cos(theta)**8 - 9.6577d5*Cos(theta)**10 + 3.3430500000000003d5*Cos(theta)**12)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=7.417612035823363d-4*(-4.29d2 + 4.5045d4*Cos(theta)**2 - 7.65765d5*Cos(theta)**4 + 4.849845d6*Cos(theta)**6 - 1.4549535d7*Cos(theta)**8 + 2.2309286999999998d7*Cos(theta)**10 - 1.6900974999999998d7*Cos(theta)**12 + 5.014575d6*Cos(theta)**14)
   case(1)
   Ylm=-1.0749141056818556d-2*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(4.29d2 - 1.4586d4*Cos(theta)**2 + 1.38567d5*Cos(theta)**4 - 5.54268d5*Cos(theta)**6 + 1.062347d6*Cos(theta)**8 - 9.6577d5*Cos(theta)**10 + 3.3430500000000003d5*Cos(theta)**12)*Sin(theta)
   case(2)
   Ylm=9.689144811888621d-3*2.7182818284590446d0**((0.d0,2.d0)*phi)*(3.3000000000000003d1 - 3.366d3*Cos(theta)**2 + 5.3295d4*Cos(theta)**4 - 2.98452d5*Cos(theta)**6 + 7.35471d5*Cos(theta)**8 - 8.171899999999999d5*Cos(theta)**10 + 3.3430500000000003d5*Cos(theta)**12)*Sin(theta)**2
   case(3)
   Ylm=-4.612955613859324d-2*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(-9.9d1 + 3.135d3*Cos(theta)**2 - 2.6334d4*Cos(theta)**4 + 8.6526d4*Cos(theta)**6 - 1.20175d5*Cos(theta)**8 + 5.8995d4*Cos(theta)**10)*Sin(theta)**3
   case(4)
   Ylm=3.2455019565917604d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*(-1.d0 + 9.5d1*Cos(theta)**2 - 1.33d3*Cos(theta)**4 + 6.118d3*Cos(theta)**6 - 1.0925d4*Cos(theta)**8 + 6.555d3*Cos(theta)**10)*Sin(theta)**4
   case(5)
   Ylm=-8.947231438933006d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(5.d0 - 1.4000000000000001d2*Cos(theta)**2 + 9.66d2*Cos(theta)**4 - 2.3000000000000003d3*Cos(theta)**6 + 1.7249999999999999d3*Cos(theta)**8)*Sin(theta)**5
   case(6)
   Ylm=3.334436284646243d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*(1.d0 - 8.4d1*Cos(theta)**2 + 9.66d2*Cos(theta)**4 - 3.22d3*Cos(theta)**6 + 3.105d3*Cos(theta)**8)*Sin(theta)**6
   case(7)
   Ylm=-6.174176267472802d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*(-7.d0 + 1.61d2*Cos(theta)**2 - 8.05d2*Cos(theta)**4 + 1.035d3*Cos(theta)**6)*Sin(theta)**7
   case(8)
   Ylm=3.4827051141890646d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*(-1.d0 + 6.8999999999999995d1*Cos(theta)**2 - 5.75d2*Cos(theta)**4 + 1.035d3*Cos(theta)**6)*Sin(theta)**8
   case(9)
   Ylm=-1.3637507176537538d0*2.7182818284590446d0**((0.d0,9.d0)*phi)*Cos(theta)*(3.d0 - 5.d1*Cos(theta)**2 + 1.35d2*Cos(theta)**4)*Sin(theta)**9
   case(10)
   Ylm=3.734785154364099d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*(1.d0 - 5.d1*Cos(theta)**2 + 2.25d2*Cos(theta)**4)*Sin(theta)**10
   case(11)
   Ylm=-3.734785154364099d0*2.7182818284590446d0**((0.d0,1.1d1)*phi)*Cos(theta)*(-1.d0 + 9.d0*Cos(theta)**2)*Sin(theta)**11
   case(12)
   Ylm=4.2288114577506475d-1*2.7182818284590446d0**((0.d0,1.2d1)*phi)*(-1.d0 + 2.7d1*Cos(theta)**2)*Sin(theta)**12
   case(13)
   Ylm=-3.107529086977258d0*2.7182818284590446d0**((0.d0,1.3d1)*phi)*Cos(theta)*Sin(theta)**13
   case(14)
   Ylm=5.87267796860102d-1*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*Sin(theta)**14
end select
case(15)
 select case(m)
   case(-15)
   Ylm=(5.969753602424045d-1*Sin(theta)**15)/2.7182818284590446d0**((0.d0,1.5d1)*phi)
   case(-14)
   Ylm=(3.2697687107953763d0*Cos(theta)*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(4.2934166569087475d-1*(-1.d0 + 2.9d1*Cos(theta)**2)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(1.3116604546845723d0*Cos(theta)*(-3.d0 + 2.9d1*Cos(theta)**2)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(3.786437582987623d-1*(1.d0 - 5.4d1*Cos(theta)**2 + 2.6100000000000003d2*Cos(theta)**4)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(8.634406161588531d-1*Cos(theta)*(5.d0 - 9.d1*Cos(theta)**2 + 2.6100000000000003d2*Cos(theta)**4)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(3.5249815546391634d-1*(-1.d0 + 7.5d1*Cos(theta)**2 - 6.75d2*Cos(theta)**4 + 1.3050000000000002d3*Cos(theta)**6)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(6.526997549224868d-1*Cos(theta)*(-7.d0 + 1.75d2*Cos(theta)**2 - 9.45d2*Cos(theta)**4 + 1.3050000000000002d3*Cos(theta)**6)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(4.81176643237967d-2*(7.d0 - 6.44d2*Cos(theta)**2 + 8.05d3*Cos(theta)**4 - 2.898d4*Cos(theta)**6 + 3.0014999999999996d4*Cos(theta)**8)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(2.256918510702296d-1*Cos(theta)*(2.1d1 - 6.44d2*Cos(theta)**2 + 4.83d3*Cos(theta)**4 - 1.242d4*Cos(theta)**6 + 1.0005d4*Cos(theta)**8)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(3.2705856424035757d-1*(-1.d0 + 1.05d2*Cos(theta)**2 - 1.61d3*Cos(theta)**4 + 8.05d3*Cos(theta)**6 - 1.5525d4*Cos(theta)**8 + 1.0005d4*Cos(theta)**10)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(4.41005678056549d-1*Cos(theta)*(-1.1d1 + 3.85d2*Cos(theta)**2 - 3.5420000000000003d3*Cos(theta)**4 + 1.2650000000000001d4*Cos(theta)**6 - 1.8975d4*Cos(theta)**8 + 1.0005d4*Cos(theta)**10)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(2.9206314511765434d-2*(1.1d1 - 1.254d3*Cos(theta)**2 + 2.1945d4*Cos(theta)**4 - 1.3459599999999998d5*Cos(theta)**6 + 3.60525d5*Cos(theta)**8 - 4.3263d5*Cos(theta)**10 + 1.90095d5*Cos(theta)**12)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.145565904735539d-2*Cos(theta)*(4.29d2 - 1.6301999999999999d4*Cos(theta)**2 + 1.7117099999999998d5*Cos(theta)**4 - 7.49892d5*Cos(theta)**6 + 1.5622749999999999d6*Cos(theta)**8 - 1.5338699999999998d6*Cos(theta)**10 + 5.702850000000001d5*Cos(theta)**12)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(7.425600850239804d-4*(-4.29d2 + 5.1051d4*Cos(theta)**2 - 9.69969d5*Cos(theta)**4 + 6.789783d6*Cos(theta)**6 - 2.2309286999999998d7*Cos(theta)**8 + 3.7182145d7*Cos(theta)**10 - 3.0421755d7*Cos(theta)**12 + 9.694844999999999d6*Cos(theta)**14)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=7.669127580949972d-4*(-6.435d3*Cos(theta) + 2.55255d5*Cos(theta)**3 - 2.909907d6*Cos(theta)**5 + 1.4549535d7*Cos(theta)**7 - 3.7182145d7*Cos(theta)**9 + 5.0702925d7*Cos(theta)**11 - 3.5102024999999997d7*Cos(theta)**13 + 9.694844999999999d6*Cos(theta)**15)
   case(1)
   Ylm=-7.425600850239804d-4*2.7182818284590446d0**((0.d0,1.d0)*phi)*(-4.29d2 + 5.1051d4*Cos(theta)**2 - 9.69969d5*Cos(theta)**4 + 6.789783d6*Cos(theta)**6 - 2.2309286999999998d7*Cos(theta)**8 + 3.7182145d7*Cos(theta)**10 - 3.0421755d7*Cos(theta)**12 + 9.694844999999999d6*Cos(theta)**14)*Sin(theta)
   case(2)
   Ylm=1.145565904735539d-2*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(4.29d2 - 1.6301999999999999d4*Cos(theta)**2 + 1.7117099999999998d5*Cos(theta)**4 - 7.49892d5*Cos(theta)**6 + 1.5622749999999999d6*Cos(theta)**8 - 1.5338699999999998d6*Cos(theta)**10 + 5.702850000000001d5*Cos(theta)**12)*Sin(theta)**2
   case(3)
   Ylm=-2.9206314511765434d-2*2.7182818284590446d0**((0.d0,3.d0)*phi)*(1.1d1 - 1.254d3*Cos(theta)**2 + 2.1945d4*Cos(theta)**4 - 1.3459599999999998d5*Cos(theta)**6 + 3.60525d5*Cos(theta)**8 - 4.3263d5*Cos(theta)**10 + 1.90095d5*Cos(theta)**12)*Sin(theta)**3
   case(4)
   Ylm=4.41005678056549d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(-1.1d1 + 3.85d2*Cos(theta)**2 - 3.5420000000000003d3*Cos(theta)**4 + 1.2650000000000001d4*Cos(theta)**6 - 1.8975d4*Cos(theta)**8 + 1.0005d4*Cos(theta)**10)*Sin(theta)**4
   case(5)
   Ylm=-3.2705856424035757d-1*2.7182818284590446d0**((0.d0,5.d0)*phi)*(-1.d0 + 1.05d2*Cos(theta)**2 - 1.61d3*Cos(theta)**4 + 8.05d3*Cos(theta)**6 - 1.5525d4*Cos(theta)**8 + 1.0005d4*Cos(theta)**10)*Sin(theta)**5
   case(6)
   Ylm=2.256918510702296d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*(2.1d1 - 6.44d2*Cos(theta)**2 + 4.83d3*Cos(theta)**4 - 1.242d4*Cos(theta)**6 + 1.0005d4*Cos(theta)**8)*Sin(theta)**6
   case(7)
   Ylm=-4.81176643237967d-2*2.7182818284590446d0**((0.d0,7.d0)*phi)*(7.d0 - 6.44d2*Cos(theta)**2 + 8.05d3*Cos(theta)**4 - 2.898d4*Cos(theta)**6 + 3.0014999999999996d4*Cos(theta)**8)*Sin(theta)**7
   case(8)
   Ylm=6.526997549224868d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*Cos(theta)*(-7.d0 + 1.75d2*Cos(theta)**2 - 9.45d2*Cos(theta)**4 + 1.3050000000000002d3*Cos(theta)**6)*Sin(theta)**8
   case(9)
   Ylm=-3.5249815546391634d-1*2.7182818284590446d0**((0.d0,9.d0)*phi)*(-1.d0 + 7.5d1*Cos(theta)**2 - 6.75d2*Cos(theta)**4 + 1.3050000000000002d3*Cos(theta)**6)*Sin(theta)**9
   case(10)
   Ylm=8.634406161588531d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*Cos(theta)*(5.d0 - 9.d1*Cos(theta)**2 + 2.6100000000000003d2*Cos(theta)**4)*Sin(theta)**10
   case(11)
   Ylm=-3.786437582987623d-1*2.7182818284590446d0**((0.d0,1.1d1)*phi)*(1.d0 - 5.4d1*Cos(theta)**2 + 2.6100000000000003d2*Cos(theta)**4)*Sin(theta)**11
   case(12)
   Ylm=1.3116604546845723d0*2.7182818284590446d0**((0.d0,1.2d1)*phi)*Cos(theta)*(-3.d0 + 2.9d1*Cos(theta)**2)*Sin(theta)**12
   case(13)
   Ylm=-4.2934166569087475d-1*2.7182818284590446d0**((0.d0,1.3d1)*phi)*(-1.d0 + 2.9d1*Cos(theta)**2)*Sin(theta)**13
   case(14)
   Ylm=3.2697687107953763d0*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*Cos(theta)*Sin(theta)**14
   case(15)
   Ylm=-5.969753602424045d-1*2.7182818284590446d0**((0.d0,1.5d1)*phi)*Sin(theta)**15
end select
case(16)
 select case(m)
   case(-16)
   Ylm=(6.062313441538353d-1*Sin(theta)**16)/2.7182818284590446d0**((0.d0,1.6d1)*phi)
   case(-15)
   Ylm=(3.4293623553521d0*Cos(theta)*Sin(theta)**15)/2.7182818284590446d0**((0.d0,1.5d1)*phi)
   case(-14)
   Ylm=(4.355294546593892d-1*(-1.d0 + 3.1d1*Cos(theta)**2)*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(1.3772650648147036d0*Cos(theta)*(-3.d0 + 3.1d1*Cos(theta)**2)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(1.2787585098510281d-1*(3.d0 - 1.7399999999999998d2*Cos(theta)**2 + 8.99d2*Cos(theta)**4)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(3.0260949470385032d-1*Cos(theta)*(1.5d1 - 2.9d2*Cos(theta)**2 + 8.99d2*Cos(theta)**4)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(7.13257419188424d-2*(-5.d0 + 4.050000000000001d2*Cos(theta)**2 - 3.915d3*Cos(theta)**4 + 8.091000000000001d3*Cos(theta)**6)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(1.3746240938998877d-1*Cos(theta)*(-3.5d1 + 9.45d2*Cos(theta)**2 - 5.481d3*Cos(theta)**4 + 8.091000000000001d3*Cos(theta)**6)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(4.860030091895121d-2*(7.d0 - 7.d2*Cos(theta)**2 + 9.45d3*Cos(theta)**4 - 3.654d4*Cos(theta)**6 + 4.0455000000000005d4*Cos(theta)**8)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(2.380918771942937d-1*Cos(theta)*(2.1d1 - 7.d2*Cos(theta)**2 + 5.67d3*Cos(theta)**4 - 1.5659999999999998d4*Cos(theta)**6 + 1.3485d4*Cos(theta)**8)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(1.5699313469606495d-2*(-2.1d1 + 2.415d3*Cos(theta)**2 - 4.025d4*Cos(theta)**4 + 2.1734999999999998d5*Cos(theta)**6 - 4.50225d5*Cos(theta)**8 + 3.10155d5*Cos(theta)**10)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(2.2202182028664117d-2*Cos(theta)*(-2.31d2 + 8.855d3*Cos(theta)**2 - 8.855d4*Cos(theta)**4 + 3.4155d5*Cos(theta)**6 - 5.50275d5*Cos(theta)**8 + 3.10155d5*Cos(theta)**10)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(2.937072610541639d-2*(1.1d1 - 1.3860000000000001d3*Cos(theta)**2 + 2.6565d4*Cos(theta)**4 - 1.7710000000000001d5*Cos(theta)**6 + 5.1232500000000005d5*Cos(theta)**8 - 6.6033d5*Cos(theta)**10 + 3.10155d5*Cos(theta)**12)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(3.642990217177658d-2*Cos(theta)*(1.43d2 - 6.006d3*Cos(theta)**2 + 6.9069d4*Cos(theta)**4 - 3.289d5*Cos(theta)**6 + 7.400250000000001d5*Cos(theta)**8 - 7.8039000000000005d5*Cos(theta)**10 + 3.10155d5*Cos(theta)**12)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(2.233660615510501d-3*(-1.43d2 + 1.9019d4*Cos(theta)**2 - 3.99399d5*Cos(theta)**4 + 3.0620589999999996d6*Cos(theta)**6 - 1.0935925d7*Cos(theta)**8 + 1.9684665000000001d7*Cos(theta)**10 - 1.7298645000000001d7*Cos(theta)**12 + 5.892945000000001d6*Cos(theta)**14)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(8.15617536617317d-4*Cos(theta)*(-6.435d3 + 2.85285d5*Cos(theta)**2 - 3.594591d6*Cos(theta)**4 + 1.9684665000000001d7*Cos(theta)**6 - 5.4679625d7*Cos(theta)**8 + 8.0528175d7*Cos(theta)**10 - 5.9879925d7*Cos(theta)**12 + 1.7678835000000002d7*Cos(theta)**14)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=4.9454077258518785d-5*(6.435d3 - 8.7516d5*Cos(theta)**2 + 1.939938d7*Cos(theta)**4 - 1.6295479199999998d8*Cos(theta)**6 + 6.6927861d8*Cos(theta)**8 - 1.4872858d9*Cos(theta)**10 + 1.8253053d9*Cos(theta)**12 - 1.1633814d9*Cos(theta)**14 + 3.00540195d8*Cos(theta)**16)
   case(1)
   Ylm=-8.15617536617317d-4*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(-6.435d3 + 2.85285d5*Cos(theta)**2 - 3.594591d6*Cos(theta)**4 + 1.9684665000000001d7*Cos(theta)**6 - 5.4679625d7*Cos(theta)**8 + 8.0528175d7*Cos(theta)**10 - 5.9879925d7*Cos(theta)**12 + 1.7678835000000002d7*Cos(theta)**14)*Sin(theta)
   case(2)
   Ylm=2.233660615510501d-3*2.7182818284590446d0**((0.d0,2.d0)*phi)*(-1.43d2 + 1.9019d4*Cos(theta)**2 - 3.99399d5*Cos(theta)**4 + 3.0620589999999996d6*Cos(theta)**6 - 1.0935925d7*Cos(theta)**8 + 1.9684665000000001d7*Cos(theta)**10 - 1.7298645000000001d7*Cos(theta)**12 + 5.892945000000001d6*Cos(theta)**14)*Sin(theta)**2
   case(3)
   Ylm=-3.642990217177658d-2*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(1.43d2 - 6.006d3*Cos(theta)**2 + 6.9069d4*Cos(theta)**4 - 3.289d5*Cos(theta)**6 + 7.400250000000001d5*Cos(theta)**8 - 7.8039000000000005d5*Cos(theta)**10 + 3.10155d5*Cos(theta)**12)*Sin(theta)**3
   case(4)
   Ylm=2.937072610541639d-2*2.7182818284590446d0**((0.d0,4.d0)*phi)*(1.1d1 - 1.3860000000000001d3*Cos(theta)**2 + 2.6565d4*Cos(theta)**4 - 1.7710000000000001d5*Cos(theta)**6 + 5.1232500000000005d5*Cos(theta)**8 - 6.6033d5*Cos(theta)**10 + 3.10155d5*Cos(theta)**12)*Sin(theta)**4
   case(5)
   Ylm=-2.2202182028664117d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(-2.31d2 + 8.855d3*Cos(theta)**2 - 8.855d4*Cos(theta)**4 + 3.4155d5*Cos(theta)**6 - 5.50275d5*Cos(theta)**8 + 3.10155d5*Cos(theta)**10)*Sin(theta)**5
   case(6)
   Ylm=1.5699313469606495d-2*2.7182818284590446d0**((0.d0,6.d0)*phi)*(-2.1d1 + 2.415d3*Cos(theta)**2 - 4.025d4*Cos(theta)**4 + 2.1734999999999998d5*Cos(theta)**6 - 4.50225d5*Cos(theta)**8 + 3.10155d5*Cos(theta)**10)*Sin(theta)**6
   case(7)
   Ylm=-2.380918771942937d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*(2.1d1 - 7.d2*Cos(theta)**2 + 5.67d3*Cos(theta)**4 - 1.5659999999999998d4*Cos(theta)**6 + 1.3485d4*Cos(theta)**8)*Sin(theta)**7
   case(8)
   Ylm=4.860030091895121d-2*2.7182818284590446d0**((0.d0,8.d0)*phi)*(7.d0 - 7.d2*Cos(theta)**2 + 9.45d3*Cos(theta)**4 - 3.654d4*Cos(theta)**6 + 4.0455000000000005d4*Cos(theta)**8)*Sin(theta)**8
   case(9)
   Ylm=-1.3746240938998877d-1*2.7182818284590446d0**((0.d0,9.d0)*phi)*Cos(theta)*(-3.5d1 + 9.45d2*Cos(theta)**2 - 5.481d3*Cos(theta)**4 + 8.091000000000001d3*Cos(theta)**6)*Sin(theta)**9
   case(10)
   Ylm=7.13257419188424d-2*2.7182818284590446d0**((0.d0,1.d1)*phi)*(-5.d0 + 4.050000000000001d2*Cos(theta)**2 - 3.915d3*Cos(theta)**4 + 8.091000000000001d3*Cos(theta)**6)*Sin(theta)**10
   case(11)
   Ylm=-3.0260949470385032d-1*2.7182818284590446d0**((0.d0,1.1d1)*phi)*Cos(theta)*(1.5d1 - 2.9d2*Cos(theta)**2 + 8.99d2*Cos(theta)**4)*Sin(theta)**11
   case(12)
   Ylm=1.2787585098510281d-1*2.7182818284590446d0**((0.d0,1.2d1)*phi)*(3.d0 - 1.7399999999999998d2*Cos(theta)**2 + 8.99d2*Cos(theta)**4)*Sin(theta)**12
   case(13)
   Ylm=-1.3772650648147036d0*2.7182818284590446d0**((0.d0,1.3d1)*phi)*Cos(theta)*(-3.d0 + 3.1d1*Cos(theta)**2)*Sin(theta)**13
   case(14)
   Ylm=4.355294546593892d-1*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*(-1.d0 + 3.1d1*Cos(theta)**2)*Sin(theta)**14
   case(15)
   Ylm=-3.4293623553521d0*2.7182818284590446d0**((0.d0,1.5d1)*phi)*Cos(theta)*Sin(theta)**15
   case(16)
   Ylm=6.062313441538353d-1*2.7182818284590446d0**((0.d0,1.6d1)*phi)*Sin(theta)**16
end select
case(17)
 select case(m)
   case(-17)
   Ylm=(6.150819049288287d-1*Sin(theta)**17)/2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)
   case(-16)
   Ylm=(3.586512999029811d0*Cos(theta)*Sin(theta)**16)/2.7182818284590446d0**((0.d0,1.6d1)*phi)
   case(-15)
   Ylm=(4.414692324673375d-1*(-1.d0 + 3.3000000000000003d1*Cos(theta)**2)*Sin(theta)**15)/2.7182818284590446d0**((0.d0,1.5d1)*phi)
   case(-14)
   Ylm=(4.325497426732421d0*Cos(theta)*(-1.d0 + 1.1d1*Cos(theta)**2)*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(3.8844113587215716d-1*(1.d0 - 6.2d1*Cos(theta)**2 + 3.41d2*Cos(theta)**4)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(3.1716085933129854d-1*Cos(theta)*(1.5d1 - 3.1d2*Cos(theta)**2 + 1.0230000000000001d3*Cos(theta)**4)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(7.213170808765836d-2*(-5.d0 + 4.35d2*Cos(theta)**2 - 4.495d3*Cos(theta)**4 + 9.889d3*Cos(theta)**6)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(1.4426341617531673d-1*Cos(theta)*(-3.5d1 + 1.0150000000000001d3*Cos(theta)**2 - 6.292999999999999d3*Cos(theta)**4 + 9.889d3*Cos(theta)**6)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(9.815882171674978d-3*(3.5d1 - 3.7800000000000002d3*Cos(theta)**2 + 5.481d4*Cos(theta)**4 - 2.26548d5*Cos(theta)**6 + 2.6700299999999997d5*Cos(theta)**8)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(1.501541242094947d-1*Cos(theta)*(3.5d1 - 1.26d3*Cos(theta)**2 + 1.0962d4*Cos(theta)**4 - 3.2363999999999997d4*Cos(theta)**6 + 2.9667d4*Cos(theta)**8)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(4.748290325698331d-2*(-7.d0 + 8.75d2*Cos(theta)**2 - 1.575d4*Cos(theta)**4 + 9.135d4*Cos(theta)**6 - 2.0227500000000003d5*Cos(theta)**8 + 1.48335d5*Cos(theta)**10)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(2.5716861974889755d-1*Cos(theta)*(-2.1d1 + 8.75d2*Cos(theta)**2 - 9.45d3*Cos(theta)**4 + 3.915d4*Cos(theta)**6 - 6.7425d4*Cos(theta)**8 + 4.0455000000000005d4*Cos(theta)**10)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(4.643919551304124d-2*(7.d0 - 9.66d2*Cos(theta)**2 + 2.0125d4*Cos(theta)**4 - 1.449d5*Cos(theta)**6 + 4.50225d5*Cos(theta)**8 - 6.2031d5*Cos(theta)**10 + 3.10155d5*Cos(theta)**12)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(5.492014393324373d-3*Cos(theta)*(1.001d3 - 4.6046d4*Cos(theta)**2 + 5.755749999999999d5*Cos(theta)**4 - 2.9600999999999997d6*Cos(theta)**6 + 7.153575d6*Cos(theta)**8 - 8.064029999999999d6*Cos(theta)**10 + 3.411705d6*Cos(theta)**12)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(2.2421054872776054d-3*(-1.43d2 + 2.1021d4*Cos(theta)**2 - 4.83483d5*Cos(theta)**4 + 4.029025d6*Cos(theta)**6 - 1.5540525d7*Cos(theta)**8 + 3.0045015000000004d7*Cos(theta)**10 - 2.8224105d7*Cos(theta)**12 + 1.0235115000000001d7*Cos(theta)**14)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(7.766881239787575d-3*Cos(theta)*(-7.1499999999999995d2 + 3.5035d4*Cos(theta)**2 - 4.83483d5*Cos(theta)**4 + 2.8778749999999995d6*Cos(theta)**6 - 8.633625d6*Cos(theta)**8 + 1.3656825000000001d7*Cos(theta)**10 - 1.0855425d7*Cos(theta)**12 + 3.411705d6*Cos(theta)**14)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(4.4546118987784125d-4*(7.1499999999999995d2 - 1.0868d5*Cos(theta)**2 + 2.66266d6*Cos(theta)**4 - 2.4496472d7*Cos(theta)**6 + 1.0935925d8*Cos(theta)**8 - 2.6246219999999996d8*Cos(theta)**10 + 3.4597290000000003d8*Cos(theta)**12 - 2.357178d8*Cos(theta)**14 + 6.4822395d7*Cos(theta)**16)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=5.093064253329884d-5*(1.09395d5*Cos(theta) - 5.54268d6*Cos(theta)**3 + 8.1477396d7*Cos(theta)**5 - 5.354228880000001d8*Cos(theta)**7 + 1.85910725d9*Cos(theta)**9 - 3.6506106d9*Cos(theta)**11 + 4.0718349d9*Cos(theta)**13 - 2.40432156d9*Cos(theta)**15 + 5.83401555d8*Cos(theta)**17)
   case(1)
   Ylm=-4.4546118987784125d-4*2.7182818284590446d0**((0.d0,1.d0)*phi)*(7.1499999999999995d2 - 1.0868d5*Cos(theta)**2 + 2.66266d6*Cos(theta)**4 - 2.4496472d7*Cos(theta)**6 + 1.0935925d8*Cos(theta)**8 - 2.6246219999999996d8*Cos(theta)**10 + 3.4597290000000003d8*Cos(theta)**12 - 2.357178d8*Cos(theta)**14 + 6.4822395d7*Cos(theta)**16)*Sin(theta)
   case(2)
   Ylm=7.766881239787575d-3*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(-7.1499999999999995d2 + 3.5035d4*Cos(theta)**2 - 4.83483d5*Cos(theta)**4 + 2.8778749999999995d6*Cos(theta)**6 - 8.633625d6*Cos(theta)**8 + 1.3656825000000001d7*Cos(theta)**10 - 1.0855425d7*Cos(theta)**12 + 3.411705d6*Cos(theta)**14)*Sin(theta)**2
   case(3)
   Ylm=-2.2421054872776054d-3*2.7182818284590446d0**((0.d0,3.d0)*phi)*(-1.43d2 + 2.1021d4*Cos(theta)**2 - 4.83483d5*Cos(theta)**4 + 4.029025d6*Cos(theta)**6 - 1.5540525d7*Cos(theta)**8 + 3.0045015000000004d7*Cos(theta)**10 - 2.8224105d7*Cos(theta)**12 + 1.0235115000000001d7*Cos(theta)**14)*Sin(theta)**3
   case(4)
   Ylm=5.492014393324373d-3*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(1.001d3 - 4.6046d4*Cos(theta)**2 + 5.755749999999999d5*Cos(theta)**4 - 2.9600999999999997d6*Cos(theta)**6 + 7.153575d6*Cos(theta)**8 - 8.064029999999999d6*Cos(theta)**10 + 3.411705d6*Cos(theta)**12)*Sin(theta)**4
   case(5)
   Ylm=-4.643919551304124d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*(7.d0 - 9.66d2*Cos(theta)**2 + 2.0125d4*Cos(theta)**4 - 1.449d5*Cos(theta)**6 + 4.50225d5*Cos(theta)**8 - 6.2031d5*Cos(theta)**10 + 3.10155d5*Cos(theta)**12)*Sin(theta)**5
   case(6)
   Ylm=2.5716861974889755d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*(-2.1d1 + 8.75d2*Cos(theta)**2 - 9.45d3*Cos(theta)**4 + 3.915d4*Cos(theta)**6 - 6.7425d4*Cos(theta)**8 + 4.0455000000000005d4*Cos(theta)**10)*Sin(theta)**6
   case(7)
   Ylm=-4.748290325698331d-2*2.7182818284590446d0**((0.d0,7.d0)*phi)*(-7.d0 + 8.75d2*Cos(theta)**2 - 1.575d4*Cos(theta)**4 + 9.135d4*Cos(theta)**6 - 2.0227500000000003d5*Cos(theta)**8 + 1.48335d5*Cos(theta)**10)*Sin(theta)**7
   case(8)
   Ylm=1.501541242094947d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*Cos(theta)*(3.5d1 - 1.26d3*Cos(theta)**2 + 1.0962d4*Cos(theta)**4 - 3.2363999999999997d4*Cos(theta)**6 + 2.9667d4*Cos(theta)**8)*Sin(theta)**8
   case(9)
   Ylm=-9.815882171674978d-3*2.7182818284590446d0**((0.d0,9.d0)*phi)*(3.5d1 - 3.7800000000000002d3*Cos(theta)**2 + 5.481d4*Cos(theta)**4 - 2.26548d5*Cos(theta)**6 + 2.6700299999999997d5*Cos(theta)**8)*Sin(theta)**9
   case(10)
   Ylm=1.4426341617531673d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*Cos(theta)*(-3.5d1 + 1.0150000000000001d3*Cos(theta)**2 - 6.292999999999999d3*Cos(theta)**4 + 9.889d3*Cos(theta)**6)*Sin(theta)**10
   case(11)
   Ylm=-7.213170808765836d-2*2.7182818284590446d0**((0.d0,1.1d1)*phi)*(-5.d0 + 4.35d2*Cos(theta)**2 - 4.495d3*Cos(theta)**4 + 9.889d3*Cos(theta)**6)*Sin(theta)**11
   case(12)
   Ylm=3.1716085933129854d-1*2.7182818284590446d0**((0.d0,1.2d1)*phi)*Cos(theta)*(1.5d1 - 3.1d2*Cos(theta)**2 + 1.0230000000000001d3*Cos(theta)**4)*Sin(theta)**12
   case(13)
   Ylm=-3.8844113587215716d-1*2.7182818284590446d0**((0.d0,1.3d1)*phi)*(1.d0 - 6.2d1*Cos(theta)**2 + 3.41d2*Cos(theta)**4)*Sin(theta)**13
   case(14)
   Ylm=4.325497426732421d0*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*Cos(theta)*(-1.d0 + 1.1d1*Cos(theta)**2)*Sin(theta)**14
   case(15)
   Ylm=-4.414692324673375d-1*2.7182818284590446d0**((0.d0,1.5d1)*phi)*(-1.d0 + 3.3000000000000003d1*Cos(theta)**2)*Sin(theta)**15
   case(16)
   Ylm=3.586512999029811d0*2.7182818284590446d0**((0.d0,1.6d1)*phi)*Cos(theta)*Sin(theta)**16
   case(17)
   Ylm=-6.150819049288287d-1*2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)*Sin(theta)**17
end select
case(18)
 select case(m)
   case(-18)
   Ylm=(6.2356619406092175d-1*Sin(theta)**18)/2.7182818284590446d0**((0.d0,1.7999999999999998d1)*phi)
   case(-17)
   Ylm=(3.741397164365531d0*Cos(theta)*Sin(theta)**17)/2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)
   case(-16)
   Ylm=(4.471824929732256d-1*(-1.d0 + 3.5d1*Cos(theta)**2)*Sin(theta)**16)/2.7182818284590446d0**((0.d0,1.6d1)*phi)
   case(-15)
   Ylm=(1.5054405987107193d0*Cos(theta)*(-3.d0 + 3.5d1*Cos(theta)**2)*Sin(theta)**15)/2.7182818284590446d0**((0.d0,1.5d1)*phi)
   case(-14)
   Ylm=(3.9309535590615665d-1*(1.d0 - 6.6000000000000005d1*Cos(theta)**2 + 3.85d2*Cos(theta)**4)*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(4.972306649191909d0*Cos(theta)*(1.d0 - 2.2d1*Cos(theta)**2 + 7.7d1*Cos(theta)**4)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(3.645872125527428d-1*(-1.d0 + 9.3d1*Cos(theta)**2 - 1.0230000000000001d3*Cos(theta)**4 + 2.387d3*Cos(theta)**6)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(1.0566741307889687d0*Cos(theta)*(-5.d0 + 1.55d2*Cos(theta)**2 - 1.0230000000000001d3*Cos(theta)**4 + 1.705d3*Cos(theta)**6)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(6.937405540452369d-2*(5.d0 - 5.8d2*Cos(theta)**2 + 8.99d3*Cos(theta)**4 - 3.9556000000000004d4*Cos(theta)**6 + 4.9445d4*Cos(theta)**8)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(1.7480618860989154d-2*Cos(theta)*(3.15d2 - 1.218d4*Cos(theta)**2 + 1.13274d5*Cos(theta)**4 - 3.56004d5*Cos(theta)**6 + 3.46115d5*Cos(theta)**8)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(4.787264634657012d-2*(-7.d0 + 9.45d2*Cos(theta)**2 - 1.827d4*Cos(theta)**4 + 1.13274d5*Cos(theta)**6 - 2.6700299999999997d5*Cos(theta)**8 + 2.07669d5*Cos(theta)**10)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(8.095999115069116d-1*Cos(theta)*(-7.d0 + 3.15d2*Cos(theta)**2 - 3.654d3*Cos(theta)**4 + 1.6181999999999999d4*Cos(theta)**6 - 2.9667d4*Cos(theta)**8 + 1.8879000000000001d4*Cos(theta)**10)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(4.674227268444126d-2*(7.d0 - 1.05d3*Cos(theta)**2 + 2.3625d4*Cos(theta)**4 - 1.827d5*Cos(theta)**6 + 6.068249999999999d5*Cos(theta)**8 - 8.9001d5*Cos(theta)**10 + 4.7197499999999994d5*Cos(theta)**12)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(6.351024226118053d-2*Cos(theta)*(9.1d1 - 4.55d3*Cos(theta)**2 + 6.1425d4*Cos(theta)**4 - 3.393d5*Cos(theta)**6 + 8.76525d5*Cos(theta)**8 - 1.05183d6*Cos(theta)**10 + 4.7197499999999994d5*Cos(theta)**12)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(2.4775012001276826d-2*(-1.3d1 + 2.093d3*Cos(theta)**2 - 5.2325d4*Cos(theta)**4 + 4.70925d5*Cos(theta)**6 - 1.9509750000000001d6*Cos(theta)**8 + 4.0320149999999995d6*Cos(theta)**10 - 4.0320149999999995d6*Cos(theta)**12 + 1.550775d6*Cos(theta)**14)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(1.363819524698825d-2*Cos(theta)*(-4.29d2 + 2.3023d4*Cos(theta)**2 - 3.45345d5*Cos(theta)**4 + 2.220075d6*Cos(theta)**6 - 7.153575d6*Cos(theta)**8 + 1.2096045d7*Cos(theta)**10 - 1.0235115000000001d7*Cos(theta)**12 + 3.411705d6*Cos(theta)**14)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(2.232073645068237d-3*(1.43d2 - 2.4024d4*Cos(theta)**2 + 6.44644d5*Cos(theta)**4 - 6.44644d6*Cos(theta)**6 + 3.108105d7*Cos(theta)**8 - 8.012004000000001d7*Cos(theta)**10 + 1.1289642d8*Cos(theta)**12 - 8.188092000000001d7*Cos(theta)**14 + 2.3881935000000003d7*Cos(theta)**16)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(4.842047577096091d-4*Cos(theta)*(1.2155d4 - 6.806799999999999d5*Cos(theta)**2 + 1.0958948d7*Cos(theta)**4 - 7.82782d7*Cos(theta)**6 + 2.9354324999999997d8*Cos(theta)**8 - 6.1910940000000005d8*Cos(theta)**10 + 7.381689d8*Cos(theta)**12 - 4.6399188d8*Cos(theta)**14 + 1.1940967500000002d8*Cos(theta)**16)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=2.618279463797645d-5*(-1.2155d4 + 2.078505d6*Cos(theta)**2 - 5.819814d7*Cos(theta)**4 + 6.2466003599999995d8*Cos(theta)**6 - 3.34639305d9*Cos(theta)**8 + 1.003917915d10*Cos(theta)**10 - 1.7644617900000001d10*Cos(theta)**12 + 1.8032411700000002d10*Cos(theta)**14 - 9.917826435d9*Cos(theta)**16 + 2.268783825d9*Cos(theta)**18)
   case(1)
   Ylm=-4.842047577096091d-4*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(1.2155d4 - 6.806799999999999d5*Cos(theta)**2 + 1.0958948d7*Cos(theta)**4 - 7.82782d7*Cos(theta)**6 + 2.9354324999999997d8*Cos(theta)**8 - 6.1910940000000005d8*Cos(theta)**10 + 7.381689d8*Cos(theta)**12 - 4.6399188d8*Cos(theta)**14 + 1.1940967500000002d8*Cos(theta)**16)*Sin(theta)
   case(2)
   Ylm=2.232073645068237d-3*2.7182818284590446d0**((0.d0,2.d0)*phi)*(1.43d2 - 2.4024d4*Cos(theta)**2 + 6.44644d5*Cos(theta)**4 - 6.44644d6*Cos(theta)**6 + 3.108105d7*Cos(theta)**8 - 8.012004000000001d7*Cos(theta)**10 + 1.1289642d8*Cos(theta)**12 - 8.188092000000001d7*Cos(theta)**14 + 2.3881935000000003d7*Cos(theta)**16)*Sin(theta)**2
   case(3)
   Ylm=-1.363819524698825d-2*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(-4.29d2 + 2.3023d4*Cos(theta)**2 - 3.45345d5*Cos(theta)**4 + 2.220075d6*Cos(theta)**6 - 7.153575d6*Cos(theta)**8 + 1.2096045d7*Cos(theta)**10 - 1.0235115000000001d7*Cos(theta)**12 + 3.411705d6*Cos(theta)**14)*Sin(theta)**3
   case(4)
   Ylm=2.4775012001276826d-2*2.7182818284590446d0**((0.d0,4.d0)*phi)*(-1.3d1 + 2.093d3*Cos(theta)**2 - 5.2325d4*Cos(theta)**4 + 4.70925d5*Cos(theta)**6 - 1.9509750000000001d6*Cos(theta)**8 + 4.0320149999999995d6*Cos(theta)**10 - 4.0320149999999995d6*Cos(theta)**12 + 1.550775d6*Cos(theta)**14)*Sin(theta)**4
   case(5)
   Ylm=-6.351024226118053d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(9.1d1 - 4.55d3*Cos(theta)**2 + 6.1425d4*Cos(theta)**4 - 3.393d5*Cos(theta)**6 + 8.76525d5*Cos(theta)**8 - 1.05183d6*Cos(theta)**10 + 4.7197499999999994d5*Cos(theta)**12)*Sin(theta)**5
   case(6)
   Ylm=4.674227268444126d-2*2.7182818284590446d0**((0.d0,6.d0)*phi)*(7.d0 - 1.05d3*Cos(theta)**2 + 2.3625d4*Cos(theta)**4 - 1.827d5*Cos(theta)**6 + 6.068249999999999d5*Cos(theta)**8 - 8.9001d5*Cos(theta)**10 + 4.7197499999999994d5*Cos(theta)**12)*Sin(theta)**6
   case(7)
   Ylm=-8.095999115069116d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*(-7.d0 + 3.15d2*Cos(theta)**2 - 3.654d3*Cos(theta)**4 + 1.6181999999999999d4*Cos(theta)**6 - 2.9667d4*Cos(theta)**8 + 1.8879000000000001d4*Cos(theta)**10)*Sin(theta)**7
   case(8)
   Ylm=4.787264634657012d-2*2.7182818284590446d0**((0.d0,8.d0)*phi)*(-7.d0 + 9.45d2*Cos(theta)**2 - 1.827d4*Cos(theta)**4 + 1.13274d5*Cos(theta)**6 - 2.6700299999999997d5*Cos(theta)**8 + 2.07669d5*Cos(theta)**10)*Sin(theta)**8
   case(9)
   Ylm=-1.7480618860989154d-2*2.7182818284590446d0**((0.d0,9.d0)*phi)*Cos(theta)*(3.15d2 - 1.218d4*Cos(theta)**2 + 1.13274d5*Cos(theta)**4 - 3.56004d5*Cos(theta)**6 + 3.46115d5*Cos(theta)**8)*Sin(theta)**9
   case(10)
   Ylm=6.937405540452369d-2*2.7182818284590446d0**((0.d0,1.d1)*phi)*(5.d0 - 5.8d2*Cos(theta)**2 + 8.99d3*Cos(theta)**4 - 3.9556000000000004d4*Cos(theta)**6 + 4.9445d4*Cos(theta)**8)*Sin(theta)**10
   case(11)
   Ylm=-1.0566741307889687d0*2.7182818284590446d0**((0.d0,1.1d1)*phi)*Cos(theta)*(-5.d0 + 1.55d2*Cos(theta)**2 - 1.0230000000000001d3*Cos(theta)**4 + 1.705d3*Cos(theta)**6)*Sin(theta)**11
   case(12)
   Ylm=3.645872125527428d-1*2.7182818284590446d0**((0.d0,1.2d1)*phi)*(-1.d0 + 9.3d1*Cos(theta)**2 - 1.0230000000000001d3*Cos(theta)**4 + 2.387d3*Cos(theta)**6)*Sin(theta)**12
   case(13)
   Ylm=-4.972306649191909d0*2.7182818284590446d0**((0.d0,1.3d1)*phi)*Cos(theta)*(1.d0 - 2.2d1*Cos(theta)**2 + 7.7d1*Cos(theta)**4)*Sin(theta)**13
   case(14)
   Ylm=3.9309535590615665d-1*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*(1.d0 - 6.6000000000000005d1*Cos(theta)**2 + 3.85d2*Cos(theta)**4)*Sin(theta)**14
   case(15)
   Ylm=-1.5054405987107193d0*2.7182818284590446d0**((0.d0,1.5d1)*phi)*Cos(theta)*(-3.d0 + 3.5d1*Cos(theta)**2)*Sin(theta)**15
   case(16)
   Ylm=4.471824929732256d-1*2.7182818284590446d0**((0.d0,1.6d1)*phi)*(-1.d0 + 3.5d1*Cos(theta)**2)*Sin(theta)**16
   case(17)
   Ylm=-3.741397164365531d0*2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)*Cos(theta)*Sin(theta)**17
   case(18)
   Ylm=6.2356619406092175d-1*2.7182818284590446d0**((0.d0,1.7999999999999998d1)*phi)*Sin(theta)**18
end select
case(19)
 select case(m)
   case(-19)
   Ylm=(6.317177321159495d-1*Sin(theta)**19)/2.7182818284590446d0**((0.d0,1.9d1)*phi)
   case(-18)
   Ylm=(3.8941696337793634d0*Cos(theta)*Sin(theta)**18)/2.7182818284590446d0**((0.d0,1.7999999999999998d1)*phi)
   case(-17)
   Ylm=(4.526880247947345d-1*(-1.d0 + 3.7d1*Cos(theta)**2)*Sin(theta)**17)/2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)
   case(-16)
   Ylm=(1.56815731784496d0*Cos(theta)*(-3.d0 + 3.7d1*Cos(theta)**2)*Sin(theta)**16)/2.7182818284590446d0**((0.d0,1.6d1)*phi)
   case(-15)
   Ylm=(1.3253348292603262d-1*(3.d0 - 2.1d2*Cos(theta)**2 + 1.295d3*Cos(theta)**4)*Sin(theta)**15)/2.7182818284590446d0**((0.d0,1.5d1)*phi)
   case(-14)
   Ylm=(1.728025201322552d0*Cos(theta)*(3.d0 - 7.d1*Cos(theta)**2 + 2.59d2*Cos(theta)**4)*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(3.6841621080251614d-1*(-1.d0 + 9.9d1*Cos(theta)**2 - 1.155d3*Cos(theta)**4 + 2.8489999999999998d3*Cos(theta)**6)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(5.513948946226d0*Cos(theta)*(-1.d0 + 3.3000000000000003d1*Cos(theta)**2 - 2.31d2*Cos(theta)**4 + 4.069999999999999d2*Cos(theta)**6)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(3.501361082216343d-1*(1.d0 - 1.24d2*Cos(theta)**2 + 2.0460000000000003d3*Cos(theta)**4 - 9.548d3*Cos(theta)**6 + 1.2617d4*Cos(theta)**8)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(1.278516297800394d-1*Cos(theta)*(4.5d1 - 1.8599999999999999d3*Cos(theta)**2 + 1.8414d4*Cos(theta)**4 - 6.138d4*Cos(theta)**6 + 6.3085d4*Cos(theta)**8)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(3.7538531052373676d-2*(-9.d0 + 1.3050000000000002d3*Cos(theta)**2 - 2.697d4*Cos(theta)**4 + 1.78002d5*Cos(theta)**6 - 4.45005d5*Cos(theta)**8 + 3.6589300000000002d5*Cos(theta)**10)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(9.411407803988352d-2*Cos(theta)*(-6.3d1 + 3.045d3*Cos(theta)**2 - 3.7758000000000003d4*Cos(theta)**4 + 1.78002d5*Cos(theta)**6 - 3.46115d5*Cos(theta)**8 + 2.32841d5*Cos(theta)**10)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(4.705703901994176d-2*(7.d0 - 1.134d3*Cos(theta)**2 + 2.7405d4*Cos(theta)**4 - 2.26548d5*Cos(theta)**6 + 8.01009d5*Cos(theta)**8 - 1.246014d6*Cos(theta)**10 + 6.98523d5*Cos(theta)**12)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(6.654870278712157d-2*Cos(theta)*(9.1d1 - 4.914d3*Cos(theta)**2 + 7.1253d4*Cos(theta)**4 - 4.20732d5*Cos(theta)**6 + 1.157013d6*Cos(theta)**8 - 1.472562d6*Cos(theta)**10 + 6.98523d5*Cos(theta)**12)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(2.4900244536365705d-2*(-1.3d1 + 2.275d3*Cos(theta)**2 - 6.1425d4*Cos(theta)**4 + 5.93775d5*Cos(theta)**6 - 2.629575d6*Cos(theta)**8 + 5.785065d6*Cos(theta)**10 - 6.135675000000001d6*Cos(theta)**12 + 2.4947250000000003d6*Cos(theta)**14)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(1.57482974060158d-1*Cos(theta)*(-3.9000000000000004d1 + 2.275d3*Cos(theta)**2 - 3.6854999999999998d4*Cos(theta)**4 + 2.54475d5*Cos(theta)**6 - 8.76525d5*Cos(theta)**8 + 1.5777450000000002d6*Cos(theta)**10 - 1.415925d6*Cos(theta)**12 + 4.989450000000001d5*Cos(theta)**14)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(8.209367515029838d-3*(3.9000000000000004d1 - 7.176d3*Cos(theta)**2 + 2.093d5*Cos(theta)**4 - 2.26044d6*Cos(theta)**6 + 1.170585d7*Cos(theta)**8 - 3.225612d7*Cos(theta)**10 + 4.838418d7*Cos(theta)**12 - 3.7218600000000004d7*Cos(theta)**14 + 1.1475735d7*Cos(theta)**16)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(8.489925769333858d-4*Cos(theta)*(7.292999999999999d3 - 4.47304d5*Cos(theta)**2 + 7.82782d6*Cos(theta)**4 - 6.038603999999999d7*Cos(theta)**6 + 2.4322155d8*Cos(theta)**8 - 5.4835404d8*Cos(theta)**10 + 6.9598782d8*Cos(theta)**12 - 4.6399188d8*Cos(theta)**14 + 1.26233085d8*Cos(theta)**16)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(1.3100239871377048d-4*(-2.431d3 + 4.59459d5*Cos(theta)**2 - 1.4090076000000002d7*Cos(theta)**4 + 1.6438422d8*Cos(theta)**6 - 9.5108013d8*Cos(theta)**8 + 3.06459153d9*Cos(theta)**10 - 5.7577174200000005d9*Cos(theta)**12 + 6.26389038d9*Cos(theta)**14 - 3.653936055d9*Cos(theta)**16 + 8.83631595d8*Cos(theta)**18)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=2.688112503031131d-5*(-2.30945d5*Cos(theta) + 1.4549535d7*Cos(theta)**3 - 2.6771144400000004d8*Cos(theta)**5 + 2.2309286999999998d9*Cos(theta)**7 - 1.003917915d10*Cos(theta)**9 + 2.646692685d10*Cos(theta)**11 - 4.20756273d10*Cos(theta)**13 + 3.967130574d10*Cos(theta)**15 - 2.0419054425d10*Cos(theta)**17 + 4.418157975d9*Cos(theta)**19)
   case(1)
   Ylm=-1.3100239871377048d-4*2.7182818284590446d0**((0.d0,1.d0)*phi)*(-2.431d3 + 4.59459d5*Cos(theta)**2 - 1.4090076000000002d7*Cos(theta)**4 + 1.6438422d8*Cos(theta)**6 - 9.5108013d8*Cos(theta)**8 + 3.06459153d9*Cos(theta)**10 - 5.7577174200000005d9*Cos(theta)**12 + 6.26389038d9*Cos(theta)**14 - 3.653936055d9*Cos(theta)**16 + 8.83631595d8*Cos(theta)**18)*Sin(theta)
   case(2)
   Ylm=8.489925769333858d-4*2.7182818284590446d0**((0.d0,2.d0)*phi)*Cos(theta)*(7.292999999999999d3 - 4.47304d5*Cos(theta)**2 + 7.82782d6*Cos(theta)**4 - 6.038603999999999d7*Cos(theta)**6 + 2.4322155d8*Cos(theta)**8 - 5.4835404d8*Cos(theta)**10 + 6.9598782d8*Cos(theta)**12 - 4.6399188d8*Cos(theta)**14 + 1.26233085d8*Cos(theta)**16)*Sin(theta)**2
   case(3)
   Ylm=-8.209367515029838d-3*2.7182818284590446d0**((0.d0,3.d0)*phi)*(3.9000000000000004d1 - 7.176d3*Cos(theta)**2 + 2.093d5*Cos(theta)**4 - 2.26044d6*Cos(theta)**6 + 1.170585d7*Cos(theta)**8 - 3.225612d7*Cos(theta)**10 + 4.838418d7*Cos(theta)**12 - 3.7218600000000004d7*Cos(theta)**14 + 1.1475735d7*Cos(theta)**16)*Sin(theta)**3
   case(4)
   Ylm=1.57482974060158d-1*2.7182818284590446d0**((0.d0,4.d0)*phi)*Cos(theta)*(-3.9000000000000004d1 + 2.275d3*Cos(theta)**2 - 3.6854999999999998d4*Cos(theta)**4 + 2.54475d5*Cos(theta)**6 - 8.76525d5*Cos(theta)**8 + 1.5777450000000002d6*Cos(theta)**10 - 1.415925d6*Cos(theta)**12 + 4.989450000000001d5*Cos(theta)**14)*Sin(theta)**4
   case(5)
   Ylm=-2.4900244536365705d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*(-1.3d1 + 2.275d3*Cos(theta)**2 - 6.1425d4*Cos(theta)**4 + 5.93775d5*Cos(theta)**6 - 2.629575d6*Cos(theta)**8 + 5.785065d6*Cos(theta)**10 - 6.135675000000001d6*Cos(theta)**12 + 2.4947250000000003d6*Cos(theta)**14)*Sin(theta)**5
   case(6)
   Ylm=6.654870278712157d-2*2.7182818284590446d0**((0.d0,6.d0)*phi)*Cos(theta)*(9.1d1 - 4.914d3*Cos(theta)**2 + 7.1253d4*Cos(theta)**4 - 4.20732d5*Cos(theta)**6 + 1.157013d6*Cos(theta)**8 - 1.472562d6*Cos(theta)**10 + 6.98523d5*Cos(theta)**12)*Sin(theta)**6
   case(7)
   Ylm=-4.705703901994176d-2*2.7182818284590446d0**((0.d0,7.d0)*phi)*(7.d0 - 1.134d3*Cos(theta)**2 + 2.7405d4*Cos(theta)**4 - 2.26548d5*Cos(theta)**6 + 8.01009d5*Cos(theta)**8 - 1.246014d6*Cos(theta)**10 + 6.98523d5*Cos(theta)**12)*Sin(theta)**7
   case(8)
   Ylm=9.411407803988352d-2*2.7182818284590446d0**((0.d0,8.d0)*phi)*Cos(theta)*(-6.3d1 + 3.045d3*Cos(theta)**2 - 3.7758000000000003d4*Cos(theta)**4 + 1.78002d5*Cos(theta)**6 - 3.46115d5*Cos(theta)**8 + 2.32841d5*Cos(theta)**10)*Sin(theta)**8
   case(9)
   Ylm=-3.7538531052373676d-2*2.7182818284590446d0**((0.d0,9.d0)*phi)*(-9.d0 + 1.3050000000000002d3*Cos(theta)**2 - 2.697d4*Cos(theta)**4 + 1.78002d5*Cos(theta)**6 - 4.45005d5*Cos(theta)**8 + 3.6589300000000002d5*Cos(theta)**10)*Sin(theta)**9
   case(10)
   Ylm=1.278516297800394d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*Cos(theta)*(4.5d1 - 1.8599999999999999d3*Cos(theta)**2 + 1.8414d4*Cos(theta)**4 - 6.138d4*Cos(theta)**6 + 6.3085d4*Cos(theta)**8)*Sin(theta)**10
   case(11)
   Ylm=-3.501361082216343d-1*2.7182818284590446d0**((0.d0,1.1d1)*phi)*(1.d0 - 1.24d2*Cos(theta)**2 + 2.0460000000000003d3*Cos(theta)**4 - 9.548d3*Cos(theta)**6 + 1.2617d4*Cos(theta)**8)*Sin(theta)**11
   case(12)
   Ylm=5.513948946226d0*2.7182818284590446d0**((0.d0,1.2d1)*phi)*Cos(theta)*(-1.d0 + 3.3000000000000003d1*Cos(theta)**2 - 2.31d2*Cos(theta)**4 + 4.069999999999999d2*Cos(theta)**6)*Sin(theta)**12
   case(13)
   Ylm=-3.6841621080251614d-1*2.7182818284590446d0**((0.d0,1.3d1)*phi)*(-1.d0 + 9.9d1*Cos(theta)**2 - 1.155d3*Cos(theta)**4 + 2.8489999999999998d3*Cos(theta)**6)*Sin(theta)**13
   case(14)
   Ylm=1.728025201322552d0*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*Cos(theta)*(3.d0 - 7.d1*Cos(theta)**2 + 2.59d2*Cos(theta)**4)*Sin(theta)**14
   case(15)
   Ylm=-1.3253348292603262d-1*2.7182818284590446d0**((0.d0,1.5d1)*phi)*(3.d0 - 2.1d2*Cos(theta)**2 + 1.295d3*Cos(theta)**4)*Sin(theta)**15
   case(16)
   Ylm=1.56815731784496d0*2.7182818284590446d0**((0.d0,1.6d1)*phi)*Cos(theta)*(-3.d0 + 3.7d1*Cos(theta)**2)*Sin(theta)**16
   case(17)
   Ylm=-4.526880247947345d-1*2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)*(-1.d0 + 3.7d1*Cos(theta)**2)*Sin(theta)**17
   case(18)
   Ylm=3.8941696337793634d0*2.7182818284590446d0**((0.d0,1.7999999999999998d1)*phi)*Cos(theta)*Sin(theta)**18
   case(19)
   Ylm=-6.317177321159495d-1*2.7182818284590446d0**((0.d0,1.9d1)*phi)*Sin(theta)**19
end select
case(20)
 select case(m)
   case(-20)
   Ylm=(6.395654582577624d-1*Sin(theta)**20)/2.7182818284590446d0**((0.d0,2.d1)*phi)
   case(-19)
   Ylm=(4.044967121727748d0*Cos(theta)*Sin(theta)**19)/2.7182818284590446d0**((0.d0,1.9d1)*phi)
   case(-18)
   Ylm=(4.580023375802296d-1*(-1.d0 + 3.9000000000000004d1*Cos(theta)**2)*Sin(theta)**18)/2.7182818284590446d0**((0.d0,1.7999999999999998d1)*phi)
   case(-17)
   Ylm=(4.890126797957373d0*Cos(theta)*(-1.d0 + 1.3d1*Cos(theta)**2)*Sin(theta)**17)/2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)
   case(-16)
   Ylm=(4.0196594668949075d-1*(1.d0 - 7.4d1*Cos(theta)**2 + 4.81d2*Cos(theta)**4)*Sin(theta)**16)/2.7182818284590446d0**((0.d0,1.6d1)*phi)
   case(-15)
   Ylm=(3.5952927257510314d-1*Cos(theta)*(1.5d1 - 3.7d2*Cos(theta)**2 + 1.443d3*Cos(theta)**4)*Sin(theta)**15)/2.7182818284590446d0**((0.d0,1.5d1)*phi)
   case(-14)
   Ylm=(3.7214815286923244d-1*(-1.d0 + 1.05d2*Cos(theta)**2 - 1.295d3*Cos(theta)**4 + 3.367d3*Cos(theta)**6)*Sin(theta)**14)/2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)
   case(-13)
   Ylm=(5.7412220779889385d0*Cos(theta)*(-1.d0 + 3.5d1*Cos(theta)**2 - 2.59d2*Cos(theta)**4 + 4.81d2*Cos(theta)**6)*Sin(theta)**13)/2.7182818284590446d0**((0.d0,1.3d1)*phi)
   case(-12)
   Ylm=(3.5334779281156066d-1*(1.d0 - 1.32d2*Cos(theta)**2 + 2.31d3*Cos(theta)**4 - 1.1396000000000002d4*Cos(theta)**6 + 1.5873000000000002d4*Cos(theta)**8)*Sin(theta)**12)/2.7182818284590446d0**((0.d0,1.2d1)*phi)
   case(-11)
   Ylm=(1.99883696331483d0*Cos(theta)*(3.d0 - 1.32d2*Cos(theta)**2 + 1.3860000000000001d3*Cos(theta)**4 - 4.884d3*Cos(theta)**6 + 5.291d3*Cos(theta)**8)*Sin(theta)**11)/2.7182818284590446d0**((0.d0,1.1d1)*phi)
   case(-10)
   Ylm=(1.1352631080451239d-1*(-3.d0 + 4.65d2*Cos(theta)**2 - 1.0230000000000001d4*Cos(theta)**4 + 7.161d4*Cos(theta)**6 - 1.89255d5*Cos(theta)**8 + 1.64021d5*Cos(theta)**10)*Sin(theta)**10)/2.7182818284590446d0**((0.d0,1.d1)*phi)
   case(-9)
   Ylm=(6.8743595021332755d-1*Cos(theta)*(-9.d0 + 4.65d2*Cos(theta)**2 - 6.138d3*Cos(theta)**4 + 3.069d4*Cos(theta)**6 - 6.3085d4*Cos(theta)**8 + 4.4733d4*Cos(theta)**10)*Sin(theta)**9)/2.7182818284590446d0**((0.d0,9.d0)*phi)
   case(-8)
   Ylm=(1.1055130486827582d-1*(3.d0 - 5.220000000000001d2*Cos(theta)**2 + 1.3485d4*Cos(theta)**4 - 1.18668d5*Cos(theta)**6 + 4.45005d5*Cos(theta)**8 - 7.3178600000000005d5*Cos(theta)**10 + 4.32419d5*Cos(theta)**12)*Sin(theta)**8)/2.7182818284590446d0**((0.d0,8.d0)*phi)
   case(-7)
   Ylm=(3.0131206709041236d-1*Cos(theta)*(2.1d1 - 1.218d3*Cos(theta)**2 + 1.8879000000000001d4*Cos(theta)**4 - 1.18668d5*Cos(theta)**6 + 3.46115d5*Cos(theta)**8 - 4.65682d5*Cos(theta)**10 + 2.32841d5*Cos(theta)**12)*Sin(theta)**7)/2.7182818284590446d0**((0.d0,7.d0)*phi)
   case(-6)
   Ylm=(3.2545422935237256d-1*(-1.d0 + 1.8900000000000001d2*Cos(theta)**2 - 5.481d3*Cos(theta)**4 + 5.6637d4*Cos(theta)**6 - 2.6700299999999997d5*Cos(theta)**8 + 6.2300699999999996d5*Cos(theta)**10 - 6.98523d5*Cos(theta)**12 + 2.99367d5*Cos(theta)**14)*Sin(theta)**6)/2.7182818284590446d0**((0.d0,6.d0)*phi)
   case(-5)
   Ylm=(9.888009307470726d-2*Cos(theta)*(-6.5d1 + 4.095d3*Cos(theta)**2 - 7.1253d4*Cos(theta)**4 + 5.25915d5*Cos(theta)**6 - 1.9283549999999998d6*Cos(theta)**8 + 3.681405d6*Cos(theta)**10 - 3.492615d6*Cos(theta)**12 + 1.297257d6*Cos(theta)**14)*Sin(theta)**5)/2.7182818284590446d0**((0.d0,5.d0)*phi)
   case(-4)
   Ylm=(2.4720023268676816d-2*(1.3d1 - 2.6d3*Cos(theta)**2 + 8.19d4*Cos(theta)**4 - 9.500399999999999d5*Cos(theta)**6 + 5.25915d6*Cos(theta)**8 - 1.542684d7*Cos(theta)**10 + 2.45427d7*Cos(theta)**12 - 1.99578d7*Cos(theta)**14 + 6.4862850000000005d6*Cos(theta)**16)*Sin(theta)**4)/2.7182818284590446d0**((0.d0,4.d0)*phi)
   case(-3)
   Ylm=(9.790588120722628d-3*Cos(theta)*(6.630000000000001d2 - 4.42d4*Cos(theta)**2 + 8.3538d5*Cos(theta)**4 - 6.9217200000000005d6*Cos(theta)**6 + 2.980185d7*Cos(theta)**8 - 7.152444d7*Cos(theta)**10 + 9.62829d7*Cos(theta)**12 - 6.785652d7*Cos(theta)**14 + 1.9458855d7*Cos(theta)**16)*Sin(theta)**3)/2.7182818284590446d0**((0.d0,3.d0)*phi)
   case(-2)
   Ylm=(1.4435434644262424d-3*(-2.21d2 + 4.5747d4*Cos(theta)**2 - 1.5249d6*Cos(theta)**4 + 1.9213740000000001d7*Cos(theta)**6 - 1.1939967d8*Cos(theta)**8 + 4.1126553d8*Cos(theta)**10 - 8.2253106d8*Cos(theta)**12 + 9.490743d8*Cos(theta)**14 - 5.85262485d8*Cos(theta)**16 + 1.4918455499999999d8*Cos(theta)**18)*Sin(theta)**2)/2.7182818284590446d0**((0.d0,2.d0)*phi)
   case(-1)
   Ylm=(1.4121203757760976d-4*Cos(theta)*(-4.6189d4 + 3.187041d6*Cos(theta)**2 - 6.374082d7*Cos(theta)**4 + 5.736673799999999d8*Cos(theta)**6 - 2.7727256700000003d9*Cos(theta)**8 + 7.81404507d9*Cos(theta)**10 - 1.3223768580000002d10*Cos(theta)**12 + 1.3223768580000002d10*Cos(theta)**14 - 7.195285845d9*Cos(theta)**16 + 1.641030105d9*Cos(theta)**18)*Sin(theta))/2.7182818284590446d0**((0.d0,1.d0)*phi)
   case(0)
   Ylm=6.8904418886600185d-6*(4.6189d4 - 9.69969d6*Cos(theta)**2 + 3.34639305d8*Cos(theta)**4 - 4.4618573999999995d9*Cos(theta)**6 + 3.011753745d10*Cos(theta)**8 - 1.1645447814d11*Cos(theta)**10 + 2.7349157745d11*Cos(theta)**12 - 3.967130574d11*Cos(theta)**14 + 3.47123925225d11*Cos(theta)**16 - 1.6789000305d11*Cos(theta)**18 + 3.4461632205d10*Cos(theta)**20)
   case(1)
   Ylm=-1.4121203757760976d-4*2.7182818284590446d0**((0.d0,1.d0)*phi)*Cos(theta)*(-4.6189d4 + 3.187041d6*Cos(theta)**2 - 6.374082d7*Cos(theta)**4 + 5.736673799999999d8*Cos(theta)**6 - 2.7727256700000003d9*Cos(theta)**8 + 7.81404507d9*Cos(theta)**10 - 1.3223768580000002d10*Cos(theta)**12 + 1.3223768580000002d10*Cos(theta)**14 - 7.195285845d9*Cos(theta)**16 + 1.641030105d9*Cos(theta)**18)*Sin(theta)
   case(2)
   Ylm=1.4435434644262424d-3*2.7182818284590446d0**((0.d0,2.d0)*phi)*(-2.21d2 + 4.5747d4*Cos(theta)**2 - 1.5249d6*Cos(theta)**4 + 1.9213740000000001d7*Cos(theta)**6 - 1.1939967d8*Cos(theta)**8 + 4.1126553d8*Cos(theta)**10 - 8.2253106d8*Cos(theta)**12 + 9.490743d8*Cos(theta)**14 - 5.85262485d8*Cos(theta)**16 + 1.4918455499999999d8*Cos(theta)**18)*Sin(theta)**2
   case(3)
   Ylm=-9.790588120722628d-3*2.7182818284590446d0**((0.d0,3.d0)*phi)*Cos(theta)*(6.630000000000001d2 - 4.42d4*Cos(theta)**2 + 8.3538d5*Cos(theta)**4 - 6.9217200000000005d6*Cos(theta)**6 + 2.980185d7*Cos(theta)**8 - 7.152444d7*Cos(theta)**10 + 9.62829d7*Cos(theta)**12 - 6.785652d7*Cos(theta)**14 + 1.9458855d7*Cos(theta)**16)*Sin(theta)**3
   case(4)
   Ylm=2.4720023268676816d-2*2.7182818284590446d0**((0.d0,4.d0)*phi)*(1.3d1 - 2.6d3*Cos(theta)**2 + 8.19d4*Cos(theta)**4 - 9.500399999999999d5*Cos(theta)**6 + 5.25915d6*Cos(theta)**8 - 1.542684d7*Cos(theta)**10 + 2.45427d7*Cos(theta)**12 - 1.99578d7*Cos(theta)**14 + 6.4862850000000005d6*Cos(theta)**16)*Sin(theta)**4
   case(5)
   Ylm=-9.888009307470726d-2*2.7182818284590446d0**((0.d0,5.d0)*phi)*Cos(theta)*(-6.5d1 + 4.095d3*Cos(theta)**2 - 7.1253d4*Cos(theta)**4 + 5.25915d5*Cos(theta)**6 - 1.9283549999999998d6*Cos(theta)**8 + 3.681405d6*Cos(theta)**10 - 3.492615d6*Cos(theta)**12 + 1.297257d6*Cos(theta)**14)*Sin(theta)**5
   case(6)
   Ylm=3.2545422935237256d-1*2.7182818284590446d0**((0.d0,6.d0)*phi)*(-1.d0 + 1.8900000000000001d2*Cos(theta)**2 - 5.481d3*Cos(theta)**4 + 5.6637d4*Cos(theta)**6 - 2.6700299999999997d5*Cos(theta)**8 + 6.2300699999999996d5*Cos(theta)**10 - 6.98523d5*Cos(theta)**12 + 2.99367d5*Cos(theta)**14)*Sin(theta)**6
   case(7)
   Ylm=-3.0131206709041236d-1*2.7182818284590446d0**((0.d0,7.d0)*phi)*Cos(theta)*(2.1d1 - 1.218d3*Cos(theta)**2 + 1.8879000000000001d4*Cos(theta)**4 - 1.18668d5*Cos(theta)**6 + 3.46115d5*Cos(theta)**8 - 4.65682d5*Cos(theta)**10 + 2.32841d5*Cos(theta)**12)*Sin(theta)**7
   case(8)
   Ylm=1.1055130486827582d-1*2.7182818284590446d0**((0.d0,8.d0)*phi)*(3.d0 - 5.220000000000001d2*Cos(theta)**2 + 1.3485d4*Cos(theta)**4 - 1.18668d5*Cos(theta)**6 + 4.45005d5*Cos(theta)**8 - 7.3178600000000005d5*Cos(theta)**10 + 4.32419d5*Cos(theta)**12)*Sin(theta)**8
   case(9)
   Ylm=-6.8743595021332755d-1*2.7182818284590446d0**((0.d0,9.d0)*phi)*Cos(theta)*(-9.d0 + 4.65d2*Cos(theta)**2 - 6.138d3*Cos(theta)**4 + 3.069d4*Cos(theta)**6 - 6.3085d4*Cos(theta)**8 + 4.4733d4*Cos(theta)**10)*Sin(theta)**9
   case(10)
   Ylm=1.1352631080451239d-1*2.7182818284590446d0**((0.d0,1.d1)*phi)*(-3.d0 + 4.65d2*Cos(theta)**2 - 1.0230000000000001d4*Cos(theta)**4 + 7.161d4*Cos(theta)**6 - 1.89255d5*Cos(theta)**8 + 1.64021d5*Cos(theta)**10)*Sin(theta)**10
   case(11)
   Ylm=-1.99883696331483d0*2.7182818284590446d0**((0.d0,1.1d1)*phi)*Cos(theta)*(3.d0 - 1.32d2*Cos(theta)**2 + 1.3860000000000001d3*Cos(theta)**4 - 4.884d3*Cos(theta)**6 + 5.291d3*Cos(theta)**8)*Sin(theta)**11
   case(12)
   Ylm=3.5334779281156066d-1*2.7182818284590446d0**((0.d0,1.2d1)*phi)*(1.d0 - 1.32d2*Cos(theta)**2 + 2.31d3*Cos(theta)**4 - 1.1396000000000002d4*Cos(theta)**6 + 1.5873000000000002d4*Cos(theta)**8)*Sin(theta)**12
   case(13)
   Ylm=-5.7412220779889385d0*2.7182818284590446d0**((0.d0,1.3d1)*phi)*Cos(theta)*(-1.d0 + 3.5d1*Cos(theta)**2 - 2.59d2*Cos(theta)**4 + 4.81d2*Cos(theta)**6)*Sin(theta)**13
   case(14)
   Ylm=3.7214815286923244d-1*2.7182818284590446d0**((0.d0,1.4000000000000001d1)*phi)*(-1.d0 + 1.05d2*Cos(theta)**2 - 1.295d3*Cos(theta)**4 + 3.367d3*Cos(theta)**6)*Sin(theta)**14
   case(15)
   Ylm=-3.5952927257510314d-1*2.7182818284590446d0**((0.d0,1.5d1)*phi)*Cos(theta)*(1.5d1 - 3.7d2*Cos(theta)**2 + 1.443d3*Cos(theta)**4)*Sin(theta)**15
   case(16)
   Ylm=4.0196594668949075d-1*2.7182818284590446d0**((0.d0,1.6d1)*phi)*(1.d0 - 7.4d1*Cos(theta)**2 + 4.81d2*Cos(theta)**4)*Sin(theta)**16
   case(17)
   Ylm=-4.890126797957373d0*2.7182818284590446d0**((0.d0,1.7000000000000002d1)*phi)*Cos(theta)*(-1.d0 + 1.3d1*Cos(theta)**2)*Sin(theta)**17
   case(18)
   Ylm=4.580023375802296d-1*2.7182818284590446d0**((0.d0,1.7999999999999998d1)*phi)*(-1.d0 + 3.9000000000000004d1*Cos(theta)**2)*Sin(theta)**18
   case(19)
   Ylm=-4.044967121727748d0*2.7182818284590446d0**((0.d0,1.9d1)*phi)*Cos(theta)*Sin(theta)**19
   case(20)
   Ylm=6.395654582577624d-1*2.7182818284590446d0**((0.d0,2.d1)*phi)*Sin(theta)**20
end select
end select
ylm_r=real(Ylm)
ylm_i=aimag(Ylm)
end subroutine
