module constant_module
  !
  ! DESCRIPTION
  !  this module contains mathematical and physical constant and 
  !  geophysical and meteorological parameter
  !
  use type_module
  implicit none
  public

  !
  ! MATHEMATICAL CONSTANT
  !  math_pi : circular constant
  !
  real(kind=dp), parameter :: math_pi  = 3.14159265358979323846_dp

  ! 
  ! GEOPHYSICAL PARAMETER
  !  earth_radius  : radius of the earth [m]
  !  earth_daysec  : seconds per day [s]
  !  earth_omega   : angular velocity of the earth [s^-1]
  !  earth_gravity : acceleration of gravity [m/s^2]
  !
  real(kind=dp), parameter :: earth_radius = 6.371e6_dp             
  real(kind=dp), parameter :: earth_daysec = 86400.0_dp             
  real(kind=dp), parameter :: earth_omega  = 2.0_dp*math_pi/earth_daysec 
  real(kind=dp), parameter :: earth_gravity = 9.80665_dp            

  ! 
  ! AIR CONSTANT
  !  air_rd    : gas constant of dry air [J K^-1 Kg^-1]
  !  air_cp    : specific heat at constant pressure [J K^-1 Kg^-1]
  !  air_cv    : specific heat at constant volume [J K^-1 Kg^-1]
  !  air_kappa : Rd/Cp
  !  air_gamma : Cp/Cv
  !
  real(kind=dp), parameter :: air_rd = 287.0_dp    
  real(kind=dp), parameter :: air_cp = 1004.0_dp   
  real(kind=dp), parameter :: air_cv =  717.0_dp   
  real(kind=dp), parameter :: air_kappa = air_rd/air_cp
  real(kind=dp), parameter :: air_gamma = air_cp/air_cv

  !
  ! INTERNATIONAL STANDART ATMOSPHERE
  !  isa_gamma : lappase late of troposphere [K/m]
  !
  real(kind=dp), parameter :: isa_gamma = 0.0065   

end module constant_module
