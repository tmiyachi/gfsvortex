module debug_module
  use type_module
  implicit none

contains
 
  subroutine write_grads(unit,fname,xn,yn,zn,svn,vnn,x,y,data,yrev)
    integer(kind=i4b), intent(in) :: unit,xn,yn,zn,svn,vnn
    character(len=*), intent(in)  :: fname
    real(kind=sp), intent(in)     :: data(xn*yn*(svn+zn*vnn)),x(xn),y(yn)
    logical, intent(in),optional  :: yrev
    integer n

    open(unit,file=fname,recl=4*xn*yn*(svn+vnn*zn),form='unformatted',access='direct')
    write(unit,rec=1) data
    close(unit)

    open(unit,file=fname//'.ctl')
    write(unit,'("dset ",a)') fname
    if (present(yrev)) then
       if (yrev) write(unit,'("options yrev")')
    end if
    write(unit,'("undef -9.99E+33")')
    write(unit,'("title debug_module")')
    write(unit,'("xdef",i6," levels")') xn
    write(unit,'(5f12.6)') x
    write(unit,'("ydef",i6," levels")') yn
    write(unit,'(5f12.6)') y
    write(unit,'("zdef",i6," linear 1 1")') zn
    write(unit, '(a)') 'tdef 1 linear 00z01jan1900 1yr'
    write(unit,'("vars",i6)') vnn+svn
    if (svn > 0) then
       do n = 1, svn
          write(unit,'("varsfc",i0," ",i3," 99 **")') n,0
       end do
    end if
    if (vnn > 0) then
       do n = 1, vnn
          write(unit,'("varlev",i0," ",i3," 99 **")') n,zn
       end do
    end if
    
    write(unit,'("endvars")')
    close(unit)

  end subroutine write_grads

end module debug_module
