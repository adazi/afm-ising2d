! ising.f90
! Author: Adam Iaizzi
! Updated 2020-03-24
! Based on a code by Anders Sandvik (see README for link)


!----------------------------------------------------------------------!
! Metropolis-algorithm simulation of the two-dimensional Ising model   !
! antiferromagnet with a uniform external field
! Reads the following from a file 'read.in':                           !
!    ll,tt, hh                                                         !
!    steps1,steps2,bins                                                !
! where: ll     = linear system size (integer)                         !
!        tt     = temperature (T/J)                                    !
!        hh     = field (h/J)                                          !
!        steps1 = # of sweeps for equilibration                        !
!        steps2 = # of sweeps per bin data collection                  !
!        bins   = # of data bins                                       !
! Random numbers seed (and integer) is read from 'seed.in'.            !
! Output files are save separated .txt files                           !
! amag.txt --> T  <m> <m^2> <m^4> (uniform magnetization, real, extensive)!
! enrg.txt --> T  <en> <en^2> (real, extensive)                       !
! smag.txt --> T  <m> <m^2> <m^4> (stagger magnetization, real, extensive)!
! spins.txt --> sigma_0 sigma_1 sigma_2 .... (spin config, integer)    !
! census.txt --> T  # of spins in each local config (see README)       !
!                                                                      !
!   See README.md for more info                                        !
!----------------------------------------------------------------------!

!----------------------!
 module systemvariables
 implicit none

 integer :: ll                       ! system length
 integer :: nn                       ! number of spins (nn=ll*ll)
 real(8) :: pflip(-4:4,-1:1)        ! flip probabilities
 integer, allocatable :: spin(:)     ! spin array

 end module systemvariables
!--------------------------!

!----------------------!
 module runSums
 implicit none
 ! array for census of local configurations
 integer :: census(-4:5) = 0
 ! x= center spin
 ! y= sum of neighboring spins
 !if x=-1, type=y. If x=+1, type=y+1
 end module runSums
!--------------------------!


!---------------!
 program ising2d
!---------------!
 use systemvariables
 use runSums
 implicit none

 integer :: i,j,m,bins,steps1,steps2,x,y
 real(8) :: crand,tt,mav,mav2,mav4,smav,smav2,smav4,en,en2,runsum
 real(8) :: hh
 integer, allocatable :: alt(:)
 ! function for calculating total flip probability
 real(8) :: sumFlips

 open(1,file='read.in',status='old')
 read(1,*)ll,tt,hh
 read(1,*)steps1,steps2,bins
 close(1)
 nn=ll*ll
 
 call initran(1)

 print*, "h=", hh, " t = ", tt
 ! set up flip probabilities
 do i=-4,4
    pflip(i,-1)=exp((2.d0*dble(i)+2*hh)/tt)
    pflip(i,0)=0.0d0
    pflip(i,+1)=exp((2.d0*dble(i)-2*hh)/tt)
 end do

 ! allocate arrays
 allocate (spin(0:nn-1))
 allocate (alt(0:nn-1))
 ! initialize random spin configuration and build checkerboard alternating sign array
 do i=0,nn-1
    spin(i)=2*int(2.d0*crand())-1
    x = mod(i,ll)
    y = i/ll
    alt(i) = (2*mod(x,2)-1)*(2*mod(y,2)-1)
 enddo

 ! equilibrate
 do i=1,steps1
    call mcstep()
 enddo

 ! do simulation, loop over bins
 do j=1,bins
    print*, "Bin ", j, " of ", bins
    ! reset magnetization sums
    mav=0.d0
    mav2=0.d0
    mav4=0.d0
    smav=0.d0
    smav2=0.d0
    smav4=0.d0
    census=0
    ! reset energy sums
    en=0.d0
    en2=0.d0
    ! do monte carlo sweeps
    do i=1,steps2
        !perfom mc sweep 
       call mcstep()
       !sum magnetization
       m=sum(spin)
       mav= mav+(1.d0*m)/steps2;
       mav2=mav2+(1.d0*m)**2/steps2;
       mav4=mav4+(1.d0*m)**4/steps2;

       !calculate energy, use runsum (real) for running sum
       runsum=0-1.0d0*m*hh
       !interactions
       do y=0,ll-1
            do x=0,ll-1
                runsum=runsum+spin(ll*y+x)*(spin(ll*mod(y+1,ll)+x)+spin(ll*y+mod(x+1,ll)))
            enddo
       enddo
       ! add energy to running sums (with field contribution)
       en=en+(1.d0*runsum)/steps2;
       en2=en2+(1.d0*runsum)**2/steps2;

       !sum staggered magnetization
       m=dot_product(spin,alt)
       smav= smav+(1.d0*m)/steps2;
       smav2=smav2+(1.d0*m)**2/steps2;
       smav4=smav4+(1.d0*m)**4/steps2;
    enddo
    
    !do census (only once per bin)
    call doCensus()

    !write data to file at the end of each bin
    
    ! uniform magnetization: T, <m>, <m^2>, <m^4>
    open(10,file='amag.txt',position='append')
    write(10,'(f8.3,3en16.6)') tt,mav,mav2,mav4
    close(10)
    ! staggered magnetization: T, <sm>, <sm^2>, <sm^4>
    open(10,file='smag.txt',position='append')
    write(10,'(f8.3,3en16.6)') tt,smav,smav2,smav4
    close(10)
    ! energy: T, <E>, <E^2>
    open(10,file='enrg.txt',position='append')
    write(10,'(f8.3,2en16.6)') tt,en,en2
    close(10)
    ! census: T, census(-4).... census(5), sum of all flip probabilities
    open(10,file='census.txt',position='append')
    write(10,'(f8.3,10i10,e10.3)') tt, census(:), sumFlips()
    close(10)
    
    !write spin configuration out to file at end of bin
    open(10,file='spins.txt',position='append')
    write(10,*) spin(:)
    close(10)
    
 enddo

 deallocate(spin)
 deallocate(alt)

!-------------------!
 end program ising2d
!-------------------!

!========================================================
!========================================================
!            EXTERNAL SUBROUTINES AND FUNCTIONS
!========================================================
!========================================================

!----------------------------------------------!
! Does census of configurations !
!----------------------------------------------!
subroutine doCensus()
!-------------------!
 use systemvariables
 use runSums
 implicit none
 
 integer :: x,y,xp1,yp1,xm1,ym1,ii,neighbors
 
 !loop over all spins
 do y=0,ll-1
    do x=0,ll-1
        !determine index of spins array
        ii = ll*y+x
        xp1 = mod(x+1,ll)
        yp1 = mod(y+1,ll)
        xm1 = mod(x+ll-1,ll)
        ym1 = mod(y+ll-1,ll)
        neighbors = spin(ll*y+xp1) + spin(ll*y+xm1) + spin(ll*yp1+x) + spin(ll*ym1+x)
        census(neighbors+(spin(ii)+2)/2) = census(neighbors+(spin(ii)+2)/2) + 1
    enddo
enddo
end subroutine doCensus
!---------------------!

!----------------------------------------------!
! Prints the configuration to stdout !
!----------------------------------------------!
 subroutine printConfig()
!-------------------!
 use systemvariables
 implicit none

 integer :: y,x

 do y=-1,ll-1
    do x=0,ll-1
        !write out line separating
        if (y==-1) then
            write(*,'(a)',advance='no') '-'
        else
            if (spin(ll*y+x) == 1) then
                write(*,'(a)',advance='no') "0"
            else
                write(*,'(a)',advance='no') "."
            end if
        end if
    end do
    write(*,*) ""
 end do

 end subroutine printConfig
!---------------------!

!----------------------------------------------!
! Carries out one Monte Carlo step, defined as !
! nn flip attempts of randomly selected spins. !
!----------------------------------------------!
 subroutine mcstep()
!-------------------!
 use systemvariables
 implicit none

 integer :: i,s,x,y,s1,s2,s3,s4
 real(8) :: temp
 real(8) :: crand

 do i=1,nn
    s=int(crand()*nn)               ! generate random site
    x=mod(s,ll); y=s/ll           ! coordinates of site s
    s1=spin(mod(x+1,ll)+y*ll)     ! spin at right neighbor of site s
    s2=spin(x+mod(y+1,ll)*ll)     ! spin at up    neighbor of site s
    s3=spin(mod(x-1+ll,ll)+y*ll)  ! spin at left  neighbor of site s
    s4=spin(x+mod(y-1+ll,ll)*ll)  ! spin at down  neighbor of site s
    temp = crand()
    if (temp<pflip(spin(s)*(s1+s2+s3+s4),spin(s))) then
        spin(s)=-spin(s)  ! flip spin with Metroplois prob
    end if
 enddo

 end subroutine mcstep
!---------------------!

!----------------------------------------------!
! calculates total probability of all possible spin flips
!----------------------------------------------!
real(8) function sumFlips()
!-------------------!
    use systemvariables
    implicit none
    integer :: s,x,y,s1,s2,s3,s4

    !calc total probability of flipping any spin
    sumFlips = 0.d0
    do s=1,nn
        x=mod(s,ll); y=s/ll           ! coordinates of site s
        s1=spin(mod(x+1,ll)+y*ll)     ! spin at right neighbor of site s
        s2=spin(x+mod(y+1,ll)*ll)     ! spin at up    neighbor of site s
        s3=spin(mod(x-1+ll,ll)+y*ll)  ! spin at left  neighbor of site s
        s4=spin(x+mod(y-1+ll,ll)*ll)  ! spin at down  neighbor of site s
        sumFlips = sumFlips + min(1.d0,pflip(spin(s)*(s1+s2+s3+s4),spin(s)))
    enddo
    
 end function sumFlips
!---------------------!

!----------------------!
 real(8) function crand()
!-----------------------------------------------------!
! 64-bit linear congruental random number generator   !
! iran64=oran64*2862933555777941757+1013904243        !
!-----------------------------------------------------!
 implicit none

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 ran64=ran64*mul64+add64
 crand=0.5d0+dmu64*dble(ran64)
 
 end function crand
!----------------!

!---------------------!
 subroutine initran(w)
!---------------------!
 implicit none

 integer(8) :: irmax
 integer(4) :: w,nb,b

 real(8)    :: dmu64
 integer(8) :: ran64,mul64,add64
 common/bran64/dmu64,ran64,mul64,add64

 irmax=2_8**31
 irmax=2*(irmax**2-1)+1
 mul64=2862933555777941757_8
 add64=1013904243
 dmu64=0.5d0/dble(irmax)

 open(10,file='seed.in',status='old')
 read(10,*)ran64
 close(10)
 if (w.ne.0) then
    open(10,file='seed.in',status='unknown')
    write(10,*)abs((ran64*mul64)/5+5265361)
    close(10)
 endif

 end subroutine initran
!----------------------!
