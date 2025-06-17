!      T.M., 2012 -> 2024
!
!      2d Ising Model using the Worm Algorithm
!
!      NOTE: update still very inefficient, should compute weights at the
!      beginning, not at every call to movepen...
!
       program ising
       implicit real*8 (a-h,o-z)
       implicit integer (i-n)
       dimension nnu(100), nnd(100)          ! n.n.s (up/down+right/left)
       dimension kbu(100,100), kbr(100,100)  ! bonds (up+right) 2 per site
	 
	 open(10,file='worm.dat')
       ! input
       read (1,*) Nth           ! number of iterations for thermalization
       read (1,*) N             ! number of iterations
       read (1,*) L             ! lattice side
       read (1,*) beta          ! in practice, beta*J
	 
       ! initialization: 
	 
       ! determine nearest neighbors (imposing periodic b.c.)
       do i = 1, L
         nnu(i) = 1 + mod(i,L)
         nnd(i) = L - mod(L-i+1,L)
       enddo
	 
       ! all empty bonds 
       do j = 1, L
         do i = 1, L
           kbu(i,j) = 0
           kbr(i,j) = 0
         enddo
       enddo
	 
       ! ira=masha: I = M = (1, 1)
       Ix = 1
       Iy = 1
       Mx = 1
       My = 1
	 
       lenw = 0         ! sum of all kb ("length" of the worm)
       L2 = L*L         ! volume of the lattice
       Lm1 = L-1        ! to check for n.n.
	 
       call srand(212)
       
       ! thermalization (each iteration "visits" all bonds)
       do kiter = 1, Nth*2*L2
         call update(Ix,Iy,Mx,My,L,nnu,nnd,kbu,kbr,lenw,beta)
       enddo
	 
       print*, lenw
	 
       ! set/reset counters to zero
       lZ = 0          ! Z estimator
       lM2 = 0         ! squared magnetization estimator = sum of all g(r)
       lg1 = 0         ! energy estimator (beta*V*g(1)/2)*1/Z
       lE = 0          ! alternative estimator (beta*sum of all kb, see PRL 2001)
	 
       ! initialize averages
       Nef = N
       avmag = 0.d0
       avmag2 = 0.d0
       avener = 0.d0
       avener2 = 0.d0
	 
       ! main loop
	 
       do kiter = 1, N
         !print*,kiter
         do jiter = 1, 2*L2
           call update(Ix,Iy,Mx,My,L,nnu,nnd,kbu,kbr,lenw,beta)
           call measure(Ix,Iy,Mx,My,Lm1,lenw,lZ,lg1,lE,lM2)
         enddo
         !print*,lenw
         if (lZ .ge. 1) then
           ! write output to fort.2
           smag = sqrt( ( dfloat(lM2+1) / dfloat(lZ) ) / dfloat(L2) )
           sener = - beta *  dfloat(lg1) / ( 2.d0 * dfloat(lZ) )
           sener2 = - beta * ( dfloat(lE) / dfloat(lZ) ) / dfloat(L2)
           write (2,*) kiter,lZ,lM2,lg1,lE,smag,sener,sener2
           !write (6,*) kiter,lZ,lM2,lg1,lE,smag,sener,sener2
           ! update averages
           avmag = avmag + smag
           avmag2 = avmag2 + smag**2
           avener = avener + sener
           avener2 = avener2 + sener**2
         else
           Nef = Nef - 1
         endif
       enddo
	 
       print*,"Nef",Nef
       avmag = avmag / ( dfloat(Nef) )
       avmag2 = sqrt( ( avmag2 / dfloat(Nef) ) - avmag**2 )
       write (6,*) "MAG", avmag, avmag2
       avener = avener / ( dfloat(Nef) )
       avener2 = sqrt( ( avener2 / dfloat(Nef) ) - avener**2 )
       write (6,*) "ENER", avener, avener2
	 
       end program
	 
	 
       subroutine update(Ix,Iy,Mx,My,L,nnu,nnd,kbu,kbr,lenw,beta)
       implicit real*8 (a-h,o-z)
       implicit integer (i-n)
       dimension nnu(100), nnd(100)                ! the neighbors
       dimension kbu(100,100), kbr(100,100)        ! the bonds
	 
       ! decide what move to make
       if (Ix .eq. Mx .and. Iy .eq. My) then
           ! move I = M somewhere else
           Ix = 1 + int( L * rand() )     ! 1 + integer in [0, L-1]
           Iy = 1 + int( L * rand() )     ! 1 + integer in [0, L-1]
           Mx = Ix
           My = Iy
           !print*, "moving masha-ira to", Ix, Iy
       endif
	 
       !print*, "let's shift masha"
       call movepen(Mx,My,nnu,nnd,kbu,kbr,lenw,beta)
	 
       return
       end
	 
	 
       subroutine measure(Ix,Iy,Mx,My,Lm1,lenw,lZ,lg1,lE,lM2)
       implicit real*8 (a-h,o-z)
       implicit integer (i-n)
       dimension kbu(100,100), kbr(100,100)              ! the bonds
	 
       ! check for CP (contrib to Z and E) and n.n. (contrib to g1)
       if (Ix .eq. Mx) then
         if (Iy .eq. My) then
           !print*, "masha=ira"
           ! count the "new" diagram as Z-type
           lZ = lZ + 1
           lE = lE + lenw
         else 
           idif = iabs(Iy-My)
           if ( idif .eq. 1 .or. idif .eq. Lm1 ) then
             lg1 = lg1 + 1
           endif
         endif
       endif
	 
       if (Iy .eq. My) then
         idif = iabs(Ix-Mx)
         if (idif .eq. 1 .or. dif .eq. Lm1) then
           lg1 = lg1 + 1
         endif
       endif

       ! count new diagram as g-type (including G(0))
       lM2 = lM2 + 1

       return
       end

       subroutine movepen(Mx,My,nnu,nnd,kbu,kbr,lenw,beta)
       implicit real*8 (a-h,o-z)
       implicit integer (i-n)
       dimension nnu(100), nnd(100)                      ! the neighbors
       dimension kbu(100,100), kbr(100,100)              ! the bonds

       ! decide if bonds will be preferably drawn or erased
       tb = tanh(beta)
       if (tb .gt. 1.d0) then
         iaux = 0                          ! drawn if empty
         paux = 1.d0                       ! otherwise erased with prob. paux
       else
         iaux = 1                          ! erased if present
         paux = tb                         ! otherwise drawn with prob. paux
       endif

       ! choose direction to move pen (AND draw or erase, accordingly)
       r = rand()
       if (r .lt. 0.25d0) then
         ! move up (draw/erase vertical bond belonging to Mx,My)
         call chbond(kbu(Mx,My),iaux,paux,lenw)
         My = nnu(My)
       else
         if (r .lt. 0.5d0) then
           ! move down (draw/erase vertical bond belonging to Mx,My-1)
           My = nnd(My)
           call chbond(kbu(Mx,My),iaux,paux,lenw)
         else
           if (r .lt. 0.75d0) then
             ! move right (draw/erase horizontal bond belonging to Mx,My)
             call chbond(kbr(Mx,My),iaux,paux,lenw)
             Mx = nnu(Mx)
           else
             ! move left (draw/erase horizontal bond belonging to Mx-1,My)
             Mx = nnd(Mx)
             call chbond(kbr(Mx,My),iaux,paux,lenw)
           endif
         endif
       endif

       return
       end


       subroutine chbond(n,i,p,lenw)
       implicit real*8 (a-h,o-z)
       implicit integer (i-n)

       if (i .eq. 0) then                ! preferably drawn
         if (n .eq. 0) then
           n = 1
           lenw = lenw + 1
         else
           raux = rand()
           if (raux .lt. p) then
             n = 0
             lenw = lenw - 1
           endif
         endif
       endif

       if (i .eq. 1) then                ! preferably erased
         if (n .eq. 1) then
           n = 0
           lenw = lenw - 1
         else
           raux = rand()
           if (raux .lt. p) then
             n = 1
             lenw = lenw + 1
           endif
         endif
       endif

       end 
