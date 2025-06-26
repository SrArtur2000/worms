program isingworm
	implicit none
	integer(8), parameter :: L = 20
	integer, parameter :: iseed1 = 10
	integer(8), parameter :: Nth = 5000, N = L**2, Na = 2*N
	
	real(8) :: T,qsi1,qsi2,E1,w,C1,C2,lmean,l2mean,beta,E2,aux,soma,lmean1,l1,l2,lmean2
	integer(8) :: j,k,k1,k2,mu,n1,n2,nu,Nb,lm2
	integer(8), dimension(0:1,2) :: snake
	integer(8), dimension(4,L,L) :: A
	integer(8), dimension(4,L) :: nn
	integer(8) :: Z,lworm,l3,l4
	
!	call random_init(.true.,.true.)
	
	open(10,file='qsi1.dat')
	open(11,file='qsi2.dat')
	open(20,file='E1.dat')
	open(21,file='E2.dat')
	open(30,file='C1.dat')
	open(31,file='C2.dat')
	open(40,file='Nb.dat')
	
	call headmettail(snake)

	do j = 1, L
	  nn(1,j) = 1 + mod(j,L)
	  nn(2,j) = 1 + mod(j,L)
	  nn(3,j) = L - mod(L-j+1,L)
	  nn(4,j) = L - mod(L-j+1,L)
	enddo
	
	T = 0.1d0
	do while (T .lt. 5.0d0)
	  beta = 1.0d0/T
	  w = tanh(beta)
	  qsi1 = 0.0d0
	  qsi2 = 0.0d0
	  lm2 = 0.0d0
	  l1 = 0.0d0
	  l2 = 0.0d0
	  l3 = 0
	  l4 = 0
	  Z = 0
	  A = 0
	  lworm = 0
	  soma = 0.0d0
	  do j = 1, Nth
	    do k = 1, N
	      call random_number(aux)
	      mu = 4*aux+1
	      call movepen(A,mu,snake,nn,l1,l2,l3,l4,lworm,Z,qsi2,beta,w,lm2,E1,soma)
	    enddo
	  enddo
	  do j = 1, Nth
	    do k = 1, N
	      call random_number(aux)
	      mu = 4*aux+1
!	      print*, mu
	      call movepen(A,mu,snake,nn,l1,l2,l3,l4,lworm,Z,qsi2,beta,w,lm2,E1,soma)
	    enddo
	    lmean = real(l1)/Z
	    if(T .eq. 0.1d0) then
	      write(2,*)j, lmean
	      write(3,*)j, Z
	    endif
	  enddo
	  print*, "Terminou a temperatura", T
	  if(Z .ne. 0) then
	    lmean = l1/soma
	    lmean1 = real(l3)/Z
	    l2mean = l2/soma
	    lmean2 = real(l4)/Z
	    E1 = lmean1*(-1.0d0/(sinh(beta)*cosh(beta)))
	    E1 = E1-Na*w
	    E1 = E1/N
	    E2 = (-Na*w-(1.0d0/w-w)*lmean1)/N
!	    E2 = real(l1)/(Z*L/2)
	    print*, lmean, l2mean, sinh(beta), lmean**2, cosh(beta), beta*w
	    C1 = (((beta*w)**2)/((sinh(beta))**4))*(2.0d0*N*((sinh(beta))**2)+lmean+l2mean-(lmean)**2 - 2*lmean*(cosh(beta))**2)/N
	    C2 = (beta**2)*((1.0d0/((sinh(beta)*cosh(beta))**2))*(lmean2 - lmean1 - lmean1**2)&
	    + (1.0d0/((cosh(beta))**2))*(Na - 2.0d0*lmean1))/N
!	    write(10,*)T, qsi1
	    write(11,*)T, qsi2
	    write(20,*)T, abs(E1-E2)
	    write(21,*)T, E2
	    write(30,*)T, abs(C1-C2)
	    write(31,*)T, C2
	    write(40,*)T, lmean, lmean1
	  endif
!	  call printsnake(A,nn,T)
	  T = T + 0.01d0
	enddo
contains

subroutine movepen(A,mu,snake,nn,l1,l2,l3,l4,lworm,Z,qsi,beta,w,lm2,E,soma)
	implicit none
	
	integer(8) :: mu,x,dir,Z,mmu,lworm,dif,lm2,l3,l4
	integer(8), dimension(0:1,2) :: snake, snakeaux
	real(8) :: r,qsi,E,beta,w,soma,l1,l2
	integer(8), dimension(4,L,L) :: A
	integer(8), dimension(4,L) :: nn
	!Move a cabeça ou o rabo da cobra numa direção mu aleatória
	call movesnake(A,mu,lworm,snake,nn,lm2,qsi,beta,w)
	!Após mover, checa se a cabeça e o rabo estão no mesmo ponto
!	print*, lworm
	if((snake(0,1) .eq. snake(1,1)) .and. (snake(0,2) .eq. snake(1,2)) ) then
	  call random_number(r)
	  Z = Z + 1
	  l1 = l1 + lworm*tanh(beta)**(lworm)
	  l2 = l2 + (lworm**2)*(tanh(beta))**(lworm)
	  l3 = l3 + lworm
	  l4 = l4 + lworm**2
	  soma = soma + (tanh(beta))**lworm
!	  lworm = 0
!	  A = 0
	  if(r .le. 0.5d0) then
	    !passo 1, remove cobra
	    call headmettail(snake)
!	    call movesnake(A,mu,lworm,snake,nn,lm2,qsi,beta,w)
	  endif
	endif
end subroutine

subroutine headmettail(snake)
	implicit none
	
	integer(8), dimension(0:1,2) :: snake
	integer :: n1,n2
	real :: a1,a2

	call random_number(a1)
	call random_number(a2)
	n1 = L*a1 + 1
	n2 = L*a2 + 1
	snake(0,:) = (/n1,n2/)
	snake(1,:) = (/n1,n2/)
end subroutine

subroutine movesnake(A,mu,lworm,snake,nn,lm2,qsi,beta,w)
	implicit none
	
	integer(8), dimension(0:1,2) :: snake, snakeaux
	integer(8) :: mu,x,dir,Z,mmu,lworm,lm2,dif
	integer(8), dimension(4,L,L) :: A
	integer(8), dimension(4,L) :: nn
	real(8), dimension(2) :: r
	real(8) :: qsi,E,beta,w

	w = tanh(beta)
	dir = mod(mu,2) + (1+(-1)**mu)
	! Simplified calculation for opposite direction
	mmu = mod(mu+1, 4) + 1
	call random_number(r)
	r(1)  = 2.0d0*r(1)
	x = 1
	if(r(2) .le. w**(1-A(mu,snake(x,1),snake(x,2))) ) then
!	  print*, "O ponto inicial é"
!	  print*, mu, snake(x,1), snake(x,2)
	  dif = mod(A(mu,snake(x,1),snake(x,2))+1,2) - A(mu,snake(x,1),snake(x,2))
	  A(mu,snake(x,1),snake(x,2)) = mod(A(mu,snake(x,1),snake(x,2))+1,2)
	  snakeaux = snake
	  snake(x,dir) = nn(mu,snake(x,dir))
!	  print*, "O ponto final é"
!	  print*, mmu, snake(x,1), snake(x,2)
	  A(mmu,snake(x,1),snake(x,2)) = mod(A(mu,snakeaux(x,1),snakeaux(x,2))+1,2)
	  lworm = lworm + dif
	  qsi = qsi + dif
	  lm2 = lm2 + 1
!	    write(1000+int(T*1000),*)snake(x,1),snake(x,2),mod(snake(x,1)-snakeaux(x,1),19),mod(snake(x,2)-snakeaux(x,2),19)
	endif
end subroutine

subroutine printsnake(A,nn,T)
	implicit none

	integer(8), dimension(4,L,L) :: A
	integer(8), dimension(4,L) :: nn
	integer(8), dimension(2) :: pos
	integer(8) :: mu1,j1,k1
	real(8) :: T

	do j1 = 1,L
	  do k1 = 1,L
	    do mu1 = 1,2
	      if (A(mu1,j1,k1) .ne. 0) then
	        pos = (/j1,k1/)
	        pos(mu1) = pos(mu1) + 1
	        write(100+int(T*100),*)j1,k1,pos(1)-j1, pos(2)-k1
!	        if(j1 .eq. L .and. mu1 .eq. 1) then
!	          write(100+int(T*100),*)j1,k1,pos(1)-j1, pos(2)-k1
!	          write(100+int(T*100),*)0,k1,pos(1)-j1, pos(2)-k1
!	        elseif (k1 .eq. L .and. mu1 .eq. 2) then
!	          write(100+int(T*100),*)j1,k1,pos(1)-j1, pos(2)-k1
!	          write(100+int(T*100),*)j1,0,pos(1)-j1, pos(2)-k1
!	        elseif(j1 .eq. 0 .and. mu .eq. 1) then
!		    write(100+int(T*100),*)L,k1,pos(1)-j1, pos(2)-k1
!	          write(100+int(T*100),*)j1,k1,pos(1)-j1, pos(2)-k1
!	        elseif(k1 .eq. 0 .and. mu .eq. 2) then
!	          write(100+int(T*100),*)j1,k1,pos(1)-j1, pos(2)-k1
!	          write(100+int(T*100),*)j1,L,pos(1)-j1, pos(2)-k1
!	        endif
	      endif
	    enddo
	  enddo
	enddo
	write(100+int(T*100),*) snake(0,1),snake(0,2) 
	write(100+int(T*100),*) snake(1,1),snake(1,2) 
end subroutine

end program
