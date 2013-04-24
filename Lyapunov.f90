program Lyapunov
implicit none
	
	real(8) :: xx, yy, zz ,kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4, kz1, kz2, kz3, kz4, h	!RK4 to solve the function x'=f(x) 
	real(8), dimension(3,3) :: ko1, ko2, ko3, ko4, o	!RK4 to solve the function T'=F(T)
	real(8), dimension(3,3) :: DF	!Jacobian Matrix of the system function
	integer :: i, j, ijk, j0, j1, j2, it	!counter
	real(8), external :: f	!function to calculate
	real(8), parameter :: a = 10., b = 28., c = 2.	!parameter of the function
	real(8), dimension(3) :: lamda, temp_lamda, lyn, cc, dd
	real(8), dimension(3,3) :: ot, QQ, RR, Rt, Qt
	integer, parameter :: ni = 3
	LOGICAL sing

	
	call random_number(xx)
	call random_number(yy)
	call random_number(zz)

	open (10, file='LorenzLyapunov.dat')
	
	do i = 1, 3	
		o(i, i) = 1
		QQ(i, i) = 1
	end do

	h=0.001
	
	DF(1, 1) = -a
	DF(1, 2) = a
	DF(1, 3) = 0
	DF(2, 2) = -1
	DF(3, 3) = -c
	
	ijk = 0
		
	do i = 1, 10000000

		kx1 = f(1, xx, yy, zz, a, b, c)
		ky1 = f(2, xx, yy, zz, a, b, c)
		kz1 = f(3, xx, yy, zz, a, b, c)	
		if (i > 5000000) then		
			DF(2, 1) = b - zz
			DF(2, 3) = -xx
			DF(3, 1) = yy
			DF(3, 2) = xx
			ko1 = matmul(DF, o)
		end if
		
		kx2 = f(1, xx + kx1 * h / 2, yy + ky1 * h / 2, zz + kz1 * h / 2, a, b, c)
		ky2 = f(2, xx + kx1 * h / 2, yy + ky1 * h / 2, zz + kz1 * h / 2, a, b, c)
		kz2 = f(3, xx + kx1 * h / 2, yy + ky1 * h / 2, zz + kz1 * h / 2, a, b, c)
		if (i > 5000000) then	
			DF(2, 1) = b - (zz + kz1 * h / 2)
			DF(2, 3) = -(xx + kx1 * h / 2)
			DF(3, 1) = yy + ky1 * h / 2
			DF(3, 2) = xx + kx1 * h / 2
			ko2 = matmul(DF, o + ko1 * h / 2)
		end if 
		
		kx3 = f(1, xx + kx2 * h / 2, yy + ky2 * h / 2, zz + kz2 * h / 2, a, b, c)
		ky3 = f(2, xx + kx2 * h / 2, yy + ky2 * h / 2, zz + kz2 * h / 2, a, b, c)
		kz3 = f(3, xx + kx2 * h / 2, yy + ky2 * h / 2, zz + kz2 * h / 2, a, b, c)
		if (i > 5000000) then		
			DF(2, 1) = b - (zz + kz2 * h / 2)
			DF(2, 3) = -(xx + kx2 * h / 2)
			DF(3, 1) = yy + ky2 * h / 2
			DF(3, 2) = xx + kx2 * h / 2
			ko3 = matmul(DF, o + ko2 * h / 2)
		end if
		
		kx4 = f(1, xx + kx3 * h, yy + ky3 * h, zz + kz3 * h, a, b, c)
		ky4 = f(2, xx + kx3 * h, yy + ky3 * h, zz + kz3 * h, a, b, c)
		kz4 = f(3, xx + kx3 * h, yy + ky3 * h, zz + kz3 * h, a, b, c)
		if (i > 5000000) then	
			DF(2, 1) = b - (zz + kz3 * h)
			DF(2, 3) = -(xx + kx3 * h)
			DF(3, 1) = yy + ky3 * h
			DF(3, 2) = xx + kx3 * h
			ko4 = matmul(DF, o + ko3 * h)
		end if

		xx = xx + h * (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6
		yy = yy + h * (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6
		zz = zz + h * (kz1 + 2 * kz2 + 2 * kz3 + kz4) / 6

		if (i > 5000000) then
			ijk = ijk + 1
			o = o + h * (ko1 + 2 * ko2 + 2 * ko3 + ko4) / 6
			
			Ot = matmul(o, QQ)
			call ddqrdcmp(Ot,ni,ni,cc,dd,sing)
			
			do j=1,ni-1
				Rt(j,j)=dd(j)
				do j1=j+1,ni
					Rt(j,j1)=Ot(j,j1)
					Ot(j,j1)=0
				end do
			end do
			j=ni
			Rt(j,j)=dd(j)
			do j=1,ni
				do j1=1,ni
					QQ(j1,j)=0
				end do
				QQ(j,j)=1
			end do

			do j=1,ni-1
				do j1=1,ni
					do j2=1,ni
						Qt(j1,j2)=-Ot(j1,j)*Ot(j2,j)/cc(j)
						if (j1 == j2) Qt(j1,j2)=Qt(j1,j2)+1
					end do
				end do
					
				QQ = matmul(QQ, Qt)
			end do

			do j=1,ni
				RR(j,j)=RR(j,j)+log(abs(Rt(j,j)))
			end do

			do it=1,ni
				do j=1,ni
					O(j,it)=0
				end do
			O(it,it)=1
			end do


			do j=1,3
			  lyn(j)=RR(j,j)
			end do

			if (mod(ijk, 10000) == 0)then
				write(10, "(101e20.12)") i * h, lyn/ijk/h
			end if

			if (mod(ijk,100000) == 0) then			
				write (*,*) i * h, lyn/ijk/h
			end if

			
		end if
		
	end do
	
		
end program





	function f(j, x, y, z, a, b, c)
	implicit none
		real(8) :: f
		real(8) :: x, y, z
		real(8) :: a, b, c
		integer :: j
		if (j == 1) then
			f = a * (y - x)
		else if (j == 2) then
			f = x * (b - z) - y
		else if (j == 3) then
			f =x * y - c * z
		end if
	end function f
	
	
	SUBROUTINE ddqrdcmp(a,n,np,c,d,sing)
	implicit none
		INTEGER n,np
		REAL(8) a(np,np),c(n),d(n)
		LOGICAL sing
		INTEGER i,j,k
		REAL(8) scale,sigma,sum,tau

		sing=.false.

		do k=1,n-1
			scale=0.
			do i=k,n
				scale=max(scale,abs(a(i,k)))
			end do
			if(scale==0.)then
				sing=.true.
				c(k)=0.
				d(k)=0.
			else
				do i=k,n
					a(i,k)=a(i,k)/scale
				end do
				sum=0.
				do i=k,n
					sum=sum+a(i,k)**2
				end do
				sigma=sign(sqrt(sum),a(k,k))
				a(k,k)=a(k,k)+sigma
				c(k)=sigma*a(k,k)
				d(k)=-scale*sigma
				do j=k+1,n
					sum=0.
					do i=k,n
						sum=sum+a(i,k)*a(i,j)
					end do
					tau=sum/c(k)
					do i=k,n
						a(i,j)=a(i,j)-tau*a(i,k)
					end do
				end do
			end if
		end do
		
		d(n)=a(n,n)
		if(d(n).eq.0.)sing=.true.
		
	END