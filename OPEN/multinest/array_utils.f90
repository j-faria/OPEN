module array_utils

	implicit none

	interface ones_like
		procedure ones_like_1d
		procedure ones_like_2d
	end interface ones_like

contains


subroutine fill_matrix_with_value(A, val)
	! Set all the elements of matrix A[m x n] to be val
	implicit none
	real(kind=8), dimension(:,:), intent(inout) :: A
	real(kind=8), intent(in) :: val

	A = val

end subroutine fill_matrix_with_value

subroutine set_diagonal_to(A, val)
	! Set diagonal elements of matrix A[m x n] to be val
	implicit none
	real(kind=8), dimension(:,:), intent(inout) :: A
	real(kind=8), intent(in) :: val

	! local variable
	integer :: i 

	FORALL (i=1:minval(shape(A)))  A(i, i) = val

end subroutine set_diagonal_to


subroutine add_value_to_diagonal(A, val)
	! Add value val to diagonal elements of matrix A[m x n]
	implicit none
	real(kind=8), dimension(:,:), intent(inout) :: A
	real(kind=8), intent(in) :: val

	! local variable
	integer :: i 

	FORALL (i=1:minval(shape(A)))  A(i,i) = A(i,i) + val

end subroutine add_value_to_diagonal


subroutine add_array_to_diagonal(A, arr)
	! Add value val to diagonal elements of matrix A[m x n]
	implicit none
	real(kind=8), dimension(:,:), intent(inout) :: A
	real(kind=8), dimension(:), intent(in) :: arr

	! local variable
	integer :: i 

	if (size(arr) /= minval(shape(A))) stop 'Size mismatch in add_array_to_diagonal'

	FORALL (i=1:minval(shape(A)))  A(i,i) = A(i,i) + arr(i)

end subroutine add_array_to_diagonal



function get_diagonal(A) result(diag)
	! Get diagonal elements of matrix A (general [m x n] dimensions)
	implicit none
	real(kind=8), dimension(:,:), intent(in) :: A
	real(kind=8), dimension(minval(shape(a))) :: diag
	! local variables
	integer :: i

	FORALL (i=1:size(diag)) diag(i) = A(i,i)

end function get_diagonal


function eye(n) result(A)
	! Create an [n x n] identity matrix
	implicit none
	integer, intent(in) :: n
	real(kind=8), dimension(n,n) :: A

	A = 0.d0
	call set_diagonal_to(A, 1.d0)
	return
	
end function eye


function ones_like_1d(arr) result(ones)
	implicit none
	real(kind=8), dimension(:), intent(in) :: arr
	real(kind=8), dimension(:), allocatable :: ones

	allocate(ones(size(arr)))
	ones = 1.d0
end function ones_like_1d

function ones_like_2d(arr) result(ones)
	implicit none
	real(kind=8), dimension(:,:), intent(in) :: arr
	real(kind=8), dimension(:,:), allocatable :: ones
	integer, dimension(2) :: dims

	dims = shape(arr)
	allocate(ones(dims(1), dims(2)))
	ones = 1.d0
end function ones_like_2d


function linspace(xmin, xmax, num) result(x)
    implicit none
    real(kind=8), intent(in) :: xmin, xmax
    integer, intent(in) :: num
    real(kind=8), dimension(num) :: x
    integer :: i
    if (num == 1) then
       if(xmin /= xmax) then
          write(0,'("ERROR: Cannot call linspace with n=1 and xmin /= xmax")')
          stop
       else
          x = xmin
       end if
    else
       do i=1,num
          x(i) = (xmax-xmin) * real(i-1, kind=8) / real(num-1, kind=8) + xmin
       end do
    end if
end function linspace

subroutine sort(A, N)
	! arguments
	integer, intent(in) :: N
	real(kind=8), dimension(N), intent(inout) :: A
	! local variables
	real(kind=8), dimension((N+1)/2) :: T

	call MergeSort(A,N,T)

	contains 
		subroutine Merge(A,NA,B,NB,C,NC)
		
			integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
			real(kind=8), intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
			real(kind=8), intent(in)     :: B(NB)
			real(kind=8), intent(in out) :: C(NC)

			integer :: I,J,K
		 
			I = 1; J = 1; K = 1;
			do while(I <= NA .and. J <= NB)
				if (A(I) <= B(J)) then
					C(K) = A(I)
					I = I+1
				else
			 		C(K) = B(J)
					J = J+1
				endif
				K = K + 1
			enddo
			do while (I <= NA)
				C(K) = A(I)
				I = I + 1
				K = K + 1
			enddo
			return
		end subroutine merge
 
		recursive subroutine MergeSort(A,N,T)
			integer, intent(in) :: N
			real(kind=8), dimension(N), intent(in out) :: A
			real(kind=8), dimension((N+1)/2), intent (out) :: T

			integer :: NA, NB
			real(kind=8) :: V
		 
			if (N < 2) return
			if (N == 2) then
				if (A(1) > A(2)) then
					V = A(1)
					A(1) = A(2)
					A(2) = V
				endif
				return
			endif      
			NA=(N+1)/2
			NB=N-NA

			call MergeSort(A,NA,T)
			call MergeSort(A(NA+1),NB,T)

			if (A(NA) > A(NA+1)) then
				T(1:NA)=A(1:NA)
				call Merge(T,NA,A(NA+1),NB,A,N)
			endif
			return
		end subroutine MergeSort
 
end subroutine sort

subroutine swap(a, b)
	implicit none
	real(kind=8), dimension(:), intent(inout) :: a, b
	real(kind=8), dimension(size(a)) :: temp
	if (size(a) /= size(b)) STOP 'Error in swap; different dimensions'
	temp = a ; a = b ; b = temp
end subroutine swap

function variance(arr)
	implicit none
	real(kind=8), intent(in), dimension(:) :: arr
	real(kind=8) :: variance, x
	x = sum(arr) / size(arr)
	variance = sum((arr-x)**2)/size(arr)
end function variance


end module array_utils