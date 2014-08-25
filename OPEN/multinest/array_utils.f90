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


end module array_utils