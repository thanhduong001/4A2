      
      subroutine generate_mesh(geom,g)

!     Create cells of the mesh to cover the domain defined by geometry curves,
!     the values of the node coordinates, x(i,j) and y(i,j) are stored in "g"

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_geometry), intent(in) :: geom
      type(t_grid), intent(inout) :: g
      real :: si_a(geom%ni_a), si_b(geom%ni_b), si(g%ni), sj(g%nj)
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT
      real :: xstart, xend, ystart, yend
      ! integer :: i, j

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the non-dimensional curve lengths of the geometry input and
!     generate linearly spaced points in i-direction at desired mesh resolution
      call dist(geom%x_a,geom%y_a,1,si_a)
      call dist(geom%x_b,geom%y_b,1,si_b)
      call linspace(0.0,1.0,si)
      call linspace(0.0,1.0,sj)

!     Interpolate the geometry curves to the required resolution in the 
!     i-direction, this allows the mesh to be refined without altering the 
!     geometry definition file, the data is stored at "j = 1" and "j = nj"
      call interp(si_a,geom%x_a,si,g%x(:,1))
      call interp(si_a,geom%y_a,si,g%y(:,1))
      call interp(si_b,geom%x_b,si,g%x(:,nj))
      call interp(si_b,geom%y_b,si,g%y(:,nj))

!     Calculate the coordinates of all the intermediate points within the mesh.
!     Create a new vector of non-dimensional spacings in the j-direction using 
!     "linspace", loop over the mesh in the i-direction and calculate the
!     intermediate coordinates from a weighted sum of the two boundaries
!     INSERT
     do i = i, ni
         xstart = g%x(i,1)
         xend = g%x(i,nj)
         ystart = g%y(i,1)
         yend = g%y(i,nj)
         call linspace(xstart, xend , g%x(i,:))
         call linspace(ystart, yend, g%y(i,:))
     end do
      ! do i = 1, ni
        
      !         g%x(i,:) = sj*g%x(i,1) + (1-sj)*g%x(i,nj)
      !         g%y(i,:) = sj*g%y(i,1) + (1-sj)*g%y(i,nj)
         
      !  end do
!      do i = 1, ni
!          call linspace(g%x(i,1), g%x(i,nj) , g%x(i,:))
!          call linspace(g%y(i,1), g%y(i,nj), g%y(i,:))
!      end do
!     In all of the test cases for the basic solver the the "j = 1" and "j = nj"
!     boundaries are walls, for the extensions you may need to return to this
!     and communicate the position of the walls to the solver in a more 
!     flexible way. The "wall" variable is an "ni * nj" logical array where 
!     "true" indicates the node is on a wall.
      g%wall = .false.
      g%wall(:,[1,g%nj]) = .true.

!     Print that the mesh has been created
      write(6,*) 'Interpolated mesh from the bounding geometry curves'
      write(6,*)

      end subroutine generate_mesh


