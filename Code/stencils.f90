
      module stencils

!     Packaging a subroutine in a module allows it to recieve the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sum_fluxes(av,flux_i,flux_j,area,prop,dcell)

!     This subroutine sums the fluxes into each cell, calculates the change in 
!     the cell property inside, distributes the change to the four nodes of the
!     cell and then adds it onto the flow property

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: flux_i(:,:), flux_j(:,:), area(:,:)
      real, intent(inout) :: prop(:,:)
      real, intent(out) :: dcell(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: dnode
      real, dimension(size(prop,1)-1,size(prop,2)-1) :: dprop_cell
      integer :: ni, nj

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Use the finite volume method to find the change in the variables "prop"
!     over the timestep "dt", save it in the array "dprop_cell"
!     INSERT
      dprop_cell = (flux_i(1:ni-1,:) &
                  - flux_i(2:ni,:) &
                  + flux_j(:,1:nj-1) &
                  - flux_j(:,2:nj)) * av%dt

!     Now distribute the changes equally to the four corners of each cell. Each 
!     interior grid point receives one quarter of the change from each of the 
!     four cells adjacent to it.
!     INSERT (internal nodes)
      dnode(2:ni-1,2:nj-1) = 0.25 * (dprop_cell(1:ni-2,1:nj-2) + &
                                    dprop_cell(2:ni-1,1:nj-2) + &
                                    dprop_cell(1:nj-2,2:ni-1) + &
                                    dprop_cell(2:ni-1,2:nj-1))

!     Bounding edge nodes do not have four adjacent cells and so must be treated
!     differently, they only recieve half the change from each of the two
!     adjacent cells. Distribute the changes for the "i = 1 & ni" edges as well
!     as the "j = 1 & nj" edges. 
!     INSERT (external nodes)
      dnode(1,:) = 0.5 * (dprop_cell(1,1:nj-2) + &
                        dprop_cell(1,2:nj-1))
      dnode(ni,:) = 0.5 * (dprop_cell(ni,1:nj-2) + &
                        dprop_cell(ni,2:nj-1))
      dnode(:,nj) = 0.5 * (dprop_cell(1:ni-2,nj) + &
                        dprop_cell(2:ni-1,nj))
      dnode(:,1) = 0.5 * (dprop_cell(1:ni-2,1) + &
                        dprop_cell(2:ni-1,1))

!     Finally distribute the changes to be to the four bounding corner points, 
!     these receive the full change from the single cell of which they form one 
!     corner.
!     INSERT
      dnode(1,1) = dprop_cell(1,1)
      dnode(1,nj) = dprop_cell(1,nj)
      dnode(ni,1) = dprop_cell(ni,1)
      dnode(ni,nj) = dprop_cell(ni,nj)

!     Update the solution by adding the changes at the nodes "dnode" to the flow
!     property "prop"
!     INSERT
      prop = prop + dnode

      end subroutine sum_fluxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine smooth_array(av,prop)

!     This subroutine smooths "prop" to stabilise the calculation, the basic 
!     solver uses second order smoothing, many improvements are possible.

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg
      integer :: ni, nj

!     Get the block size and store locally for convenience
      ni = size(prop,1); nj = size(prop,2)

!     Calculate the average values at the nodes in the interior region of the
!     mesh, use the four neighbouring nodes in the plus and minus i and 
!     j-directions.
!     INSERT
      prop_avg(2:ni-1,2:nj-1)=(prop(1:ni-2,2:nj-1) &
                              + prop(3:ni,2:nj-1) &
                              + prop(2:ni-1,1:nj-2) &
                              + prop(2:ni-1,3:nj))*0.25

!     Edge values are also averaged in both the i and j-directions. Parallel to
!     the boundary the averaging is centred, the averages of two nodes are taken
!     either side of the current point. Perpendicular to the boundary the
!     algorithm is one-sided, the value at the current point is extrapolated
!     from the values at two nodes away from the boundary point.
!     INSERT
      prop_avg(:,1) = (prop(1:ni-2,1) &
                        + prop(3:ni,1) &
                        + 2*prop(2:ni-1,2) &
                        - prop(2:ni-1,3))/4
      prop_avg(:,nj) = (prop(1:ni-2,nj) &
                        + prop(3:ni,nj) &
                        + 2*prop(2:ni-1,nj-1) &
                        - prop(2:ni-1,nj-2))/4
      prop_avg(1,:) = (prop(1,1:nj-2) &
                        + prop(1,3:nj) &
                        + 2*prop(2,2:nj-1) &
                        - prop(3,2:nj-1))/4
      prop_avg(ni,:) = (prop(ni,1:nj-2) &
                        + prop(ni,3:nj) &
                        + 2*prop(ni-1,2:nj-1) &
                        - prop(ni-2,2:nj-1))/4

!     The corner values are not currently smoothed
      prop_avg([1,ni],[1,nj]) = prop([1,ni],[1,nj])
      ! prop_avg(1,1) = (prop(2,1) + prop(1,2))/2
      ! prop_avg(1,nj) = (prop(2,nj) + prop(1,nj-1))/2
      ! prop_avg(ni,1) = (prop(ni-1,1) + prop(ni,2))/2
      ! prop_avg(ni,nj) = (prop(ni-1,nj) + prop(ni,nj-1))/2

!     Now apply the artificial viscosity by smoothing "prop" towards "prop_avg",
!     take (1-sfac) * the calculated value of the property + sfac * the average 
!     of the surrounding values. 
!     INSERT
      prop = av%sfac*prop_avg + (1-av%sfac)*prop

      end subroutine smooth_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module stencils


