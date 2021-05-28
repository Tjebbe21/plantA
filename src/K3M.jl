using StaticArrays
using Images


"""
Neighborhood weight calculation

bit_img: Binary image with 0 as background and 1 the figure
N:       Neighbourhood bit values matrix
(x,y):   Index of pixel
"""
@inline function weight(bit_img,N,x,y)
    weight = 0
    for i ∈ 1:3, j ∈ 1:3
        weight += N[i,j] * bit_img[x+i-2,y+j-2]
    end
    return weight
end

"""
Make boundaries of bit_img background

bit_img:          Binary image with 0 as background and 1 the figure
boundary_width :  Number of pixels from boundray that should be background (0)
"""
function boundary_background(bit_img,boundary_width)
     for i in 1:boundary_width
        bit_img[i,:]       .= 0 # top
        bit_img[:,i]       .= 0 # left
        bit_img[end+1-i,:] .= 0 # bottom
        bit_img[:,end+1-i] .= 0 # right
    end
end

"""
phase 0 : Marking borders

border:  Bookeeping array of size bit_img, 1: px is border, 0 if not.
bit_img: Binary image with 0 as background and 1 the figure
N:       Neighbourhood bit values matrix
"""
function phase0!(border,bit_img,N,A0)
    # Check boundaries of bit_img are background
    # assert_boundary_background(bit_img,3)

    n, m = size(bit_img)
    for x ∈ 2:n-1, y ∈ 2:m-1
        border[x,y] = weight(bit_img,N,x,y) ∈ A0
    end 
end

"""
phase 1-5: Delete borders
"""
function phaseᵢ!(border,bit_img,N,A)
    n,m = size(border)
    for x in 2:n-1, y in 2:m-1
        if border[x,y] == 1
            if weight(bit_img,N,x,y) ∈ A
                bit_img[x,y] = 0
            end
        end
    end
end

"""
Unmark the remaining borders
"""
@inline function phase6!(border)
    border .= 0
end

"""
K3M algorithm
See Saeed et al. (2010) Appl. Math. Compt. Sci. DOI: 10.2478/v10006-010-0024-4

bit_img:   Binary image with 0 as background and 1 the figure 
nr_iters:  Number of iterations (to do: stop creterion)
"""
function K3M!(bit_img,nr_iters)
    # Staticarrays
    N    = SA[128 1 2; 64 0 4; 32 16 8];
    A0   = SA[3, 6, 7, 12, 14, 15, 24, 28, 30, 31, 48, 56, 60, 62, 63, 96, 112, 120, 124, 126, 127, 129, 131, 135,143, 159, 191, 192, 193, 195, 199, 207, 223, 224, 225, 227, 231, 239, 240, 241, 243, 247, 248, 249, 251, 252, 253, 254];
    A1   = SA[7, 14, 28, 56, 112, 131, 193, 224];
    A2   = SA[7, 14, 15, 28, 30, 56, 60, 112, 120, 131, 135, 193, 195, 224, 225, 240];
    A3   = SA[7, 14, 15, 28, 30, 31, 56, 60, 62, 112, 120, 124, 131, 135, 143, 193, 195, 199, 224, 225, 227, 240, 241, 248];
    A4   = SA[7, 14, 15, 28, 30, 31, 56, 60, 62, 63, 112, 120, 124, 126, 131, 135, 143, 159, 193, 195, 199, 207, 224, 225, 227, 231, 240, 241, 243, 248, 249, 252];
    A5   = SA[7, 14, 15, 28, 30, 31, 56, 60, 62, 63, 112, 120, 124, 126, 131, 135, 143, 159, 191, 193, 195, 199, 207, 224, 225, 227, 231, 239, 240, 241, 243, 248, 249, 251, 252, 254];
    A1px = SA[3, 6, 7, 12, 14, 15, 24, 28, 30, 31, 48, 56, 60, 62, 63, 96, 112, 120, 124, 126, 127, 129, 131, 135, 143, 159, 191, 192, 193, 195, 199, 207, 223, 224, 225, 227, 231, 239, 240, 241, 243, 247, 248, 249, 251, 252, 253, 254];
    
    n,m      = size(bit_img)       # img size
    border   = zeros(n,m)          # Bookeeping array
    bit_imgs = zeros(n,m,nr_iters+5) # Save output of every iteration
    borders  = zeros(n,m,nr_iters+5) # Save border of every iteration

    # Make boundary (of 1 px wide) part of the background
    boundary_background(bit_img, 1)
    
    bit_imgs[:,:,1:5] .= bit_img   # few frames of start for gif

    for i in 5:nr_iters+5
        bit_imgs[:,:,i] .= bit_img   # save bit_im
        
        phase0!(border,bit_img,N,A0) # marking borders
        borders[:,:,i] .= border     # save borders

        phaseᵢ!(border,bit_img,N,A1) # deleting pixels having 3 sticking neighbours
        phaseᵢ!(border,bit_img,N,A2) # deleting pixels having 3 or 4 sticking neighbours
        phaseᵢ!(border,bit_img,N,A3) # deleting pixels having 3, 4 or 5 sticking neighbours
        phaseᵢ!(border,bit_img,N,A4) # deleting pixels having 3, 4, 5 or 6 sticking neighbours
        phaseᵢ!(border,bit_img,N,A5) # deleting pixels having 3, 4, 5, 6 or 7 sticking neighbours
        phase6!(border)              # reset border selection
    end
    phaseᵢ!(border,bit_img,N,A1px)   # thinning to a one-pixel width skeleton
    
    return bit_imgs,borders
end

"""
input: binary image
output: border of image
"""
function border(bit_img)
    N  = SA[128 1 2; 64 0 4; 32 16 8];
    A0 = SA[3, 6, 7, 12, 14, 15, 24, 28, 30, 31, 48, 56, 60, 62, 63, 96, 112, 120, 124, 126, 127, 129, 131, 135,143, 159, 191, 192, 193, 195, 199, 207, 223, 224, 225, 227, 231, 239, 240, 241, 243, 247, 248, 249, 251, 252, 253, 254];
    
    n, m   = size(bit_img)       # img size
    border = zeros(n,m)          # Bookeeping array
    phase0!(border,bit_img,N,A0)  # marking border 
    return border
end



#-----------------------------------------
# Trial

"""
phase 1-5: Delete borders
without if statement in double for loop.

resulting K3M was a bit slower with this...
maybe better if indexes array is pre allocated?
"""
function phaseᵢ_pif!(border,bit_img,N,A)
    n,m = size(border)
    indexes = Tuple.([id for id in CartesianIndices(bit_img) if 
                (bit_img[id] == 1) &  # select indices of points that are in the skeleton
                (1 < id[1] < n) &     # exclude boundaries
                (1 < id[2] < m)
    ])

    for (x,y) in indexes
        if weight(bit_img,N,x,y) ∈ A
            bit_img[x,y] = 0
        end
    end
end
