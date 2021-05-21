using StaticArrays
using DSP 
using Printf
using Combinatorics
using AngleBetweenVectors # https://github.com/JeffreySarnoff/AngleBetweenVectors.jl

"""
calculate weight to determine if point is a branchingpoint

bit_img: Binary image with 0 as background and 1 the figure 
(x,y):   Index of pixel
"""
function branch_weight(bit_img,x,y)
    BSA = SA[0 1 0; 1 0 1; 0 1 0]
    weight = 0
    for i ∈ 1:3
        for j ∈ 1:3
            weight += BSA[i,j] * bit_img[x+i-2,y+j-2]
        end
    end
    return weight
end

"""
detect branchingpoints.
Set elements in zero array branchingpoints to 1 if that point is a branching point in skelet

branchingpoints: array filled with zeros
skelet: 1 px wide binary image 
"""
function branchingpoints!(branchingpoints,skelet)
    n,m = size(skelet)
    for x in 1:n, y in 1:m
        if skelet[x,y] == 1
            if branch_weight(skelet,x,y) ≥ 3
                branchingpoints[x,y] = 1
            end
        end
    end
end

function branchingpointsv2!(branchingpoints, skelet)
    for I ∈ findall(skelet .== 1)
        if branch_weight(skelet, I.I...) ≥ 3
            branchingpoints[I] = 1
        end
    end
end

"""
returns list of CartesianIndices of branchingpoints

skelet: 1 px wide binary image
"""
function idx_of_bp(skelet)
    branchingpoints = zero(skelet);
    @time branchingpointsv2!(branchingpoints,skelet)

    idx_of_bp = [id for id in CartesianIndices(branchingpoints) if branchingpoints[id] == 1]
    return idx_of_bp
end



"""
Paints nr_steps of branch in direction in color
idx is starting branch point.

returns 
- (direction, (0,0)) if there is no branch in direction
- (direction, CartesianIndex) of branch after nr_steps or end of branch or next branching point
"""
function few_steps_in_branch(rgb_skelet,idx,nr_steps,direction,color)
    # initialize previous idx to branch point to avoid waking backwards
    previous_idx = idx

    # Take first step in specified direction
    if direction == "right"
        idx = idx + CartesianIndex(0,1)
    elseif direction == "left"
        idx = idx + CartesianIndex(0,-1)
    elseif direction == "top"
        idx = idx + CartesianIndex(-1,0)
    elseif direction == "bottom"
        idx = idx + CartesianIndex(1,0)
    end

    # Return (direction, (0,0)) if there is no branch in direction (for type stabilty)
    if skelet[idx] != 1
        println("no branch in direction: " * direction)
        return (direction, CartesianIndex(0,0))
    end

    # Paint first step
    rgb_skelet[idx] = color
    
    # Take and paint nr_steps steps on the branch
    for _ in 1:nr_steps
        # define possible walking directions
        below = idx + CartesianIndex(1,0)
        top   = idx + CartesianIndex(-1,0)
        right = idx + CartesianIndex(0,1)
        left  = idx + CartesianIndex(0,-1)

        # select next direction
        directions = [dir for dir in [right,left,top,below] if (skelet[dir] == 1) && (dir != previous_idx)]
       
        # Set previous_idx to avoid waking backwards
        previous_idx = idx

        # Directions should have 1 length. Otherwise we encountered a new branchpoint on the branch or the end
        if length(directions) == 1
            idx = directions[1]
            rgb_skelet[idx] = color # paint branch
        elseif length(directions) > 1
            println("multiple directions")
            return idx
        elseif length(directions) == 0
            println("end of the branch")
            return (direction,idx)
        end
    end

    # return last coordinate
    return (direction,idx)
end

"""
Prints angles between branches at branchingpoint idx.
Also colors the branches

idx is CartesianIndex of branchingpoint
nr_steps is length of branch to calculate angles
"""
function angles_branchingpoint!(rgb_skelet,idx,nr_steps)
    # Find end coordinates in the branches
    coord_right  = few_steps_in_branch(rgb_skelet,idx,nr_steps,"right",  colorant"#ffbc0a")
    coord_left   = few_steps_in_branch(rgb_skelet,idx,nr_steps,"left",   colorant"#c200fb")
    coord_top    = few_steps_in_branch(rgb_skelet,idx,nr_steps,"top",    colorant"#ec0868")
    coord_bottom = few_steps_in_branch(rgb_skelet,idx,nr_steps,"bottom", colorant"#fc2f00")

    # Translate coordinates s.t. idx is origin (if branch is there)
    # The != (0,0) is a choice in few_steps_in_branch.
    vecs = [(d, Float64.(Tuple(c  - idx))) for (d,c) in [coord_right,coord_left,coord_top,coord_bottom] if c != CartesianIndex(0,0)]

    # Calculate all possible pairs of branches
    pers = combinations(vecs,2)
    
    # For every pair of branches calculate and print the angle between them
    for p in pers
        p1, p2 = p
        d1, c1 = p1
        d2, c2 = p2
        phi = angle(c1,c2) * 180/ π
        @printf("angle between %6s and %6s = %3.1f deg\n",d1,d2,phi)
    end
end

"""
Decomposes the tree

Inputs:
    skelet: the skelet matrix
    bp: the branch points

Returns:
    dictionary:
        keys: (origin, target) origin/target are brnaching points/tip of the branch
        values: the path connecting the origin and the target as a list
"""
function decompose_skelet(skelet, bp)

    book_skelet = copy(skelet)
    book_skelet[bp].= 0

    # define neighbors
    neighbors = [CartesianIndex(x,y) for (x,y) in [ (0,1),(1,0),(0,-1),(-1,0) ] ]

    branches = Dict()
    for origin in bp
        # See the potential neighbors or a bp For each we make a branch
        for n in neighbors
            if book_skelet[origin + n] == 1
                current_node = origin+n
                edge = [current_node]
                book_skelet[current_node] = 0
                search = true
                while search
                    search = false
                    for nn in neighbors
                        if book_skelet[current_node+nn] == 1
                            search = true
                            current_node += nn
                            append!(edge, [current_node])
                            book_skelet[current_node] = 0
                        end
                    end
                end
                target = edge[end]
                for nn in neighbors
                    if (edge[end] + nn) in bp
                        target += nn
                    end
                end


                branches[(origin, target)] = edge    
            end                    
        end
    end
    return branches
end

#----------------------------------

"""
Paint branchingpoints and neighbors
returns rbg image

skelet: 1 px wide binary image

color_bp:     color of branchingpoints
color_bp_nbh: color of neighbor of branchingpoint
"""
function color_branching_points(skelet,color_bp,color_bp_nbh)
    n,m = size(skelet)
    skelet_colored = RGB.(skelet,skelet,skelet)
    for x in 1:n, y in 1:m
        if skelet[x,y] == 1
            if branch_weight(skelet,x,y) ≥ 3
                skelet_colored[x,y] = color_bp
                for i in -1:1, j in -1:1
                    skelet_colored[x+i,y+j] = color_bp_nbh
                end
            end
        end
    end
    return skelet_colored
end
