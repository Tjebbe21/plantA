using StaticArrays
using DSP 

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

function detect_branching_points(skelet)
    n,m = size(skelet)
    skelet_colored = RGB.(0,skelet,0)
    for x in 1:n
        for y in 1:m
            if skelet[x,y] == 1
                if branch_weight(skelet,x,y) ≥ 3
                    skelet_colored[x,y] = RGB(1,0,0)
                end
            end
        end
    end
    return skelet_colored
end

function detect_branching_pointsv2(skelet)
    skelet_colored = RGB.(0,skelet,0)
    for i ∈ findall(skelet .== 1)
        if branch_weight(skelet, i.I...) ≥ 3
            skelet_colored[i] = RGB(1,0,0)
        end
    end
    return skelet_colored
end

function detect_branching_points_conv(skelet)
    kernel = [ 0 1 0; 1 10 1; 0 1 0]
    conv_skelet = @view conv(kernel, skelet)[2:end-1, 2:end-1]
    skelet_colored = RGB.(0,skelet, 0)

    for i ∈ findall(conv_skelet .≥ 13)
        skelet_colored[i] = RGB(1,0,0)
    end
    return skelet_colored
end