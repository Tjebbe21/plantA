using Pkg
Pkg.activate(".")

using Images, ImageView
using Revise

includet("src/K3M.jl")             # K3M!()
includet("src/branchingpoints.jl") # angles_branchingpoint()

##

# load image
img = load("img/plant.jpg")
img_crop = img[1000:end-360,900:end-1100] # crop
# img_r = imresize(img_crop, ratio=0.5) # resize
# imgg1 = Gray.(img_r) # convert to grayscale
imgg1 = Gray.(img_crop) # convert to grayscale
imgg1 = imgg1 .< 0.70; # threshold




nr_iters = 50
bimg = copy(imgg1);

@time results,borders = K3M!(bimg,nr_iters);

save("output/k3m.gif",results)
save("output/k3m_b.gif",borders)

skelet = results[:,:,end];
imshow(skelet)

##

# Zwart wit skelet, maar met rgb type
rgb_skelet = RGB.(skelet,skelet,skelet)
# indices van branchingpoints
bp = idx_of_bp(skelet)

##

@time branches = decompose_skelet(skelet, bp);

##

col = distinguishable_colors(length(branches), [RGB(1,1,1),RGB(.1,.1,.1)])
for (i, edge) ∈ enumerate(values(branches))
    c = col[i]
    rgb_skelet[edge] .= c
end
imshow(rgb_skelet)

##
# 1 kind of in the middle
idx = bp[3]

# eerste stap in hoeken uitrekenen (en rgb_skelet krijgt gekleurde branches
angles_branchingpoint!(rgb_skelet,idx,10)


# Kleur branchpoints (en de px er om heen)
colored_bp = color_branching_points(skelet,RGB(1,0,0),RGB(1,0,0))
save("output/colored_bp.jpg",colored_bp)


# Poging tot bp labellen, maar kan het niet opslaan
guidict = imshow(colored_bp)
h_ofset = 0
for i in 1:length(bp)
    y,x = Tuple(bp[i])
    annotate!(guidict, AnnotationText(x+h_ofset, y-5, string(i), color=RGB(1,0,0), fontsize=5));
end
#canvas = guidict["gui"]["canvas"]

##

imshow(guidict)

##
# op github staat iets als dit, maar werkt niet..
using Cairo
function write_to_png(guidict, filename)
    canvas = guidict["gui"]["canvas"]
    ctx = getgc(canvas)
    Cairo.write_to_png(ctx.surface, filename)
end
write_to_png(guidict,"annotated")
