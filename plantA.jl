using Images
using Pkg
Pkg.activate(".")

using Revise
using Random, ImageView


includet("K3M.jl")


# load image
img = load("plant.jpg")
img_crop = img[1000:end-360,900:end-1100] # crop
# img_r = imresize(img_crop, ratio=0.5) # resize
# imgg1 = Gray.(img_r) # convert to grayscale
imgg1 = Gray.(img_crop) # convert to grayscale
imgg1 = imgg1 .< 0.70; # threshold

# img_A = load("letterA.png")
# imgg1 = Gray.(img_A) # convert to grayscale
# imgg1 = imgg1 .< 0.45; # threshold

nr_iters = 50
bimg = copy(imgg1);
results,borders = K3M!(bimg,nr_iters);
save("k3m.gif",results)
save("k3m_b.gif",borders)

skelet = results[:,:,end];


skelet_green = RGB.(0,skelet,0)
save("skelet_green.jpg",skelet_green)

