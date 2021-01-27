using Statistics
using LinearAlgebra

function align!(mobile::Pose, target::Pose, selection::AbstractSelection)
    
    # 1) Apply selection
    mobile_coords = copy(mobile.state.x.coords[:, selection(mobile).content])
    target_coords = copy(target.state.x.coords[:, selection(target).content])

    # 2) Get centroid and rotation matrix

    # 2.1) Center each centroid on [0, 0, 0]
    c_mobile   = mean(mobile_coords, dims = 2)
    c_target   = mean(target_coords, dims = 2)
    mobile_coords   .-= c_mobile
    target_coords   .-= c_target

    # 2.2) Obtain correlation matrix
    cm         = mobile_coords * target_coords'
    u, d, vt   = svd(cm)
    rot_matrix = vt * u'

    # 3) Apply transformation (rotation+ translation)
    centered_mobile = mobile.state.x.coords .- c_mobile
    mobile.state.x[:,1:end] = (centered_mobile' * rot_matrix')' .+ c_target

    return mobile
end

# function get_centroid_rot_matrix(mobile::Matrix{Float64}, target::Matrix{Float64})

#     # 1) Center each centroid on [0, 0, 0]
#     c_mobile   = mean(mobile, dims = 2)
#     c_target   = mean(target, dims = 2)
#     mobile   .-= c_mobile
#     target   .-= c_target

#     # 2) Obtain correlation matrix
#     cm         = mobile * target'
#     u, d, vt   = svd(cm)
#     rot_matrix = vt * u'

#     return c_mobile, c_target, rot_matrix
# end