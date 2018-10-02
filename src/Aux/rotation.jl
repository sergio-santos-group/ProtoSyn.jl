#TO DO:
# 1) DOCUMENT THE FUNCTION
function rotation_matrix_from_axis_angle(axis::Vector{Float64}, angle::Float64)
    q0 = cos(0.5 * angle)
    q1, q2, q3 = sin(0.5 * angle) * axis ./ norm(axis)
    [1-2*q2*q2-2*q3*q3   2*q1*q2-2*q0*q3   2*q1*q3+2*q0*q2;
       2*q2*q1+2*q0*q3 1-2*q3*q3-2*q1*q1   2*q2*q3-2*q0*q1;
       2*q3*q1-2*q0*q2   2*q3*q2+2*q0*q1 1-2*q1*q1-2*q2*q2]
end