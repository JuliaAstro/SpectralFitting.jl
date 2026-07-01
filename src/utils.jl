function _convolve_implementation!(
    output::AbstractVector{T},
    vec_A::AbstractVector{T},
    kernel::AbstractVector{T},
) where {T<:Number}
    # Based on https://discourse.julialang.org/t/97658/15
    J = length(vec_A)
    K = length(kernel)
    @assert length(output) == J + K - 1 "Output is $(length(output)); should be $(J + K - 1)"

    # do the kernel's side first
    for i = 1:(K-1)
        total = zero(T)
        for k = 1:K
            ib = (i >= k)
            oa = ib ? vec_A[i-k+1] : zero(T)
            total += kernel[k] * oa
        end
        output[i] = total
    end
    # now the middle
    for i = K:(J-1)
        total = zero(T)
        for k = 1:K
            oa = vec_A[i-k+1]
            total += kernel[k] * oa
        end
        output[i] = total
    end
    # and finally the end
    for i = J:(J+K-1)
        total = zero(T)
        for k = 1:K
            ib = (i < J + k)
            oa = ib ? vec_A[i-k+1] : zero(T)
            total += kernel[k] * oa
        end
        output[i] = total
    end
    output
end

function convolve!(output, A, kernel)
    if length(kernel) <= length(A)
        _convolve_implementation!(output, A, kernel)
    else
        _convolve_implementation!(output, kernel, A)
    end
end
function convolve(A, kernel)
    output = zeros(eltype(A), length(A) + length(kernel) - 1)
    convolve!(output, A, A_domain, kernel, kernel_domain)
end

"""
Assumes A is binned in X1 and kernel is binned in X2. Output will also be binned on X1
"""
function _convolve_irregular_grid!(output, A, X1, kernel, X2)
    @assert length(X1) == length(A) + 1
    @assert length(X2) == length(kernel) + 1

    function _kernel_func(x)
        i1 = searchsortedfirst(X2, x)
        if i1 >= length(X2) || i1 <= 1
            return zero(x)
        end
        w = (x - X2[i1-1]) / (X2[i1] - X2[i1-1])
        w * kernel[i1] + (1 - w) * kernel[i1-1]
    end

    fill!(output, 0)
    for i in eachindex(output)
        for j in eachindex(output)
            x = X1[i] / X1[j]
            k = A[j] * _kernel_func(x)
            if k > 0
                output[i] += k
            end
        end
    end
end

function convolve!(output, A, A_domain, kernel, kernel_domain)
    _convolve_irregular_grid!(output, A, A_domain, kernel, kernel_domain)
end
function convolve(A, A_domain, kernel, kernel_domain)
    output = zeros(eltype(A), length(A) + length(kernel) - 1)
    convolve!(output, A, A_domain, kernel, kernel_domain)
end

# poisson statistics and related utility functionso

"""
    count_error(k, σ)

Gives the error on `k` mean photon counts to a significance of `σ` standard deviations
about the mean.

Derived from likelihood of binomial distributions being the beta function.
"""
function count_error(k, σ)
    p = Distributions.cdf(Distributions.Normal(), σ)
    kₑ = gamma_inc_inv(k + 1, p, 1 - p)
    abs(k - kₑ)
end
