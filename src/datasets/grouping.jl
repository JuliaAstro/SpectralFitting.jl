struct GroupingIterator
    grouping::Vector{Int}
end

function Base.length(grouping::GroupingIterator)
    return count(==(1), grouping.grouping)
end

const GroupingIteratorItem = Union{Nothing, Tuple{Int, Int, Int}}
const GroupingIteratorState = Tuple{Int, Union{Nothing, Int}}
const GroupingIteratorItemAndState = Tuple{GroupingIteratorItem, GroupingIteratorState}

Base.eltype(::GroupingIterator) = Tuple{Int, Int, Int}

function Base.iterate(grouping::GroupingIterator)
    i1 = findfirst(==(1), grouping.grouping)
    if isnothing(i1)
        return nothing
    end
    return iterate(grouping, (0, i1))
end

function Base.iterate(grouping::GroupingIterator, state)
    i::Int, i1::Int = state
    if i > lastindex(grouping.grouping) || i1 > lastindex(grouping.grouping)
        return nothing
    end
    index = findnext(==(1), grouping.grouping, i1 + 1)
    i2::Int = if isnothing(index)
        lastindex(grouping.grouping)
    else
        index - 1
    end
    return (i + 1, i1, i2), (i + 1, i2 + 1)
end

function regroup!(vector::Vector, grouping::Vector{Int})
    itt = GroupingIterator(grouping)
    for grp in itt
        regroup_vector!(vector, grp)
    end
    resize!(vector, length(itt))
    return vector
end

regroup_vector!(vector::Vector, grouping::Tuple{Int, Int, Int}) =
    vector[grouping[1]] = @views sum(vector[grouping[2]:grouping[3]])
regroup_vector!(output::Vector, vector::Vector, grouping::Tuple{Int, Int, Int}) =
    output[grouping[1]] = @views sum(vector[grouping[2]:grouping[3]])

# function _regroup_warnings(data::SpectralDataset)
#     if !all(==(0), data.spectrum.quality[data.mask])
#         @warn "Dataset contains bad channels, and regrouping may be erroneous. Proceed only if you know what you're doing.\nBad channels masked out by using `mask_bad_channels!`."
#     end
# end

const BAD_QUALITY = 1
const GOOD_QUALITY = 0
function regroup_quality_vector!(quality::Vector{Int}, grouping::Tuple{Int, Int, Int})
    qs = @views quality[grouping[2]:grouping[3]]
    if all(==(GOOD_QUALITY), qs)
        quality[grouping[1]] = GOOD_QUALITY
    else
        quality[grouping[1]] = BAD_QUALITY
    end
    return
end

function warn_bad_channels(bad_channels)
    warns = String[]
    N = length(bad_channels)
    j = 1
    for (i, Δ) in enumerate(diff(bad_channels))
        if Δ != 1
            if i - j == 0
                push!(warns, "$i")
            else
                push!(warns, "$j-$i")
            end
            j = i
        end
    end
    if N - j == 0
        push!(warns, "$j")
    else
        push!(warns, "$j-$N")
    end
    chans = join(warns, ", ")
    @warn "Grouped channels $chans contain bad quality channels."
    return
end

export regroup!
