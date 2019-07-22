function make_function_dataframe(symbols::Array{Symbol,1}, function_array::Union{Array{Any,1},Array{Function,1},Array{NamedArray,1},Array{Union{NamedArray,Function},1}}, types::Array{Symbol,1}) where T<:Real
    if sum(isa.(function_array, Union{Real,Function,NamedArray})) < length(function_array) error("Function array can only contain reals, Functions or NamedArrays") end
    return DataFrame(name = symbols,
                   func = function_array,
                   type = types)
end
