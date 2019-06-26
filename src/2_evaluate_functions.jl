function evaluate_function(func::Union{NamedArray,Function},x_val::Real,func_name::Symbol,x_name::Symbol)
    new_dd = DataFrame()
    new_dd[x_name] = x_val
    feval = isa(func, NamedArray) ? func : func(x_val)
    if isa(feval, NamedArray)
        for name in names(feval)[1]
            new_dd[Symbol(name)] = feval[name]
        end
    else
        error("FunctionEnvelope doesn't know how to interpret what is being returned by function.")
    end
    new_dd[:name] = func_name
    return new_dd
end

function evaluate_function_and_add_to_dd(dd_eval::DataFrame, func::Union{NamedArray,Function},x_val::Real,func_name::Symbol,x_name::Symbol, output_for_envelope::Symbol) where T<:Real
    new_dd = evaluate_function(func,x_val,func_name,x_name)
    full_dd = vcat(dd_eval,new_dd, cols=:union)
    result = new_dd[1,output_for_envelope]
    return result, sort(full_dd, x_name)
end

function get_function_value(dd_eval::DataFrame, dd_function::DataFrame, x_val::Real, func_name::Symbol, x_name::Symbol, output_for_envelope::Symbol) where T<:Real
    if size(dd_eval)[1] > 0
        trying_point_in_cached_data = dd_eval[(dd_eval[:name] .== func_name) .& (dd_eval[x_name] .== x_val),:]
        if size(trying_point_in_cached_data)[1] > 1
            error("Duplication of same evaluation in dd_eval. Logical problem in code probably")
        elseif size(trying_point_in_cached_data)[1] == 1
            return trying_point_in_cached_data[1, output_for_envelope], dd_eval
        end
    end
    # In this case we don't have a saved result and so must evaluate.
    func = dd_function[dd_function[:name] .== func_name,:func][1]
    return evaluate_function_and_add_to_dd(dd_eval, func, x_val,func_name,x_name, output_for_envelope)
end
