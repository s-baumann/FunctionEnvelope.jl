function initial_evaluations(dd_function::DataFrame, limits::Tuple{T,T}, step_right_to_left::Bool, x_name::Symbol, output_for_envelope::Symbol; do_all_at_ends::Bool = false) where T<:Real
    from_side = limits[1+convert(Int,step_right_to_left)]
    to_side = limits[2-convert(Int,step_right_to_left)]
    dd_eval = DataFrame()
    for i in 1:size(dd_function)[1]
        func_name = dd_function[i,:name]
        type = dd_function[i,:type]
        # Everything gets done on the from side
        # On the to side we need decreasings if we step right to left and increasings if we step left to right.
        output, dd_eval = get_function_value(dd_eval, dd_function, from_side, func_name, x_name, output_for_envelope)
        #output, dd = evaluate_function_and_add_to_dd(dd, func,from_side,name,x_name, output_for_envelope)
        if ((type == :Decreasing) & step_right_to_left) | ((type == :Increasing) & (!step_right_to_left)) | do_all_at_ends
            output, dd_eval = get_function_value(dd_eval, dd_function, to_side, func_name, x_name, output_for_envelope)
            #output, dd = evaluate_function_and_add_to_dd(dd, func,to_side,name,x_name, output_for_envelope)
        end
    end
    return dd_eval
end

function maximum_possible_value(target::Real, dd_function::DataFrame, dd_eval::DataFrame,
                                output_for_envelope::Symbol, x_name::Symbol)
    dic = Dict{Symbol,typeof(dd_eval[output_for_envelope][1])}()
    for i in 1:size(dd_function)[1]
        name = dd_function[i,:name]
        type = dd_function[i,:type]
        reduced_evals = dd_eval[dd_eval[:name] .== name,:]
        if type == :Constant
            dic[name] = dd_function[i,:func][String(output_for_envelope)]
        elseif type == :Increasing
            index     = searchsortedlast(reduced_evals[x_name], target)+1
            if index > length(reduced_evals[x_name])
                dic[name] = Inf
            else
                dic[name] = reduced_evals[index,output_for_envelope]
            end
            dic[name] = reduced_evals[index,output_for_envelope]
        elseif type == :Decreasing
            index     = searchsortedfirst(reduced_evals[x_name], target)-1
            if index < 1
                dic[name] = Inf
            else
                dic[name] = reduced_evals[index,output_for_envelope]
            end
        else
            dic[name] = Inf
        end
    end
    return dic
end

function get_func(dd_function::DataFrame, func_name::Symbol)
    return dd_function[dd_function[:name] .== func_name,:func][1]
end

function functions_potentially_above_mark(marc::Real, max_vals::Dict{Symbol,T}) where T<:Real
    result = Array{Symbol,1}()
    for k in keys(max_vals)
        if max_vals[k] > marc
            push!(result, k)
        end
    end
    return result
end

function top_function_at_point(dd_eval::DataFrame, x_name::Symbol, output_for_envelope::Symbol, point::Real, tol = 1e-14)
    dd_at_point = dd_eval[abs.(dd_eval[x_name] .- point) .< tol,:]
    len = size(dd_at_point)[1]
    if len == 0 error("No functions evaluated at this point") end
    return sort(dd_at_point, output_for_envelope)[len, :name]
end

"""
This crawls from one side to the other to get the top function value at each point in a grid.
"""
function crawler(dd_function::DataFrame, limits::Tuple{T,T}, stepp::R, x_name::Symbol, output_for_envelope::Symbol;
                 step_right_to_left::Bool = true, do_all_at_ends::Bool = false) where T<:Real where R<:Real
    grid = limits[1]:step:limits[2]
    if step_right_to_left grid = reverse(grid) end
    return crawler(dd_function, grid, x_name, output_for_envelope; do_all_at_ends = do_all_at_ends)
end
function crawler(dd_function::DataFrame, grid::Union{Array,StepRangeLen}, x_name::Symbol, output_for_envelope::Symbol; do_all_at_ends::Bool = false)
    limits = (min(grid[1], grid[length(grid)]), max(grid[1], grid[length(grid)]))
    step_right_to_left = issorted(reverse(grid))
    if !step_right_to_left & !issorted(grid) error("The input grid must be strictly increasing or strictly decreasing") end
    dd_eval = initial_evaluations(dd_function, limits, step_right_to_left, x_name, output_for_envelope; do_all_at_ends = do_all_at_ends)
    top_function_at_last_point = top_function_at_point(dd_eval, x_name, output_for_envelope, grid[1])
    for i in 2:length(grid)
        target = grid[i]
        max_vals =  maximum_possible_value(target, dd_function, dd_eval, output_for_envelope, x_name)
        result, dd_eval = get_function_value(dd_eval, dd_function, target, top_function_at_last_point, x_name, output_for_envelope)
        #result, dd_eval = evaluate_function_and_add_to_dd(dd_eval, get_func(dd_function, top_function_at_last_point), target, top_function_at_last_point, x_name, output_for_envelope)
        other_functions_to_evaluate = functions_potentially_above_mark(result, max_vals)
        if length(other_functions_to_evaluate) > 0
            for f in other_functions_to_evaluate
                _, dd_eval = get_function_value(dd_eval, dd_function, target, f, x_name, output_for_envelope)
                #_, dd_eval = evaluate_function_and_add_to_dd(dd_eval, get_func(dd_function, f), target, f, x_name, output_for_envelope)
            end
        end
        top_function = top_function_at_point(dd_eval, x_name, output_for_envelope, grid[i])
        top_function_at_last_point = top_function
    end
    return dd_eval
end

# Intercept finders
"""
get_intercept_with_constant(dd_eval::DataFrame, func_name::Symbol, type::Symbol, constant::Real, dd_function::DataFrame, interval::Tuple{R,R}, x_name::Symbol, output_for_envelope::Symbol,
                                     bracketing_parameter::Real, max_interval_width::Real; recursion_count = 1, recursion_limit = 100) where R <:Real
Get an intercept between a function and a constant.
"""
function get_intercept_with_constant(dd_eval::DataFrame, func_name::Symbol, type::Symbol, constant::Real, dd_function::DataFrame, interval::Tuple{R,R}, x_name::Symbol, output_for_envelope::Symbol,
                                     bracketing_parameter::Real, max_interval_width::Real; recursion_count = 1, recursion_limit = 100) where R <:Real
    if recursion_count > recursion_limit
        error("We have gone over the recursion limit in the function get_intercept_with_constant. This can happen if bracketing_parameter is too close to 1.")
    end
    @assert (type == :Increasing) | (type == :Decreasing) "No constant or unspecified functions can be used in the get intercept with constant function."
    # Making the first function
    dd_1 = dd_eval[dd_eval[:,:name] .== func_name , :]
    s1 =  Schumaker(dd_1[:,x_name], dd_1[:,output_for_envelope])
    # Finding roots
    roots = find_roots(s1; root_value = constant)[:roots]
    root_in_interval = roots[(roots .>= interval[1]-10*eps()) .& (roots .<= interval[2]+10*eps())]
    num_roots = length(root_in_interval)
    if num_roots == 0
        error("No roots were found. Potentially more evaluations of this function are needed.")
    elseif num_roots > 1
        error("This package does not support multiple roots in an interval.")
    end
    if abs(interval[2] - interval[1]) < max_interval_width
        return root_in_interval[1], dd_eval
    end
    # Finding point to try.
    point_to_try =  (1-bracketing_parameter) * (sum(interval)/2) + bracketing_parameter * (root_in_interval[1])
    # Now trying this point and saving it in dd_eval.
    result, dd_eval = get_function_value(dd_eval, dd_function, point_to_try, func_name, x_name, output_for_envelope)
    # Now using the result to make a smaller interval.
    if ((result > constant) & (type == :Increasing)) | ((result <= constant) & (type == :Decreasing))  # Intercept is left
        new_interval = Tuple{R,R}([interval[1], point_to_try])
        return get_intercept_with_constant(dd_eval, func_name, type, constant, dd_function, new_interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width; recursion_count = recursion_count + 1, recursion_limit = recursion_limit)
    else
        new_interval = Tuple{R,R}([point_to_try, interval[2]])
        return get_intercept_with_constant(dd_eval, func_name, type, constant, dd_function, new_interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width; recursion_count = recursion_count + 1, recursion_limit = recursion_limit)
    end
end

"""
get_intercept(dd_eval::DataFrame, increasing_function::Symbol, decreasing_function::Symbol, dd_function::DataFrame, interval::Tuple{R,R}, x_name::Symbol,
              output_for_envelope::Symbol, bracketing_parameter::Real, max_interval_width::Real; recursion_count = 1, recursion_limit = 100) where R <:Real
Get an intercept between two functions.
"""
function get_intercept(dd_eval::DataFrame, increasing_function::Symbol, decreasing_function::Symbol, dd_function::DataFrame, interval::Tuple{R,R}, x_name::Symbol,
                       output_for_envelope::Symbol, bracketing_parameter::Real, max_interval_width::Real; recursion_count = 1, recursion_limit = 100) where R <:Real
    if recursion_count > recursion_limit
        error("We have gone over the recursion limit in the function get_intercept. This can happen if bracketing_parameter is too close to 1.")
    end
    # Making the first function
    dd_1 = dd_eval[dd_eval[:,:name] .== increasing_function , :]
    s1 =  Schumaker(dd_1[:,x_name], dd_1[:,output_for_envelope])
    # Making the second function
    dd_2 = dd_eval[dd_eval[:,:name] .== decreasing_function , :]
    s2 =  Schumaker(dd_2[:,x_name], dd_2[:,output_for_envelope])
    # Finding roots
    roots = get_intersection_points(s1, s2)
    root_in_interval = roots[(roots .>= interval[1]-10*eps()) .& (roots .<= interval[2]+10*eps())]
    point_to_try = (sum(interval)/2)
    if length(root_in_interval) == 1 # 2 or 0 roots can happen sometimes if the spline is a bit inaccurate.
        point_to_try =  (1-bracketing_parameter) * (sum(interval)/2) + bracketing_parameter * (root_in_interval[1])
        if abs(interval[2] - interval[1]) < max_interval_width
            return root_in_interval[1], dd_eval
        end
    else
        if abs(interval[2] - interval[1]) < max_interval_width
            return (sum(interval)/2), dd_eval
        end
    end
    # Evaluation of function.
    result_inc, dd_eval = get_function_value(dd_eval, dd_function, point_to_try, increasing_function, x_name, output_for_envelope)
    #result1, dd_eval = evaluate_function_and_add_to_dd(dd_eval, get_func(dd_function, increasing_function), point_to_try, increasing_function, x_name, output_for_envelope)
    result_dec, dd_eval = get_function_value(dd_eval, dd_function, point_to_try, decreasing_function, x_name, output_for_envelope)
    #result2, dd_eval = evaluate_function_and_add_to_dd(dd_eval, get_func(dd_function, decreasing_function), point_to_try, decreasing_function, x_name, output_for_envelope)
    if result_inc > result_dec # We are to the right of the intercept as the increasing function is higher.
        new_interval = Tuple{R,R}([interval[1], point_to_try])
        return get_intercept(dd_eval, increasing_function, decreasing_function, dd_function, new_interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width; recursion_count = recursion_count + 1, recursion_limit = recursion_limit)
    else
        new_interval = Tuple{R,R}([point_to_try, interval[2]])
        return get_intercept(dd_eval, increasing_function, decreasing_function, dd_function, new_interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width; recursion_count = recursion_count + 1, recursion_limit = recursion_limit)
    end
end

"""
    get_upper_envelope(dd_eval::DataFrame, dd_function::DataFrame, x_name::Symbol, output_for_envelope::Symbol, bracketing_parameter::Real, max_interval_width::Real;
                             grid::Array{R,1} = sort(unique(dd_eval[x_name]))) where R<:Real
"""
function get_upper_envelope(dd_eval::DataFrame, dd_function::DataFrame, x_name::Symbol, output_for_envelope::Symbol, bracketing_parameter::Real, max_interval_width::Real;
                                      grid::Union{StepRangeLen,Array{R,1}} = sort(unique(dd_eval[x_name]))) where R<:Real
    interval_dd = DataFrame(interval_start = grid[1], interval_end = NaN, func = top_function_at_point(dd_eval, x_name, output_for_envelope, grid[1]))
    i = 1 # This is the row index for interval_dd
    for j in 2:length(grid)
        # We need to ensure that the previous function has been evaluated at this point. Otherwise we might get a dodgy "change"
        _, dd_eval = get_function_value(dd_eval, dd_function, grid[j], interval_dd[i, :func], x_name, output_for_envelope)
        top_func = top_function_at_point(dd_eval, x_name, output_for_envelope, grid[j])
        if top_func != interval_dd[i, :func]
            old_top_func = interval_dd[i, :func]
            old_func_type = dd_function[dd_function[:,:name] .== old_top_func,:type][1]
            top_func_type = dd_function[dd_function[:,:name] .== top_func    ,:type][1]
            intercept = NaN
            interval = Tuple([grid[j-1], grid[j]])
            if (old_func_type == :Constant) && (top_func_type == :Constant)
                error("This does not make sense. How did they cross?")
            elseif (old_func_type == :Constant)
                old_func_const = dd_function[dd_function[:,:name] .== old_top_func,:func][1][String(output_for_envelope)]
                intercept, dd_eval = get_intercept_with_constant(dd_eval, top_func, top_func_type, old_func_const, dd_function, interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width)
            elseif (top_func_type == :Constant)
                top_func_const = dd_function[dd_function[:,:name] .== top_func,:func][1][String(output_for_envelope)]
                intercept, dd_eval = get_intercept_with_constant(dd_eval, old_top_func, old_func_type, top_func_const, dd_function, interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width)
            else
                intercept, dd_eval = get_intercept(dd_eval, top_func, old_top_func, dd_function, interval, x_name, output_for_envelope, bracketing_parameter, max_interval_width)
            end
            interval_dd[i, :interval_end] = intercept
            i+=1
            new_interval = DataFrame(interval_start = intercept, interval_end = NaN, func = top_func)
            interval_dd = vcat(interval_dd, new_interval)
        end
    end
    interval_dd[i,:interval_end] = grid[length(grid)]
    interval_dd[:length] = interval_dd[:interval_end] .- interval_dd[:interval_start]
    return interval_dd, dd_eval
end
