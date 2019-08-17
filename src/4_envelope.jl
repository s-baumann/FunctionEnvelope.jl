struct Envelope
    dd_interval_::DataFrame
    dd_evals_::DataFrame
    dd_function_::DataFrame
    x_name_::Symbol
    output_for_envelope_::Symbol
    splines_::Dict{Symbol,Dict{Symbol,Schumaker}}
    chance_of_function_::Multinomial
    funcs_::Array{Symbol,1}
    function Envelope(dd_interval::DataFrame, dd_eval::DataFrame, dd_function::DataFrame, x_name::Symbol, output_for_envelope::Symbol)
         func_names = dd_function[!,:name]
         val_names  = names(dd_eval)[(names(dd_eval) .!= :name) .& (names(dd_eval) .!= x_name)]
         splines = Dict{Symbol,Dict{Symbol,Schumaker}}()
         for name in val_names
             name_d = Dict{Symbol,Schumaker}()
             for func in func_names
                 data = dd_eval[dd_eval[!,:name] .== func,:]
                 if sum(ismissing.(data[!,name])) > 0 continue end
                 schum = Schumaker(data[!,x_name], data[!,name])
                 name_d[func] = schum
             end
             splines[name] = name_d
         end
         mass_by_func = by(dd_interval, :func, mass = :length => sum)
         total_mass   = sum(mass_by_func[!,:mass])
         mass_by_func[!,:mass] .= mass_by_func[!,:mass] ./ total_mass
         chance_of_function_ = Multinomial(1, mass_by_func[!,:mass])
         return new(dd_interval, dd_eval, dd_function, x_name, output_for_envelope, splines, chance_of_function_, mass_by_func[!,:func])
    end

    function Envelope(dd_function::DataFrame, grid::Union{Array,StepRangeLen}, output_for_envelope::Symbol,  x_name::Symbol = :x;
                      do_all_at_ends::Bool = true, bracketing_parameter::Real = 0.8, max_interval_width::Real = 0.01)
        dd_eval = crawler(dd_function, grid, x_name, output_for_envelope; do_all_at_ends = do_all_at_ends)
        dd_interval, dd_eval = get_upper_envelope(dd_eval, dd_function, x_name, output_for_envelope, bracketing_parameter, max_interval_width; grid = sort(grid))
        return Envelope(dd_interval, dd_eval, dd_function, x_name, output_for_envelope)
    end
end

import SchumakerSpline.evaluate, SchumakerSpline.evaluate_integral

"""
    evaluate_integral(env::Envelope, lhs::Real, rhs::Real, var::Symbol = env.output_for_envelope_)
Evaluates the integral of an envelope between two spots.
"""
function evaluate_integral(env::Envelope, lhs::Real, rhs::Real, variable::Symbol = env.output_for_envelope_)
    first_func_index = max(1, searchsortedlast(env.dd_interval_[!,:interval_start], lhs))
    last_func_index  = searchsortedlast(env.dd_interval_[!,:interval_start], rhs)
    spl_start_name   = env.dd_interval_[first_func_index, :func]
    integral         = 0.0
    if first_func_index == last_func_index
        #  If all is in the same interval.
        integral  += evaluate_integral(env.splines_[variable][spl_start_name], lhs, rhs)
        return integral
    else
        # Start and end components
        integral  += evaluate_integral(env.splines_[variable][spl_start_name], lhs, env.dd_interval_[first_func_index,:interval_end])
        spl_end_name = env.dd_interval_[last_func_index, :func]
        integral  += evaluate_integral(env.splines_[variable][spl_end_name], env.dd_interval_[last_func_index,:interval_start], rhs)
    end
    if last_func_index - first_func_index == 1
        # If start and end components are all they are.
        return integral
    end
    other_intervals = (first_func_index+1):1:(last_func_index-1)
    for i in other_intervals
        spl_name = env.dd_function_[i, :name]
        integral  += evaluate_integral(env.splines_[variable][spl_name], env.dd_interval_[i,:interval_start], env.dd_interval_[i,:interval_end])
    end
    return integral
end

function top_function_at_point(env::Envelope, x::Real)
    return env.dd_interval_[searchsortedlast(env.dd_interval_[!,:interval_start], x),:func]
end

"""
    evaluate(env::Envelope, x::Real, output_for_envelope::Symbol = env.output_for_envelope_; derivative::Integer = 0)
Evaluates an upper envelope at a particular coordinate. It can also estimate other inputs on the same intervals and derivatives.
"""
function evaluate(env::Envelope, x::Real, output_for_envelope::Symbol = env.output_for_envelope_; derivative::Integer = 0)
    top_spline_name = top_function_at_point(env, x)
    top_spline = env.splines_[output_for_envelope][top_spline_name]
    return evaluate(top_spline, x, derivative)
end
function (env::Envelope)(x::Real, derivative::Integer = 0)
    return evaluate(env, x, env.output_for_envelope_; derivative = derivative)
end





function sample_envelope_given_func(env::Envelope, func_name::Symbol, var::Symbol, twister::MersenneTwister; derivative::Integer = 0)
    eval_quantile = rand(twister,1)[1]
    x, y = sample_envelope_given_func(env, func_name, var, eval_quantile; derivative = 0)
    return x, y
end

function sample_envelope_given_func(env::Envelope, func_name::Symbol, var::Symbol, quantile::Real; derivative::Integer = 0)
    func_dd = env.dd_interval_[env.dd_interval_[:func] .== func_name,:]
    func_mass = sum(func_dd[:interval_end] .- func_dd[:interval_start])
    distance_to_point = func_mass * quantile
    final_x = NaN
    for i in 1:size(func_dd)[1]
        interval_width = func_dd[i,:interval_end] .- func_dd[i,:interval_start]
        if distance_to_point < interval_width
            final_x = func_dd[i,:interval_start] + distance_to_point
            break
        end
        distance_to_point = distance_to_point - interval_width
    end
    final_y = evaluate(env.splines_[var][func_name], final_x, derivative)
    return final_x, final_y
end

function sample_envelope(env::Envelope, var::Symbol, twister::MersenneTwister; derivative::Integer = 0)
    func_name = env.funcs_[findall(rand(twister, env.chance_of_function_) .== 1)[1]]
    x, y = sample_envelope_given_func(env, func_name, var, twister; derivative = 0)
    return func_name, x, y
end
