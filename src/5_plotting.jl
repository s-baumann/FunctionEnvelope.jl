import Plots.plot

function plot(env::Envelope, variable::Symbol = env.output_for_envelope_; points::Bool = true, legend::Symbol = :best)
    func_names = env.dd_function_[!,:name]
    intervals = size(env.dd_interval_)[1]
    plt = plot(; xlabel = string(env.x_name_), ylabel = string(variable))
    for i in 1:intervals
        int_func = env.dd_interval_[i, :func]
        int_ends = (env.dd_interval_[i, :interval_start] , env.dd_interval_[i, :interval_end])
        for j in 1:length(func_names)
            func = func_names[j]
            style = func == int_func ? :solid : :dash
            plot_options = (label = string(func), color = j, linestyle = style)
            plt = plot(env.splines_[variable][func], int_ends;derivs = false, plot_options = plot_options, plt = plt)
        end
    end
    if points
        plt = plot!(plt, env.dd_evals_[!,env.x_name_], env.dd_evals_[!,variable], seriestype=:scatter; color = :black, label = "points", legend=:none)
    end
    plt = plot!(plt; legend=legend)
    return plt
end
