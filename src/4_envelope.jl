struct Envelope
    dd_interval_::DataFrame
    dd_evals_::DataFrame
    dd_function_::DataFrame
    x_name::Symbol
    splines_::Dict{Symbol,Dict{Symbol,Schumaker}}
    function Envelope(dd_interval::DataFrame, dd_eval::DataFrame, dd_function::DataFrame, x_name::Symbol)
         func_names = dd_function[:name]
         val_names  = names(dd_eval)[(names(dd_eval) .!= :name) .& (names(dd_eval) .!= x_name)]
         splines = Dict{Symbol,Dict{Symbol,Schumaker}}()
         for name in val_names
             name_d = Dict{Symbol,Schumaker}()
             for func in func_names
                 data = dd_eval[dd_eval[:name] .== func,:]
                 schum = Schumaker(data[x_name], data[name])
                 name_d[func] = schum
             end
             splines[name] = name_d
         end
         return new(dd_interval, dd_eval, dd_function, x_name, splines)
    end
end
