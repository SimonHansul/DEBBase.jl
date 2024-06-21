
@recipe function f(p::Priors; q_limits = 1e-5)
    param_names = p.params
    num_params = length(param_names)
    layout --> (1,num_params)

    for (i,param) in enumerate(param_names)
        @series begin
            dist = p.priors[i]
            seriestype := :path
            subplot := i
            limits = quantile.(dist, (q_limits, 1-q_limits))
            xlabel --> param
            color --> :black
            if i == 1 
                ylabel --> "Density"
            end
            x = range(limits..., length = 100)
            y = pdf.(dist, x)
            x,y
        end
    end
end

# TODO: add lower diagonal
@recipe function f(res::SMCResult; q_limits = 1e-5)
    param_names = fieldnames(typeof(res.priors))
    num_params = length(param_names)
    layout := (num_params, num_params)

    idx = 0
    for (i,param1) in enumerate(param_names)
        for (j,param2) in enumerate(param_names)
            idx += 1
            if i == j
                @series begin # plot prior
                    dist = getfield(res.priors, param1)
                    seriestype := :path
                    subplot := idx
                    limits = quantile.(dist, (q_limits, 1-q_limits))
                    xlabel --> param1
                    color --> :black
                    linestyle --> :dash
                    x = range(limits..., length = 100)
                    y = pdf.(dist, x)
                    x,y
                end
                @series begin # plot marginal accepted kde
                    U = DEBABC.wkde(res.accepted, param1)
                    subplot := idx
                    color --> :gray
                    fill --> true
                    fillalpha --> .5
                    xmin = minimum(res.accepted[:,param1]) * 0.1
                    xmax = maximum(res.accepted[:,param1]) * 1.1
                    x = range(xmin, xmax, length = 100) #(minimum(accepted[:,param1] .* 0.01), maximum(accepted[:,param1]))
                    y = [pdf(U, xi) for xi in x]
                    x,y
                end
                @series begin # rugplot of accepted values
                    subplot := idx
                    seriestype := :scatter
                    markershape := :vline
                    color --> :gray
                    y = zeros(nrow(res.accepted))
                    x = res.accepted[:,param1]

                    x,y
                end
            elseif i>j # below diagonal : 2d KDE
                # FIXME: plotting of binkde causes cryptic error message
                subplot := idx
                seriestype := :contour
                U = DEBABC.wkde(res.accepted, param1, param2)
                #xlim = (minimum(accepted[:,param1] * 1), maximum(accepted[:,param1]) * 1.1)
                #ylim = (minimum(accepted[:,param2] * 1), maximum(accepted[:,param2]) * 1.1)
                #x = range(xlim..., length = 100)
                #y = range(ylim..., length = 100)
                #z = hcat([[pdf(U, xi,yi) for xi in y] for yi in y]...) |> Matrix
                (collect(U.x), collect(U.y), collect(U.density'))
            end
        end
    end
end
