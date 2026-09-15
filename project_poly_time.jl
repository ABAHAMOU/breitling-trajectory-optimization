using JuMP, Cbc, Printf

function nextCity(x::Matrix{Int64}, i::Int64, n::Int64)
    for j in 1:n
        if x[i,j] == 1
            return j
        end
    end
    return -1
end

function findPath(x::Matrix{Int64}, start::Int64, n::Int64)
    path = [start]
    current = start
    visited = falses(n)
    visited[start] = true
    while true
        next = nextCity(x, current, n)
        if next == -1 || visited[next]
            break
        end
        visited[next] = true
        push!(path, next)
        current = next
    end
    return path
end

function readInstance(filename::String)
    open(filename, "r") do f
        lines = [strip(line) for line in readlines(f) if strip(line) != ""]
        idx = 1
        n       = parse(Int64, lines[idx]); idx += 1
        s       = parse(Int64, lines[idx]); idx += 1
        t       = parse(Int64, lines[idx]); idx += 1
        Vmin    = parse(Int64, lines[idx]); idx += 1
        _m      = parse(Int64, lines[idx]); idx += 1
        regions = parse.(Int64, split(lines[idx])); idx += 1
        Rmax    = parse(Int64, lines[idx]);          idx += 1
        coords  = Matrix{Int64}(undef, n, 2)
        for i in 1:n
            parts = parse.(Int64, split(lines[idx])); idx += 1
            coords[i, 1] = parts[1]
            coords[i, 2] = parts[2]
        end
        d = Matrix{Int64}(undef, n, n)
        for i in 1:n, j in 1:n
            dx = coords[i,1] - coords[j,1]
            dy = coords[i,2] - coords[j,2]
            d[i,j] = round(Int64, sqrt(dx^2 + dy^2))
        end
        return n, s, t, Vmin, regions, Rmax, coords, d
    end
end

function aerodromeProblem(d::Matrix{Int64}, s::Int64, t::Int64,
                          Vmin::Int64, regions::Vector{Int64}, Rmax::Int64,
                          time_limit::Float64 = 300.0)
    n = size(d, 1)
    m = maximum(regions)

    model = Model(Cbc.Optimizer)
    set_optimizer_attribute(model, "seconds", time_limit)

    # ------------------------------------------------------------------
    # Modèle ILP — formulation MTZ (nombre polynomial de contraintes)
    #
    # Variables :
    #   x[i,j] ∈ {0,1}  : 1 si l'arc i→j est emprunté
    #   y[i]   ∈ {0,1}  : 1 si l'aérodrome i est visité
    #   u[i]   ∈ [1,n]  : ordre de visite de i (variable de position MTZ)
    #
    # Objectif : minimiser la distance totale
    #   min  Σ_{i,j} d[i,j] * x[i,j]
    # ------------------------------------------------------------------
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, y[1:n], Bin)
    @variable(model, 1 <= u[1:n] <= n)

    # Arcs invalides : boucle ou distance > Rmax
    for i in 1:n, j in 1:n
        if i == j || d[i,j] > Rmax
            @constraint(model, x[i,j] == 0)
        end
    end

    @objective(model, Min, sum(d[i,j] * x[i,j] for i in 1:n, j in 1:n))

    # Un seul départ depuis s, une seule arrivée en t, pas de retour sur s ou départ depuis t
    @constraint(model, sum(x[s,j] for j in 1:n) == 1)
    @constraint(model, sum(x[i,s] for i in 1:n) == 0)
    @constraint(model, sum(x[i,t] for i in 1:n) == 1)
    @constraint(model, sum(x[t,j] for j in 1:n) == 0)

    # Conservation du flot : si i est visité, il a exactement un arc entrant et un sortant
    for i in 1:n
        if i != s && i != t
            @constraint(model, sum(x[i,j] for j in 1:n) == y[i])
            @constraint(model, sum(x[j,i] for j in 1:n) == y[i])
        end
    end

    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)

    # Nombre minimum d'aérodromes visités
    @constraint(model, sum(y[i] for i in 1:n) >= Vmin)

    # Chaque région visitée au moins une fois
    for r in 1:m
        @constraint(model, sum(y[i] for i in 1:n if regions[i] == r) >= 1)
    end

    # ------------------------------------------------------------------
    # Contraintes MTZ — élimination des sous-tours en O(n²)
    #
    # Pour chaque arc x[i,j] actif entre deux nœuds visités,
    # on impose u[i] < u[j] (j est visité après i dans le chemin).
    # Sous forme linéaire :
    #   u[i] - u[j] + n*x[i,j]  <=  (3n-1) - n*y[i] - n*y[j]
    #
    # Si y[i]=y[j]=1 (les deux nœuds visités) : contrainte MTZ classique
    # Si y[i]=0 ou y[j]=0 (nœud non visité)  : contrainte inactive (membre droit grand)
    # ------------------------------------------------------------------
    @constraint(model, u[s] == 1)
    for i in 1:n, j in 1:n
        if i != j && j != s
            @constraint(model,
                u[i] - u[j] + n * x[i,j] <= (3*n - 1) - n*y[i] - n*y[j]
            )
        end
    end

    set_silent(model)
    optimize!(model)

    status = string(termination_status(model))

    if termination_status(model) ∉ [MOI.OPTIMAL, MOI.FEASIBLE_POINT]
        return nothing, Inf, status
    end

    xsol = round.(Int64, value.(x))
    ysol = round.(Int64, value.(y))
    dist = objective_value(model)

    path = findPath(xsol, s, n)
    println("  Chemin   : ", path)
    println("  Visités  : ", sum(ysol), " aérodromes")
    println("  Régions  : ", unique([regions[i] for i in path if regions[i] != 0]))

    return model, dist, status
end

function solveInstance(filename::String; time_limit::Float64 = 300.0)
    println("\n========================================")
    println("Instance : $filename")
    println("========================================")
    if !isfile(filename)
        println("  ERREUR: fichier introuvable.")
        return
    end
    n, s, t, Vmin, regions, Rmax, coords, d = readInstance(filename)
    println("  n=$n  s=$s  t=$t  Vmin=$Vmin  Rmax=$Rmax")
    t_start = time()
    model, dist, status = aerodromeProblem(d, s, t, Vmin, regions, Rmax, time_limit)
    t_elapsed = time() - t_start
    println("----------------------------------------")
    println("  Statut          : $status")
    println("  Distance min    : $(dist < Inf ? round(Int, dist) : "N/A")")
    @printf("  Temps CPU       : %.3f s\n", t_elapsed)
    println("========================================")
end

function main()
    instances = [
        "instance_6_1.txt",

    ]
    if length(ARGS) >= 1
        solveInstance(ARGS[1])
    else
        for f in instances
            solveInstance(f)
        end
    end
end

main()
