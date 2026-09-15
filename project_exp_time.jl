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

# Détecte un sous-tour isolé dans la solution x.
# Retourne :valid si le chemin s->t est correct,
# :subtour si un cycle déconnecté existe,
# :incomplete si le chemin depuis s ne rejoint pas t.
function findSubtour(x::Matrix{Int64}, n::Int64, s::Int64, t::Int64)
    visited = falses(n)
    path_s = findPath(x, s, n)
    for node in path_s
        visited[node] = true
    end
    for i in 1:n
        if !visited[i] && sum(x[i,:]) > 0
            subtour = findPath(x, i, n)
            if length(subtour) >= 2
                return subtour, :subtour
            end
        end
    end
    if !isempty(path_s) && path_s[end] == t
        return path_s, :valid
    else
        return path_s, :incomplete
    end
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
    # Modèle ILP — élimination de sous-tours (nombre exponentiel de contraintes)
    #
    # Variables :
    #   x[i,j] ∈ {0,1}  : 1 si l'arc i→j est emprunté
    #   y[i]   ∈ {0,1}  : 1 si l'aérodrome i est visité
    #
    # Objectif : minimiser la distance totale
    #   min  Σ_{i,j} d[i,j] * x[i,j]
    # ------------------------------------------------------------------
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, y[1:n], Bin)

    # Arcs invalides : boucle ou distance > Rmax
    for i in 1:n, j in 1:n
        if i == j || d[i,j] > Rmax
            @constraint(model, x[i,j] == 0)
        end
    end

    @objective(model, Min, sum(d[i,j] * x[i,j] for i in 1:n, j in 1:n))

    # Un seul départ depuis s, une seule arrivée en t, pas de retour sur s ou départ depuis t
    @constraint(model, sum(x[s,j] for j in 1:n) == 1)
    @constraint(model, sum(x[i,t] for i in 1:n) == 1)
    @constraint(model, sum(x[i,s] for i in 1:n) == 0)
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

    set_silent(model)
    optimize!(model)

    if termination_status(model) ∉ [MOI.OPTIMAL, MOI.FEASIBLE_POINT]
        return nothing, Inf, string(termination_status(model)), 0
    end

    xsol = round.(Int64, value.(x))
    ysol = round.(Int64, value.(y))
    iteration = 0

    # ------------------------------------------------------------------
    # Élimination itérative de sous-tours
    #
    # Si la solution contient un sous-tour isolé S (cycle déconnecté du chemin s→t),
    # on ajoute la contrainte :
    #   Σ_{i ∈ S, j ∈ S} x[i,j]  <=  |S| - 1
    # qui interdit que tous les arcs internes à S soient actifs simultanément.
    # On résout à nouveau et on répète jusqu'à obtenir un chemin valide s→t.
    # ------------------------------------------------------------------
    while true
        tour, status_tour = findSubtour(xsol, n, s, t)

        if status_tour == :valid
            break

        elseif status_tour == :subtour
            iteration += 1
            println("  Itération $iteration: sous-tour éliminé (taille $(length(tour))): $tour")
            @constraint(model, sum(x[i,j] for i in tour, j in tour) <= length(tour) - 1)

        elseif status_tour == :incomplete
            # Le chemin depuis s ne rejoint pas t mais il n'y a pas de sous-tour isolé.
            # On interdit la solution courante entière (no-good cut).
            println("  Itération $(iteration+1): chemin incomplet, ajout no-good")
            active_arcs = [(i,j) for i in 1:n, j in 1:n if xsol[i,j] == 1]
            if isempty(active_arcs)
                break
            end
            iteration += 1
            @constraint(model, sum(x[i,j] for (i,j) in active_arcs) <= length(active_arcs) - 1)
        end

        optimize!(model)

        st = termination_status(model)
        if st ∉ [MOI.OPTIMAL, MOI.FEASIBLE_POINT]
            println("  Plus de solution (statut: $st).")
            return nothing, Inf, string(st), iteration
        end

        xsol = round.(Int64, value.(x))
        ysol = round.(Int64, value.(y))
    end

    dist   = objective_value(model)
    status = string(termination_status(model))
    path   = findPath(xsol, s, n)

    println("  Chemin   : ", path)
    println("  Visités  : ", sum(ysol), " aérodromes")
    println("  Régions  : ", unique([regions[i] for i in path if regions[i] != 0]))

    return model, dist, status, iteration
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
    model, dist, status, iters = aerodromeProblem(d, s, t, Vmin, regions, Rmax, time_limit)
    t_elapsed = time() - t_start
    println("----------------------------------------")
    println("  Statut          : $status")
    println("  Distance min    : $(dist < Inf ? round(Int, dist) : "N/A")")
    println("  Itérations      : $iters")
    @printf("  Temps CPU       : %.3f s\n", t_elapsed)
    println("========================================")
end

function main()
    #il faut insérer le fichier pour l'utiliser
    instances = [
        
        "custom_maps/instance_custom_1.txt",
        "custom_maps/instance_custom_2.txt",
        "custom_maps/instance_custom_3.txt",

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
