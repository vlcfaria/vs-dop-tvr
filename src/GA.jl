module GA

import Random
include("TimeFunctions.jl")
include("AcceleratedDubins.jl")
include("Helper.jl")
include("Vns.jl")

struct GAParams{T1,T2,T3,T4}
    popsize::Int
    cxpb::Float64
    mutpb::Float64
    ngen::Int
    elite_percentage::Float64

    selection::T1
    crossover::T2
    mutation::T3
    op_params::T4

    perturb_chance::Float64 #TODO add chance to perturb waypoints -> random waypoint or fastest speed or highest score
end

mutable struct Individual{T}
    genome::T
    fitness::Float64
    valid::Bool
end

function evolution(params::GAParams, fitness_fun, initial_genomes)
    #Calculate number of elite individuals
    elite_size = ceil(Int64, params.elite_percentage * params.popsize)
    offspring_size = params.popsize - elite_size

    #Evaluate initial individuals
    population = map(genome -> Individual(genome, fitness_fun(params.op_params, genome), true), initial_genomes)
    local_search(population, params)
    sort!(population, by=x->x.fitness, rev=true)

    for gen in 1:params.ngen #Run generations
        #Generate our offspring
        offspring = params.selection(population, offspring_size)
        #Vary individuals
        offspring = varAnd(offspring, params, params.cxpb, params.mutpb)

        #Re-evaluate invalid fitnesses
        for i in eachindex(offspring)
            if !offspring[i].valid
                offspring[i].fitness = fitness_fun(params.op_params, offspring[i].genome)
                offspring[i].valid = true
            end
        end

        #Reconstruct the population from elite individuals + offspring
        #Population is sorted by fitness, just pick the first as elite
        population[elite_size+1:end] = offspring

        local_search(population, params)

        sort!(population, by=x->x.fitness, rev=true) #Sort again
        println((gen, population[1].fitness))
    end

    return population
end

function local_search(population, params, l_max=3)
    #Apply one full iteration of VNS local searches
    len = length(population[1].genome)
    tmax = params.op_params.tmax

    @Threads.threads for ind_idx in eachindex(population)
        ind = population[ind_idx]

        best_sequence = ind.genome
        best_score, best_time, _ = Helper.calculate_seq_results(params.op_params, best_sequence)
        l = 1
        while l <= l_max
            local_sequence = deepcopy(best_sequence)
            local_score, local_time, local_limit_idx = Helper.calculate_seq_results(params.op_params, local_sequence)
            #Search
            for _ in 1:len^2
                search_seq, _ = Vns.search(deepcopy(local_sequence), params.op_params.graph, l, local_limit_idx)
                search_score, search_time, search_limit_idx = Helper.calculate_seq_results(params.op_params, search_seq)

                # Check if searched solution is better OR equal -> Higher score within tmax OR same score, lower time
                if (search_score > local_score && search_time <= tmax) || (search_score == local_score && search_time <= local_time)
                    local_sequence = search_seq
                    local_time = search_time
                    local_score = search_score
                    local_limit_idx = search_limit_idx
                end
            end

            #Higher score found through local search OR equal score with lower time
            if (local_time <= tmax && local_score > best_score) || (local_score == best_score && local_time < best_time)
                best_time = local_time
                best_sequence = local_sequence
                best_score = local_score
                l = 1
            else
                l += 1
            end
        end

        #Override individual genome
        ind.genome = best_sequence
        ind.fitness = best_score
    end
end

#Calc fitness
function ind_fitness(op_params, seq)
    return Helper.calculate_seq_results(op_params,seq)[1]
end

#Apply both crossover and mutation on individuals
function varAnd(offspring, params::GAParams, cxpb, mutpb)
    #Assumes offspring was already copied
    l = length(offspring)
    for i in 2:2:l
        if rand() < cxpb
            #TODO this can be optimized a little bit, since we are copying again in the function
            offspring[i-1].genome, offspring[i].genome = params.crossover(offspring[i-1].genome, offspring[i].genome, params)
            offspring[i-1].valid = offspring[i].valid = false
        end
    end

    for i in eachindex(offspring)
        if rand() < mutpb
            #Not all mutations are in-place, so assign
            offspring[i].genome = params.mutation(offspring[i].genome, params)
            offspring[i].valid = false
        end
    end

    return offspring
end

#Tournament selection
function selection_tournament(inds, k, tournsize=2)
    chosen = []
    for _ in 1:k
        # Choose tournsize number aspirants
        aspirants = rand(inds,tournsize)
        # Choose copy of the one with the max value
        push!(chosen, deepcopy(argmax(x->x.fitness, aspirants))) 
    end

    return chosen
end

#Order crossover, adapted from Evolutionary.jl
function order_crossover(v1::T, v2::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    s = length(v1)
    from, to = rand(rng, 1:s, 2) #You actually can switch index 1 here
    from, to = from > to ? (to, from)  : (from, to)
    c1 = fill((0,0,0), s)
    c2 = fill((0,0,0), s)
    # Swap
    c1[from:to] = v2[from:to]
    c2[from:to] = v1[from:to]
    # Also mantain a set of the first value in the tuple
    in_c1, in_c2 = Set([x[1] for x in c1]), Set([x[1] for x in c2])
    
    # Fill in from parents -> Get first element
    k = 1
    j = 1 
    #We also dont need to check if to+1 is in bounds, because if from == 1 and to+1 > len(v1), for loop below doesnt run

    for i in vcat(1:from-1, to+1:s)
        while in(v1[k][1],in_c1) #while value is already in new child
            k = k+1 # No need to wrap around here
        end
        c1[i] = v1[k] #Set value
        push!(in_c1, v1[k][1])

        while in(v2[j][1],in_c2) #Do the same for c2
            j = j+1
        end
        c2[i] = v2[j]
        push!(in_c2, v2[j][1])
    end
    return c1, c2
end

function mutation(recombinant::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    inversion_mutation = (x,y) -> displacement_mutation(x,y,true)
    method = rand(rng, [displacement_mutation, inversion_mutation, waypoint_perturb])

    return method(recombinant, params)
end

#Also acts as the inversion mutation operator
function displacement_mutation(recombinant::T, params, rev=false, rng=Random.default_rng()) where {T <: AbstractVector}
    #Select a subtour at random, and insert it into a random place
    start, stop = rand(rng, 2:length(recombinant), 2)
    if stop < start
        start,stop = stop,start
    end

    #Extract subtour
    subtour = recombinant[start:stop]
    if rev reverse!(subtour) end

    h,s = params.op_params.graph.num_headings, params.op_params.graph.num_speeds

    for i in eachindex(subtour)
        if rand() < params.perturb_chance
            subtour[i] = (subtour[i][1], rand(1:s), rand(1:h))
        end
    end

    removed = vcat(recombinant[1:start-1], recombinant[stop+1:end])

    insertion = rand(rng, 1:length(removed))
    return vcat(removed[1:insertion], subtour, removed[insertion+1:end])
end

#Takes a random point and inserts elsewhere
function insertion_mutation(recombinant::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)
    from, to = rand(rng, 2:l, 2) #Position 1 cant be swapped
    val = recombinant[from]
    deleteat!(recombinant, from)
    
    if (rand() < params.perturb_chance)
        h,s = params.op_params.graph.num_headings, params.op_params.graph.num_speeds
        val = (val[1], rand(1:s), rand(1:h))
    end
    insert!(recombinant, to, val)
    return recombinant
end

#Changes waypoints of all positions according to chance
function waypoint_perturb(recombinant::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)

    h,s = params.op_params.graph.num_headings, params.op_params.graph.num_speeds
    for idx in eachindex(recombinant)
        if rand() < params.perturb_chance
            recombinant[idx] = (recombinant[idx][1], rand(1:s), rand(1:h))
        end
    end

    return recombinant
end

end