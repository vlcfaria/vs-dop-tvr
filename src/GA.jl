module GA

import Random
include("TimeFunctions.jl")
include("AcceleratedDubins.jl")
include("Helper.jl")

struct GAParams{T1,T2,T3,T4}
    popsize::Int
    cxpb::Float64
    mutpb::Float64
    ngen::Int

    selection::T1
    crossover::T2
    mutation!::T3

    op_params::T4
    perturb_chance::Float64 #TODO add chance to perturb waypoints -> random waypoint or fastest speed or highest score
end

mutable struct Individual{T}
    genome::T
    fitness::Float64
    valid::Bool
end

function evolution(params::GAParams, fitness_fun, initial_genomes)
    #Evaluate initial individuals
    population = map(genome -> Individual(genome, fitness_fun(params.op_params, genome), true), initial_genomes)
    sort!(population, by=x->x.fitness, rev=true)

    #TODO insert elitism, maybe also dont select from offspring and population?
    for gen in 1:params.ngen #Run generations
        #Generate our offspring
        offspring = params.selection(population, params.popsize)
        #Vary individuals
        offspring = varAnd(population, params, params.cxpb, params.mutpb)

        #Re-evaluate invalid fitnesses
        for i in eachindex(offspring)
            if !offspring[i].valid
                offspring[i].fitness = fitness_fun(params.op_params, offspring[i].genome)
                offspring[i].valid = true
            end
        end

        #Reconstruct the population from elite individuals + offspring
        #TODO

    end

    return population
end

#Calc fitness
function ind_fitness(op_params, seq)
    return Helper.calculate_seq_results(op_params,seq)[1]
end

#Apply both crossover and mutation on individuals
function varAnd(population, params::GAParams, cxpb, mutpb)
    offspring = [deepcopy(ind) for ind in population]
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
            params.mutation!(offspring[i].genome, params)
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
        aspirants = [rand(inds) for _ in 1:tournsize]
        # Choose the one with the max value
        push!(chosen, argmax(x->x.fitness, aspirants))
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
    method = rand(rng, [insertion, swap2])

    return method(recombinant, params)
end

#Takes a random point and inserts elsewhere
function insertion(recombinant::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)
    from, to = rand(rng, 2:l, 2) #Position 1 cant be swapped
    val = recombinant[from]
    deleteat!(recombinant, from)
    return insert!(recombinant, to, val)
end

#Swaps 2 points
function swap2(recombinant::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)
    p1,p2 = rand(rng, 2:l, 2) #Position 1 cant be swapped
    recombinant[p1], recombinant[p2] = recombinant[p2], recombinant[p1]
    return recombinant
end

#Changes waypoints
function waypoint_perturb(recombinant::T, params, rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)

    h,s = params.op_params.graph.num_headings, params.op_params.graph.num_speeds
    for i in eachindex(recombinant)
        if rand() < params.perturb_chance
            recombinant[idx] = (recombinant[idx][1], rand(1:s), rand(1:h))
        end
    end
end

end