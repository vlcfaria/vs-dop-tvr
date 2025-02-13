module GA

import Random

struct GAParams{T1,T2,T3}
    popsize::Int
    cxpb::Float64
    mutpb::Float64
    ngen::Int

    selection::T1
    crossover::T2
    mutation::T3
end

mutable struct Individual{T}
    genome::T
    fitness::Float64
    valid::Bool
end

function evolution(params::GAParams, fitness_fun, initial_genomes)
    #Evaluate initial individuals
    population = map(genome -> Individual(genome, fitness_fun(genome), true), initial_genomes)

    for gen in 1:ngen #Run generations
        #Generate our offspring
        offspring = varAnd(population, params, cxpb, mutpb)

        #Re-evaluate invalid fitnesses
        for i in eachindex(offspring)
            if !offspring[i].valid
                offspring[i].fitness = fitness_fun(offspring[i].genome)
                offspring[i].valid = true
            end
        end

        #Select next generation, select from population AND offspring!
        population = params.select(vcat(population, offspring), params.popsize)
    end

    return population
end

#Apply both crossover and mutation on individuals
function varAnd(population, params::GAParams, cxpb, mutpb)
    offspring = [deepcopy(ind) for ind in population]
    l = length(offspring)
    for i in 2:2:l
        if rand() < cxpb
            offspring[i-1].genome, offspring[i].genome = params.crossover(offspring[i-1].genome, offspring[i].genome)
            offspring[i-1].valid = offspring[i].valid = false
        end
    end

    for i in eachindex(offspring)
        if rand() < mutpb
            offspring[i].genome = params.mutation(offspring[i].genome)
            offspring[i].valid = false
        end
    end

    return offspring
end

#Order crossover, adapted from Evolutionary.jl
function order_crossover(v1::T, v2::T; rng=Random.default_rng()) where {T <: AbstractVector}
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
    k = from == 1 ? to+1 : 1 #child1 index before modifications
    j = from == 1 ? to+1 : 1 #child2 index before modifications

    for i in vcat(1:from-1, to+1:s)
        while in(v1[k][1],in_c1) #while value is already in new child
            k = k+1 > s ? 1 : k+1
        end
        c1[i] = v1[k] #Set value
        push!(in_c1, v1[k][1])

        while in(v2[j][1],in_c2) #Do the same for c2
            j = j+1 > s ? 1 : j+1
        end
        c2[i] = v2[j]
        push!(in_c2, v2[j][1])
    end
    return c1, c2
end

function mutation(recombinant::T; rng=Random.default_rng()) where {T <: AbstractVector}
    method = rand(rng, [insertion, swap2])

    println(method)
    return method(recombinant)
end

#Takes a random point and inserts elsewhere
function insertion(recombinant::T; rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)
    from, to = rand(rng, 2:l, 2) #Position 1 cant be swapped
    val = recombinant[from]
    deleteat!(recombinant, from)
    return insert!(recombinant, to, val)
end

#Swaps 2 points
function swap2(recombinant::T; rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)
    p1,p2 = rand(rng, 2:l, 2) #Position 1 cant be swapped
    evol.swap!(recombinant, p1, p2)
    return recombinant
end

#Changes waypoints
function waypoint_perturb(recombinant::T, rng=Random.default_rng()) where {T <: AbstractVector}
    l = length(recombinant)
    idx = rand(rng, 1:l)
    recombinant[idx] = (recombinant[idx][1], 1, 1) #TODO find a way to pass this as param
end

end