using OctreeBH
using StaticArrays
using .Threads #for parallel treewalk
using PyPlot
plt.style.use("dark_background")

const N = 3 #spatial dimension
const T = Float64

const Npart = 2000 #number of particles
const ANGLE = 0.6 #node opening criteria
const FAC_GRAV = 0.05 #timestep control
const boxsizes = @SVector(ones(N)) * Inf  #open B.C.
const hsml = zeros(Npart) #dummy

const R0 = 1.0 #scale radius
const Rmax = 5.0 #cutoff radius
const softening = 5*R0 / sqrt(Npart) #gravitational softening length
const M_tot = 1.0 #total mass
const m_gas = M_tot / Npart
const mass = ones(Npart) .* m_gas #particle mass
const G_const = 1.0 #in code units

const tend = 10. #end time
const dt_dump = 0.02 #time between snapshots

using Random

function set_spherical_cloud(Npart, R0, Rmax)
    Random.seed!(1114);
    X = Vector{SVector{3,T}}(undef,Npart)
    V = [zero(SVector{N,T}) for _ in 1:Npart]
    for i in 1:Npart
        x,y,z = 0.,0.,0.
        r = 0.
        while(true)
            m = rand();
            r = R0 * sqrt(m) / (1. - sqrt(m)) #importance sampling
            if(r<Rmax) break end
        end
        phi = ( 2 * rand() - 1 ) * pi;
        cos_theta = 2 * rand() - 1;

        x = r * sqrt( 1. - cos_theta^2 ) * cos(phi);
        y = r * sqrt( 1. - cos_theta^2 ) * sin(phi);
        z = r * cos_theta;

        X[i] = SVector{3,T}(x,y,z)
    end
    return X,V
end


function calc_grav_acc(X)
    part = [Data{N,T}(SVector(X[i]), i, hsml[i], mass[i]) for i in eachindex(X)]

    Xmin = @SVector [minimum(getindex.(X,i)) for i in 1:N]
    Xmax = @SVector [maximum(getindex.(X,i)) for i in 1:N]

    topnode_length = (Xmax - Xmin) * 1.01  #actual length of tree
    center = 0.5 * (Xmax + Xmin)

    #build the tree
    tree = buildtree(part, center, topnode_length);

    acc = Vector{SVector{N,T}}(undef,Npart)

    @threads for i in eachindex(X)
        ga = GravTreeGather{N,T}()
        gravity_treewalk!(ga,X[i],tree,ANGLE,softening,boxsizes)
        acc[i] = G_const .* ga.acc
    end
    return acc
end

function plot_particle_configuration(X,count_dump)
    clf()
    fig, ax = PyPlot.subplots(1, 1, figsize=(7,6))
    ax.plot(getindex.(X,1), getindex.(X,2), ".", ms=1, c="white")
    ax.axis(R0.*[-1,1,-1,1])
    if count_dump < 10
        num = "000"*string(count_dump)
    elseif count_dump < 100
        num = "00"*string(count_dump)
    elseif count_dump < 1000
        num = "0"*string(count_dump)
    else
        num = string(count_dump)
    end
    savefig("plot_"*num*".png")
end

norm(a::SVector{N,T}) where{N,T} = sqrt(sum(a.^2))

function run()

    #randomly distributing particles
    #X = [@SVector rand(N) for _ in 1:Npart]
    #V = [zero(SVector{N,T}) for _ in 1:Npart]
    X,V = set_spherical_cloud(Npart, R0, Rmax)

    acc = calc_grav_acc(X)

    count_dump::Int64 = 0
    t = 0.
    t_output = 0.
    println("Start the calculation...")
    @time while t < tend
        if t >= t_output
            @show t
            plot_particle_configuration(X,count_dump)
            t_output += dt_dump
            count_dump += 1
        end

        dt = minimum( @. FAC_GRAV * sqrt(softening / norm(acc) ) )
        dt = tend - t < dt ? tend - t : dt

        #use the kick-drift-kick integration (i.e. midpoint method)
        V = V + 0.5 * acc * dt #1st kick
        X = X + V * dt #drift
        acc = calc_grav_acc(X) #calculate force using X^(n+1)
        V = V + 0.5 * acc * dt #2nd kick
        t += dt
    end
    return X,V
end

X,V = run();
0
