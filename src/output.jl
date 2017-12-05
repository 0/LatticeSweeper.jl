"""
    SweepOutput

Output generator for sweep information.
"""
abstract type SweepOutput end

"""
    init(::SweepOutput)

Initialize the `SweepOutput` before sweeping.
"""
init(::SweepOutput) = nothing

"""
    step(::SweepOutput, ::SweepState)

Update the `SweepOutput` with `SweepState` after a single step.
"""
step(::SweepOutput, ::SweepState) = nothing

"""
    sweep(::SweepOutput, ::SweepState)

Update the `SweepOutput` with `SweepState` after an entire sweep.
"""
sweep(::SweepOutput, ::SweepState) = nothing

"""
    done(::SweepOutput)

Finalize the `SweepOutput` before sweeping.
"""
done(::SweepOutput) = nothing


"""
    SweepOutputFile

`SweepOutput` that writes to a file.
"""
struct SweepOutputFile <: SweepOutput
    io::IOStream
end

"""
    SweepOutputFile(path::String)

Create a `SweepOutputFile` that will write to `path`.
"""
SweepOutputFile(path::String) = SweepOutputFile(open(path, "w"))

function init(so::SweepOutputFile)
    println(so.io, "# sweep_num energy dH2 SvN duration")
    flush(so.io)

    nothing
end

function sweep(so::SweepOutputFile, state::SweepState)
    items = Any[state.sweep_num, state.energy, state.dH2,
                S_vn(state.middle_eigvals), state.duration]
    println(so.io, join(items, " "))
    flush(so.io)

    nothing
end


"""
    SweepOutputDynamic

`SweepOutput` that shows sweep information in real time.
"""
mutable struct SweepOutputDynamic <: SweepOutput
    "Minimum time (in seconds) between outputs."
    time_step::Float64
    "Time of last output (in seconds)."
    time_prev::Float64
end

"""
    SweepOutputDynamic(time_step::Float64=0.1)

Create a `SweepOutputDynamic` that will output no more than once every
`time_step` seconds.
"""
SweepOutputDynamic(time_step::Float64=0.1) = SweepOutputDynamic(time_step, 0.0)

function init(::SweepOutputDynamic)
    println()

    for _ in 1:6
        println()
    end

    nothing
end

function make_string{L}(::SweepOutputDynamic, state::SweepState{L};
                        sweep=false)
    io = IOBuffer()

    for _ in 1:6
        # Go up, clear line.
        print(io, "\u1b[A\u1b[K")
    end

    # Don't exceed 80 columns.
    if L > 78
        scaling = Int(ceil((L-1)/77))
        L_eff = Int(ceil((L-1)/scaling))+1
        site_eff = div(state.site-1, scaling)+1
    else
        L_eff = L
        site_eff = state.site
    end

    print(io, "[")
    for j in 1:(site_eff-1)
        print_with_color(:blue, io, "-")
    end
    if sweep
        print_with_color(:blue, io, "--")
    else
        if state.dir == Right
            print_with_color(:red, io, ">>")
        elseif state.dir == Left
            print_with_color(:red, io, "<<")
        end
    end
    for j in (site_eff+2):L_eff
        print_with_color(:blue, io, "-")
    end
    print(io, "]\n")

    if sweep
        @printf(io, "done sweep %d\n", state.sweep_num)
    else
        @printf(io, "on sweep %d at site %d of %d\n",
                state.sweep_num, state.site, L)
    end
    if isnan(state.duration)
        println(io, "")
    else
        @printf(io, "previous sweep took %.3f s\n", state.duration)
    end
    if isnan(state.energy)
        println(io, "")
    else
        @printf(io, "energy is %.5g\n", state.energy)
    end
    if isempty(state.middle_eigvals)
        println(io, "")
    else
        @printf(io, "von Neumann entropy across middle bond is %.5g\n",
                S_vn(state.middle_eigvals))
    end
    if isnan(state.dH2)
        println(io, "")
    else
        @printf(io, "Hamiltonian fluctuation (dH2) is %.5g\n", state.dH2)
    end

    io |> take! |> String
end

function step{L}(so::SweepOutputDynamic, state::SweepState{L})
    time_now = time()
    time_now >= so.time_prev + so.time_step || return

    print(make_string(so, state; sweep=false))

    so.time_prev = time_now

    nothing
end

function sweep{L}(so::SweepOutputDynamic, state::SweepState{L})
    time_now = time()

    print(make_string(so, state; sweep=true))

    # Don't rate limit the end-of-sweep updates, but still update the time.
    so.time_prev = time_now

    nothing
end

function done(::SweepOutputDynamic)
    println()

    nothing
end
