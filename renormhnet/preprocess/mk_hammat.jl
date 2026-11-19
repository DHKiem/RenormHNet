using DataFrames
using CSV
using ArgParse
#using LinearAlgebra
#using DelimitedFiles


#outputname = "Spin_Matrix_symmed.csv"

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--atom", "-a"
            help = "a list of the atom"
            nargs = '*'
            arg_type = Int
        "--radius", "-r"
            help = "a cutoff radius"
            arg_type = Float64        
        "--outputname", "-o"
            help = ""
            arg_type = String
            default = "./Spin_Matrix.csv"
    end
    return parse_args(s)
end

#=
function get_lattice(;latticefile = "lat.dat")
    lat = readdlm(latticefile)
    return lat
end

function get_positions(;positionfile = "site_positions.dat")
    position = readdlm(positionfile)
    return position
end
=#

parsed_args = parse_commandline()
outputname = parsed_args["outputname"]
atoms = parsed_args["atom"]
#radius = parsed_args["radius"]
#latticevec = get_lattice()
#atompos = get_positions()

rows = DataFrame(Rx=Int[], Ry=Int[], Rz=Int[],
                 Atom1=Int[], Atom2=Int[], J=Float64[])

for rx in -1:1, ry in -1:1, rz in -1:1
    for a1 in atoms, a2 in atoms
        if rx == 0 && ry == 0 && rz == 0 && a1 == a2
            continue
        end
        #=
        R = rx*latticevec[:, 1] .+ ry*latticevec[:, 2] .+ rz*latticevec[:, 3]
        distance = norm(atompos[:,a2] .- atompos[:,a1] .+ R)

        if distance > radius
            continue
        end
        =#
        #println((rx, ry, rz, a1, a2, 0.0), "  distance: ", distance)
        #println("Rx, Ry, Rz, Atom1, Atom2, J", rx, " ", ry, " ", rz , " ", a1, " ", a2, " " , 0)
        push!(rows, (rx, ry, rz, a1, a2, 0.0))
    end
end
CSV.write(outputname, rows)