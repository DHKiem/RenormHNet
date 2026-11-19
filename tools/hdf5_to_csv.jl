```
Author: Do Hoon Kiem
works for Jx
```
using HDF5
using DataFrames
using CSV
using ArgParse
using DelimitedFiles
using LinearAlgebra


@info "This code works for Jx only."
s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "Input file name.h5"
        arg_type = String
        default = "hamiltonians_pred.h5"         
    "--output", "-o"
        help = "Output file name.csv"
        arg_type = String
        default = "Hamiltonian.csv"    
    "--cutoff"
        help = "Cutoff distance for interactions (optional, default: no cutoff), then it reads lat.dat and site_positions.dat"
        arg_type = Float64
end

parsed_args = parse_args(s)
filename = parsed_args["input"]
outputname = parsed_args["output"]
cutoff = get(parsed_args, "cutoff", nothing)
println("Input file: $filename")
println("Output file: $outputname")

file_dir = dirname(abspath(filename))
println("file directory: $file_dir")
function get_lattice(;latticefile = "lat.dat")
    lat = readdlm(file_dir * "/" * latticefile)
    return lat
end
lat = get_lattice()

function get_positions(;positionfile = "site_positions.dat")
    position = readdlm(file_dir * "/" * positionfile)
        return position
end
position = get_positions()

if cutoff !== nothing
    println("Cutoff distance: $cutoff")
end



DF = []
h5open(filename, "r") do f
    for (key, obj) in pairs(f)
                
        clean_str = strip(key, ['[', ']'])
        key = parse.(Int, split(clean_str, ","))
        
        if obj isa HDF5.Dataset
            data = read(obj)
            # This part allocates one spin basis per one atom.
            J = data[1,1]
            #println(key," ", J)

            # Calculate the distance between the two atoms
            atom1_pos = position[:, key[4]]
            atom2_pos = position[:, key[5]]    
            position_vec = key[1] * lat[:,1] + key[2] * lat[:,2] + key[3] * lat[:,3] + (atom2_pos - atom1_pos)
            distance = norm(position_vec)
            if distance == 0.0 
                continue
            end
            if cutoff !== nothing
                if distance > cutoff 
                    continue
                end
                push!(DF, (Rx=key[1], Ry=key[2], Rz=key[3], Atom1=key[4], Atom2=key[5], J=J, Distance=distance))
            else
                push!(DF, (Rx=key[1], Ry=key[2], Rz=key[3], Atom1=key[4], Atom2=key[5], J=J, Distance=distance))
            end
            
        end
    end
end
DF = DataFrame(DF)
println(first(DF,5))
CSV.write(outputname, DF)

