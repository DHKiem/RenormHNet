# Author: D.H. Kiem
# Last update: 2025-Nov-05

using HDF5
using LinearAlgebra
using CSV
using DataFrames
using ArgParse
using DelimitedFiles


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input_dir", "-i"
            help = ""
            arg_type = String
            default = "./"
        "--output_dir", "-o"
            help = ""
            arg_type = String
            default = "./output"
    end
    return parse_args(s)
end

#        "--mode"
#            help = ""
#            arg_type = String
#            default = "nothing"

function get_num_atoms(;elementfile = "element.dat")
    n = countlines(elementfile)
    return n 
end

function get_lattice(;latticefile = "lat.dat")
    lat = readdlm(latticefile)
    return lat
end

function get_positions(;positionfile = "site_positions.dat")
    position = readdlm(positionfile)
    return position
end


function symmed_spin_matrix_to_hdf5(targetfile::String)
    DF2 = CSV.read(targetfile, DataFrames.DataFrame)
    select!(DF2, [:Rx, :Ry, :Rz, :Atom1, :Atom2, :J])

    Atom_list = unique(DF2[:, :Atom1])
    nonmag_list = filter(x -> !(x in Atom_list), collect(1:get_num_atoms()))
    # to avoid error and guarantee the nearest site is itself.
    #for atom in Atom_list
    for atom in 1:get_num_atoms()
        #println(atom)
        push!(DF2, (Rx=0, Ry=0, Rz=0 ,Atom1=atom, Atom2=atom, J=0))
    end
    for atom1 in Atom_list, atom2 in nonmag_list
        for rxi in minimum(DF2.Rx):maximum(DF2.Rx), ryi in minimum(DF2.Ry):maximum(DF2.Ry), rzi in minimum(DF2.Rz):maximum(DF2.Rz)
            push!(DF2, (Rx=rxi, Ry=ryi, Rz=rzi ,Atom1=atom1, Atom2=atom2, J=0))
            push!(DF2, (Rx=rxi, Ry=ryi, Rz=rzi ,Atom1=atom2, Atom2=atom1, J=0))
        end
    end
    for atom1 in nonmag_list, atom2 in nonmag_list
        for rxi in minimum(DF2.Rx):maximum(DF2.Rx), ryi in minimum(DF2.Ry):maximum(DF2.Ry), rzi in minimum(DF2.Rz):maximum(DF2.Rz)
            if atom1 == atom2 
                continue
            end
            push!(DF2, (Rx=rxi, Ry=ryi, Rz=rzi ,Atom1=atom1, Atom2=atom2, J=0))
        end
    end    

    #R_list = Set{Vector{Int64}}() 
    R_list = unique(DF2[:, [:Rx, :Ry, :Rz]])

    #println(Rlist)

    open("R_list.dat", "w") do f
        for row in eachrow(R_list)
            println(f, join(row, " "))
        end
    end

    spinful = false
    hamiltonians = Dict{Array{Int64, 1}, Array{Float64, 2}}()


    g = groupby(DF2, [:Rx, :Ry, :Rz, :Atom1, :Atom2])


    for sub in g
        #println(sub)
        M = reshape(sub.J, 1, 1)

        key = [sub.Rx[1], sub.Ry[1], sub.Rz[1], sub.Atom1[1], sub.Atom2[1] ]
        hamiltonians[key] = M
    end
    return hamiltonians
    # save to hdf5 format 


end

function symmed_spin_matrix_to_hdf5_mntesb(targetfile::String)
    DF2 = CSV.read(targetfile, DataFrames.DataFrame)

    Atom_list = unique(DF2[:, :Atom1])
    nonmag_list = filter(x -> !(x in Atom_list), collect(1:get_num_atoms()))
    # to avoid error and guarantee the nearest site is itself.
    #for atom in Atom_list
    #for atom in 1:get_num_atoms()
    
    for atom in Atom_list
        println(atom)
        push!(DF2, (Rx=0, Ry=0, Rz=0 ,Atom1=atom, Atom2=atom, J=0))
    end
    

    #=
    for atom in 1:get_num_atoms()
        println(atom)
        push!(DF2, (Rx=0, Ry=0, Rz=0 ,Atom1=atom, Atom2=atom, J=0))
    end
    =#
    
    
    #=
    for atom1 in Atom_list, atom2 in nonmag_list
        for rxi in minimum(DF2.Rx):maximum(DF2.Rx), ryi in minimum(DF2.Ry):maximum(DF2.Ry), rzi in minimum(DF2.Rz):maximum(DF2.Rz)
            push!(DF2, (Rx=rxi, Ry=ryi, Rz=rzi ,Atom1=atom1, Atom2=atom2, J=0))
            push!(DF2, (Rx=rxi, Ry=ryi, Rz=rzi ,Atom1=atom2, Atom2=atom1, J=0))
        end
    end
    for atom1 in nonmag_list, atom2 in nonmag_list
        for rxi in minimum(DF2.Rx):maximum(DF2.Rx), ryi in minimum(DF2.Ry):maximum(DF2.Ry), rzi in minimum(DF2.Rz):maximum(DF2.Rz)
            if atom1 == atom2 
                continue
            end
            push!(DF2, (Rx=rxi, Ry=ryi, Rz=rzi ,Atom1=atom1, Atom2=atom2, J=0))
        end
    end    
    =#
    

    #R_list = Set{Vector{Int64}}() 
    R_list = unique(DF2[:, [:Rx, :Ry, :Rz]])

    #println(Rlist)

    open("R_list.dat", "w") do f
        for row in eachrow(R_list)
            println(f, join(row, " "))
        end
    end

    spinful = false
    hamiltonians = Dict{Array{Int64, 1}, Array{Float64, 2}}()


    g = groupby(DF2, [:Rx, :Ry, :Rz, :Atom1, :Atom2])


    for sub in g
        #println(sub)
        M = reshape(sub.J, 1, 1)

        key = [sub.Rx[1], sub.Ry[1], sub.Rz[1], sub.Atom1[1], sub.Atom2[1] ]
        hamiltonians[key] = M
    end
    return hamiltonians
    # save to hdf5 format 


end




parsed_args = parse_commandline()
parsed_args["input_dir"] = abspath(parsed_args["input_dir"])
cd(parsed_args["output_dir"])


targetfile = joinpath(parsed_args["input_dir"], "Spin_Matrix.csv")         
hamiltonians = symmed_spin_matrix_to_hdf5(targetfile)
#if parsed_args["mode"] == "skipte"
#    hamiltonians = symmed_spin_matrix_to_hdf5_mntesb(targetfile)    
#else
#    hamiltonians = symmed_spin_matrix_to_hdf5(targetfile)
#end
#println(hamiltonians)
h5open("Jhamiltonians.h5", "w") do fid
    for (key, hamiltonian) in hamiltonians
        write(fid, string(key), permutedims(hamiltonian))
    end
end