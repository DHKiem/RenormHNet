#using StaticArrays
using LinearAlgebra
#using HDF5
using JSON
using DelimitedFiles
#using Statistics
using ArgParse

#= cmd from preprocess.py 
This part is a subprocess in preprocess.py
julia --project={project} openmx_get_data_for_jx.jl 
--input_dir ${raw}
--output_dir ${preprocessed} 
--save_overlap false

get_data function will be run first. 

patch
read out not scfout, but not writing orbital_types.dat, R_list.dat, hamiltonian.h5


=#

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
        "--if_DM", "-d"
            help = ""
            arg_type = Bool
            default = false
       "--save_overlap", "-s"
            help = ""
            arg_type = Bool
            default = false
    end
    return parse_args(s)
end
parsed_args = parse_commandline()


# @info string("get data from: ", parsed_args["input_dir"])
periodic_table = JSON.parsefile(joinpath(@__DIR__, "periodic_table.json"))


# The function parse_openmx below is come from "https://github.com/HopTB/HopTB.jl"
# parse_openmx has been also reedited from DeepH-Pack to read only .out instead of .scfout
function parse_openmx(filepath::String; return_DM::Bool = false)
    # define some helper functions for mixed structure of OpenMX binary data file.
    #=
    #function multiread(::Type{T}, f, size)::Vector{T} where T
    #    ret = Vector{T}(undef, size)
    #    read!(f, ret);ret
    #end
    #multiread(f, size) = multiread(Int32, f, size)

    #function read_mixed_matrix(::Type{T}, f, dims::Vector{<:Integer}) where T
    #    ret::Vector{Vector{T}} = []
    #    for i = dims; t = Vector{T}(undef, i);read!(f, t);push!(ret, t); end; ret
    #end

    function read_matrix_in_mixed_matrix(::Type{T}, f, spins, atomnum, FNAN, natn, Total_NumOrbs) where T
        ret = Vector{Vector{Vector{Matrix{T}}}}(undef, spins)
        for spin = 1:spins;t_spin = Vector{Vector{Matrix{T}}}(undef, atomnum)
            for ai = 1:atomnum;t_ai = Vector{Matrix{T}}(undef, FNAN[ai])
                for aj_inner = 1:FNAN[ai]
                    t = Matrix{T}(undef, Total_NumOrbs[natn[ai][aj_inner]], Total_NumOrbs[ai])
                    read!(f, t);t_ai[aj_inner] = t
                end;t_spin[ai] = t_ai
            end;ret[spin] = t_spin
        end;return ret
    end
    read_matrix_in_mixed_matrix(f, spins, atomnum, FNAN, natn, Total_NumOrbs) = read_matrix_in_mixed_matrix(Float64, f, spins, atomnum, FNAN, natn, Total_NumOrbs)
    
    read_3d_vecs(::Type{T}, f, num) where T = reshape(multiread(T, f, 4 * num), 4, Int(num))[2:4,:]
    read_3d_vecs(f, num) = read_3d_vecs(Float64, f, num)
    # End of helper functions

    bound_multiread(T, size) = multiread(T, f, size)
    bound_multiread(size) = multiread(f, size)
    bound_read_mixed_matrix() = read_mixed_matrix(Int32, f, FNAN)
    bound_read_matrix_in_mixed_matrix(spins) = read_matrix_in_mixed_matrix(f, spins, atomnum, FNAN, natn, Total_NumOrbs)
    bound_read_3d_vecs(num) = read_3d_vecs(f, num)
    bound_read_3d_vecs(::Type{T}, num) where T = read_3d_vecs(T, f, num)
    # End of bound helper functions
    for binary
    =#
    f = open(filepath)
    
    #=
    atomnum, SpinP_switch, Catomnum, Latomnum, Ratomnum, TCpyCell, order_max = bound_multiread(7)
    println("atomnum: ", atomnum) 
    println("SpinP_switch: ", SpinP_switch) 
    
    @assert (SpinP_switch >> 2) == 3 "DeepH-pack only supports OpenMX v3.9. Please check your OpenMX version"
    SpinP_switch &= 0x03

    atv, atv_ijk = bound_read_3d_vecs.([Float64,Int32], TCpyCell + 1)

    Total_NumOrbs, FNAN = bound_multiread.([atomnum,atomnum])
    FNAN .+= 1
    natn = bound_read_mixed_matrix()
    ncn = ((x)->x .+ 1).(bound_read_mixed_matrix()) # These is to fix that atv and atv_ijk starts from 0 in original C code.

    tv, rtv, Gxyz = bound_read_3d_vecs.([3,3,atomnum])

    Hk = bound_read_matrix_in_mixed_matrix(SpinP_switch + 1)
    iHk = SpinP_switch == 3 ? bound_read_matrix_in_mixed_matrix(3) : nothing
    OLP = bound_read_matrix_in_mixed_matrix(1)[1]
    OLP_r = []
    for dir in 1:3, order in 1:order_max
        t = bound_read_matrix_in_mixed_matrix(1)[1]
        if order == 1 push!(OLP_r, t) end
    end
    OLP_p = bound_read_matrix_in_mixed_matrix(3)
    DM = bound_read_matrix_in_mixed_matrix(SpinP_switch + 1)
    iDM = bound_read_matrix_in_mixed_matrix(2)
    solver = bound_multiread(1)[1]
    chem_p, E_temp = bound_multiread(Float64, 2)
    dipole_moment_core, dipole_moment_background = bound_multiread.(Float64, [3,3])
    Valence_Electrons, Total_SpinS = bound_multiread(Float64, 2)
    dummy_blocks = bound_multiread(1)[1]
    for i in 1:dummy_blocks
        bound_multiread(UInt8, 256)
    end
    #### Above the part of that read from scfout. Which will not be required in the end.
    =#

    lines = readlines(filepath)

    function findkey(key)
        k = lowercase(key)
        for l in lines
            w = split(l)
            if !isempty(w) && lowercase(w[1]) == k
                return length(w) >= 2 ? w[2] : nothing
            end
        end
        return nothing
    end

    atomnum_str = findkey("Atoms.Number")
    atomnum_str === nothing && error("Atoms.Number not found in $filepath")
    atomnum = parse(Int, atomnum_str)
    unit    = something(findkey("Atoms.UnitVectors.Unit"), "Ang")

    i0 = findfirst(l -> strip(l) == "<Atoms.UnitVectors", lines)

    tv = zeros(3, 3)
    for k in 1:3
        tv[:, k] = parse.(Float64, split(lines[i0 + k])[1:3])
    end
    lowercase(unit) == "au" || (tv ./= 0.529177249)   
    
    f = open(filepath)   


    # jx mode does not read .scfout but only read .out 
    function strip1(s::Vector{UInt8})
        startpos = 0
        for i = 1:length(s) + 1
            if i > length(s) || s[i] & 0x80 != 0 || !isspace(Char(s[i] & 0x7f))
                startpos = i
                break
            end
        end
        return s[startpos:end]
    end
    function startswith1(s::Vector{UInt8}, prefix::Vector{UInt8})
        return length(s) >= length(prefix) && s[1:length(prefix)] == prefix
    end
    function isnum(s::Char)
        if s >= '1' && s <= '9'
            return true
        else
            return false
        end
    end

    function isorb(s::Char)
        if s in ['s','p','d','f']
            return true
        else
            return false
        end
    end

    function orbital_types_str2num(str::AbstractString)
        orbs = split(str, isnum, keepempty = false)
        nums = map(x->parse(Int, x), split(str, isorb, keepempty = false))
        orb2l = Dict("s" => 0, "p" => 1, "d" => 2, "f" => 3)
        @assert length(orbs) == length(nums)
        orbital_types = Array{Int64,1}(undef, 0)
        for (orb, num) in zip(orbs, nums)
            for i = 1:num
                push!(orbital_types, orb2l[orb])
            end
        end
        return orbital_types
    end

    function find_target_line(target_line::String)
        target_line_UInt8 = Vector{UInt8}(target_line)
        while !startswith1(strip1(Vector{UInt8}(readline(f))), target_line_UInt8)
            if eof(f)
                error(string(target_line, "not found. Please check if the .out file was appended to the end of .scfout file!"))
            end
        end
    end


    
#     @info "get Chemical potential (Hartree)"
    find_target_line("(see also PRB 72, 045121(2005) for the energy contributions)")
    readline(f)
    readline(f)
    readline(f)
    str = split(readline(f))
    @assert "Chemical" == str[1]
    @assert "potential" == str[2]
    @assert "(Hartree)" == str[3]
    ev2Hartree = 0.036749324533634074
    fermi_level = parse(Float64, str[length(str)])/ev2Hartree

    # @info "get Chemical potential (Hartree)"
    # find_target_line("Eigenvalues (Hartree)")
    # for i = 1:2;@assert readline(f) == "***********************************************************";end
    # readline(f)
    # str = split(readline(f))
    # ev2Hartree = 0.036749324533634074
    # fermi_level = parse(Float64, str[length(str)])/ev2Hartree
    

#     @info "get Fractional coordinates & orbital types:"
    find_target_line("Fractional coordinates of the final structure")
    target_line = Vector{UInt8}("Fractional coordinates of the final structure")
    for i = 1:2;@assert readline(f) == "***********************************************************";end
    @assert readline(f) == ""
    orbital_types = Array{Array{Int64,1},1}(undef, 0) #orbital_types
    element = Array{Int64,1}(undef, 0) #orbital_types
    atom_frac_pos = zeros(3, atomnum) #Fractional coordinates
    for i = 1:atomnum
        str = readline(f)
        element_str = split(str)[2]
        #push!(orbital_types, orbital_types_element[element_str]) #original
        push!(orbital_types, [0]) # for jx Heisenberg
        m = match(r"^\s*\d+\s+\w+\s+([0-9+-.Ee]+)\s+([0-9+-.Ee]+)\s+([0-9+-.Ee]+)", str)
        push!(element, periodic_table[element_str]["Atomic no"])
        atom_frac_pos[:,i] = ((x)->parse(Float64, x)).(m.captures)
    end
    atom_pos = tv * atom_frac_pos
    close(f)

    # use the atom_pos to fix
    # TODO: Persuade wangc to accept the following code, which seems hopeless and meaningless.
    """
    for axis = 1:3
        ((x2, y2, z)->((x, y)->x .+= z * y).(x2, y2)).(OLP_r[axis], OLP, atom_pos[axis,:])
    end
    """

    #return element, atomnum, atv, atv_ijk,  FNAN, natn, ncn, tv, orbital_types, fermi_level, atom_pos
    return element, atomnum, tv, orbital_types, fermi_level, atom_pos
end #parse_openmx

function get_data(filepath_out::String, Rcut::Float64; if_DM::Bool = false)
    #element, nsites, SpinP_switch, atv, atv_ijk, Total_NumOrbs, FNAN, natn, ncn, lat, Hk, iHk, OLP, OLP_r, orbital_types, fermi_level, site_positions, DM = parse_openmx(filepath_scfout; return_DM=if_DM)
    element, nsites, lat, orbital_types, fermi_level, site_positions = parse_openmx(filepath_out; return_DM=if_DM)

    site_positions .*= 0.529177249 # Bohr to Ang
    lat .*= 0.529177249 # Bohr to Ang
    
    # get neighbour_site
    nsites = size(site_positions, 2)

    spinful = false

    #return element, overlaps, density_matrixs, hamiltonians, fermi_level, orbital_types, lat, site_positions, spinful, R_list
    return element, fermi_level, orbital_types, lat, site_positions, spinful
end # get_data

parsed_args["input_dir"] = abspath(parsed_args["input_dir"])
mkpath(parsed_args["output_dir"])
cd(parsed_args["output_dir"])

#element, overlaps, density_matrixs, hamiltonians, fermi_level, orbital_types, lat, site_positions, spinful, R_list = get_data(joinpath(parsed_args["input_dir"], "openmx.scfout"), -1.0; if_DM=parsed_args["if_DM"])
element,  fermi_level, orbital_types, lat, site_positions, spinful = get_data(joinpath(parsed_args["input_dir"], "openmx.out"), -1.0; if_DM=parsed_args["if_DM"])


info_dict = Dict(
    "fermi_level" => fermi_level,
    "isspinful" => spinful
    )
open("info.json", "w") do f
    JSON.print(f, info_dict, 4)
end
open("site_positions.dat", "w") do f
    writedlm(f, site_positions)
end

open("lat.dat", "w") do f
    writedlm(f, lat)
end
rlat = 2pi * inv(lat)'
open("rlat.dat", "w") do f
    writedlm(f, rlat)
end

open("orbital_types.dat", "w") do f
    writedlm(f, orbital_types)
end

open("element.dat", "w") do f
    writedlm(f, element)
end
