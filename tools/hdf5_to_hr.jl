
using HDF5
#using DataFrames
#using CSV
using ArgParse
using Printf
using ProgressMeter

s = ArgParseSettings()
@add_arg_table s begin
    "--input", "-i"
        help = "Input file name.h5"
        arg_type = String
        default = "hamiltonians_pred.h5"         
    "--output", "-o"
        help = "Output file name _hr.dat"
        arg_type = String
        default = "Hamiltonian_hr.dat"    
    "--overlapnc"
        help = "Whether to convert the overlap matrix as noncollinear. If true, the overlap matrix will be converted to noncollinear form by duplicating the collinear overlap matrix for spin up and spin down. This is useful when the Hamiltonian is noncollinear but the overlap matrix is collinear."
        action = :store_true
        #default = false
    "--noncollinear"
        help = "Whether the Hamiltonian is noncollinear or not."
        action = :store_true # Stores true if the flag is present
        #default = false
end

parsed_args = parse_args(s)
filename = parsed_args["input"]
outputname = parsed_args["output"]
println("Input file: $filename")
noncollinear = 0
if parsed_args["noncollinear"] && parsed_args["overlapnc"]
    @warn "Both --noncollinear and --overlapnc options are provided."
end

if parsed_args["noncollinear"] || parsed_args["overlapnc"]
    noncollinear = 2    
else
    noncollinear = 1
end

if parsed_args["noncollinear"]
    println("Noncollinear Hamiltonian will be processed.")
else
    println("Collinear/Unpoloarized Hamiltonian will be processed.")
end

if parsed_args["overlapnc"]
    println("Overlap matrix will be converted to noncollinear form by duplicating the collinear overlap matrix for spin up and spin down.")
    println("This option is exclusive with --noncollinear.")
else
    println("Overlap matrix will be processed as spinless.")
end



function get_start_orbital_index(;orbitals = "orbital_types.dat")
    
    start_orbital_index = Int[0]
    open(orbitals, "r") do f
        for line in eachline(f)
            num_orbital_per_atom = 0
            for each_orbital in split(line)
                l = parse(Int, each_orbital)
                #print(num_orbital_per_atom, " ", each_orbital, " ", each_orbital *2 +1, "\n")
                num_orbital_per_atom += (2*l +1 )*noncollinear
            end
            push!(start_orbital_index, start_orbital_index[end] + num_orbital_per_atom)
        end
    end
    return start_orbital_index
end



@inline function get_orbital_index(start_orbital_index, atom1, atom2, index1, index2)

    orbital1 = start_orbital_index[atom1] + index1
    orbital2 = start_orbital_index[atom2] + index2

    return orbital1, orbital2
end

start_orbital_index = get_start_orbital_index()
println("Start orbital index: ", start_orbital_index)
println("Total number of orbitals: ", start_orbital_index[end])
println("Total number of atoms: ", length(start_orbital_index)-1)

h5open(filename, "r") do fin
open(outputname, "w") do fout

    Rvecs = Set{NTuple{3,Int}}()
    write(fout, "Converted # Rx Ry Rz Orb1 Orb2 H_real H_imag\n")

    norbital = start_orbital_index[end]
    #natom = 0 
    @info "Processing RVec and num_orbit..."
    for (key, obj) in pairs(fin)
        clean_str = strip(key, ['[', ']'])
        key = parse.(Int, split(clean_str, ","))
        push!(Rvecs, (key[1], key[2], key[3]))
    end
    println("Number of orbitals: $norbital")
    println("Number of Rvecs: ", length(Rvecs))
    write(fout, "$norbital\n")
    write(fout, "$(length(Rvecs))")
    for i in 1:length(Rvecs)
        if i % 15 == 1 
            write(fout, "\n")
        end
        write(fout, "1 ")
    end
    write(fout, "\n")

    #txt = ""

    @showprogress for (key, obj) in pairs(fin)
        
        clean_str = strip(key, ['[', ']'])
        key = parse.(Int, split(clean_str, ","))
        push!(Rvecs, (key[1], key[2], key[3]))
        if obj isa HDF5.Dataset
            data = read(obj)
            (ORB2,ORB1) = size(data)
            #println("Atoms: ", key[4], " ", key[5], " size: ", ORB1, " ", ORB2)
            @inbounds for j in 1:ORB2, i in 1:ORB1
                    J = data[j,i]
                    #key[4] and key[5] are atoms. 
                    orbital1, orbital2 = get_orbital_index(start_orbital_index,key[4], key[5], i, j)
                    #println("atom1: ", key[4], " atom2: ", key[5], " orbital1: ", orbital1, " orbital2: ", orbital2, " J: ", J)
                    #txt *= @sprintf("%d %d %d %d %d %.8f %.8f\n", key[1], key[2], key[3], orbital1, orbital2, real(J), imag(J))
                    #if orbital1 > norbital || orbital2 > norbital
                    #    @error "Orbital index exceeds the total number of orbitals. Check the orbital_types.dat file and the input Hamiltonian."
                    #    println(key[1], key[2], key[3], " atom1: ", key[4], " atom2: ", key[5], " i ", i, " j ", j, " orbital1: ", orbital1, " orbital2: ", orbital2, " J: ", J)
                    #    continue
                    #end
                    @printf(fout, "%d %d %d %d %d  %.8f %.8f\n", key[1], key[2], key[3],  orbital1, orbital2, real(J), imag(J))
                    if parsed_args["overlapnc"]
                        # If the overlap matrix is collinear but the Hamiltonian is noncollinear, we need to duplicate the overlap matrix for spin up and spin down.
                        # This means that for each entry in the collinear overlap matrix, we will write two entries in the output file: one for spin up and one for spin down.
                        # ORB1 and ORB2 are the number of orbitals for atom1 and atom2 without spin. 
                        @printf(fout, "%d %d %d %d %d  %.8f %.8f\n", key[1], key[2], key[3],  orbital1+ORB1 , orbital2+ORB2,  real(J), imag(J))
                    end
            end

            # This works only for the case one orbital per atom. 
            #J = data[1,1]            
            #txt *= @sprintf("%d %d %d %d %d %.8f %.8f\n", key[1], key[2], key[3], key[4], key[5], real(J), imag(J))

        end #if
    end # for



    #write(fout, txt)

end
end

println("Output file: $outputname")
