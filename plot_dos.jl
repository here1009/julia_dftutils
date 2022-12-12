#--types--
mutable struct DosInput
    is_partial_dos::Int64
    is_interpolation_method::Int64
    energy_smearing_eV::Float64
    number_of_dos_grid::Int64
    interpolation_grid::Vector{Int64}
    is_exist_occupation_file::Int64
end
#--functions--
function read_dos_input(f)
    open(f,"r") do file
        p1=parse(Int64,readline(file))
        p2=parse(Int64,readline(file))
        p3=parse(Float64,readline(file))
        p4=parse(Int64,readline(file))
        p5=parse.(Int64,rsplit(readline(file),[',',' ']))
        p6=parse(Int64,readline(file))
    end
    return DosInput(1,1,2.2,2,[2,1,3],1)
end

#--main progress--

in=read_dos_input("DOS.input")

println(in.is_partial_dos)
println(in.is_interpolation_method)
println(in.energy_smearing_eV)
println(in.number_of_dos_grid)
println(in.interpolation_grid)
println(in.is_exist_occupation_file)