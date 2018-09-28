using Profile
using Dorothy
using Formats

profilepath = @__DIR__
datapath = joinpath(@__DIR__, "..", "..", "data")

readf("$(datapath)/PLC.pdb")

@profile (for i = 1:2; readf("$(datapath)/PLC.pdb"); end)

open("$(profilepath)/pdbread.txt", "w") do io
    Profile.print(IOContext(io, :displaysize => (24, 500)))
end
