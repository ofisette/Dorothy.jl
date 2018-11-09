using Dorothy
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

function xtcloop(traj::MolecularTrajectory, frame::MolecularModel, n::Integer)
	for i = 1:n
		seekstart(traj)
		while !eof(traj)
			read!(traj, frame)
		end
	end
end

traj = streamf(IOBuffer(read("$(datapath)/eMHC.xtc")))
frame = MolecularModel()
println("XTC loop (first):")
@time xtcloop(traj, frame, 1)
println("XTC loop (1):")
@time xtcloop(traj, frame, 1)
println("XTC loop (5):")
@time xtcloop(traj, frame, 5)
println("XTC loop (50):")
@time xtcloop(traj, frame, 50)
# println("XTC loop (500):")
# @time xtcloop(traj, frame, 500)
