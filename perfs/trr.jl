using Dorothy
using Formats
using FormatStreams

datapath = joinpath(@__DIR__, "..", "data")

function trrloop(traj::MolecularTrajectory, frame::MolecularModel, n::Integer)
	for i = 1:n
		seekstart(traj)
		while !eof(traj)
			read!(traj, frame)
		end
	end
end

traj = streamf(IOBuffer(read("$(datapath)/eMHC.trr")))
frame = MolecularModel()
println("TRR loop (first):")
@time trrloop(traj, frame, 1)
println("TRR loop (1):")
@time trrloop(traj, frame, 1)
println("TRR loop (5):")
@time trrloop(traj, frame, 5)
println("TRR loop (50):")
@time trrloop(traj, frame, 50)
println("TRR loop (500):")
@time trrloop(traj, frame, 500)
