### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 2e27f842-b928-11ed-29ac-d1db772859fc
begin
	import Pkg
	Pkg.add("CairoMakie")
	Pkg.add("CSV")
	Pkg.add("DataFrames")
	using CairoMakie
	using CSV
	using DataFrames
end

# ╔═╡ 906b08fa-da9a-4e2f-8096-3ee44d6a0b64
r = (1:.2:3)

# ╔═╡ 4fb9d49a-6b96-4fa9-99ac-b2ba7c909527
begin
	# Read in all the data
	images = []
	for i in 1:length(r)
		filename = "result_r$(i).csv"
		df = DataFrame(CSV.File(filename))
		push!(images, Matrix(df))
	end

	dx = 0.05
	dt = 0.005
end

# ╔═╡ a539e3fc-402d-40f6-b95e-08e6f910232f
begin
	test = zeros(400,11)
	for i in 1:11
		test[:,i] .= images[i][:,1500]
	end
end

# ╔═╡ d07fe0ad-647b-4f56-8f59-2460d47c2c2e
begin
	f = Figure(resolution=(800,700), fontsize=20, font = "Serif")
	x = dx * (1:400)
	t = dt * (1:2000)
	
	ax1 = Axis(f[1,1], xlabel="Position", ylabel="Time", title="ϵ=0.15") 
	ax2 = Axis(f[1,2], xlabel="Position", ylabel="Time", title="ϵ=0.1")
	ax3 = Axis(f[2,1], xlabel="Position", ylabel="Arrival Rate", title="Timeslice")
	ax4 = Axis(f[2,2], xlabel="Position", ylabel="Time", title="ϵ=0.05")
	l1 = lines!(ax3, x, test[:,11])
	l2 = lines!(ax3, x, test[:,6])
	l3 = lines!(ax3, x, test[:,1])
	xlims!(ax3,0,x[end])
	axislegend(ax3, [l1,l2,l3],["ϵ=0.15","ϵ=0.10","ϵ=0.05"])

	hm1 = heatmap!(ax1, x[1:2:end], t[1:2:end], images[11][1:2:end,1:2:end], colorrange=(0,5))
	hm2 = heatmap!(ax2, x[1:2:end], t[1:2:end], images[6][1:2:end,1:2:end], colorrange=(0,5))
	hm3 = heatmap!(ax4, x[1:2:end], t[1:2:end], images[1][1:2:end,1:2:end], colorrange=(0,5))

	lines!(ax1, [0,20],[1500*dt,1500*dt], color=:red, linestyle="--")
	lines!(ax2, [0,20],[1500*dt,1500*dt], color=:red, linestyle="--")
	lines!(ax4, [0,20],[1500*dt,1500*dt], color=:red, linestyle="--")

	Colorbar(f[:,3], hm1, label = "Arrival Rate")

	Label(f[1,1,Bottom()], "(a)", padding = (0,0,10,60))
	Label(f[1,2,Bottom()], "(b)", padding = (0,0,10,60))
	Label(f[2,1,Bottom()], "(c)", padding = (0,0,10,60))
	Label(f[2,2,Bottom()], "(d)", padding = (0,0,10,60))

	f
end

# ╔═╡ 9b7e4281-a055-4c1f-89b8-61c263f3aa81
size(images[1])

# ╔═╡ 339007d0-87a5-4abb-b1dc-4290e306a1ca
size(x)

# ╔═╡ 7632a877-69ec-4dcc-a255-7711e4df66e7
Makie.to_font("Serif")

# ╔═╡ de1e3f2e-7181-45a3-b288-01c3b62d921c
save("converge.pdf",f)

# ╔═╡ Cell order:
# ╠═2e27f842-b928-11ed-29ac-d1db772859fc
# ╠═906b08fa-da9a-4e2f-8096-3ee44d6a0b64
# ╠═4fb9d49a-6b96-4fa9-99ac-b2ba7c909527
# ╠═a539e3fc-402d-40f6-b95e-08e6f910232f
# ╠═d07fe0ad-647b-4f56-8f59-2460d47c2c2e
# ╠═9b7e4281-a055-4c1f-89b8-61c263f3aa81
# ╠═339007d0-87a5-4abb-b1dc-4290e306a1ca
# ╠═7632a877-69ec-4dcc-a255-7711e4df66e7
# ╠═de1e3f2e-7181-45a3-b288-01c3b62d921c
