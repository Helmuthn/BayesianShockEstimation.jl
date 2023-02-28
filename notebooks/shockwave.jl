### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 981bb432-8ba4-11ed-1108-bd491a5d58cc
begin
	import Pkg
	Pkg.add(url="https://github.com/Helmuthn/naumer_Dimensionality_2022.jl")
	Pkg.add("CairoMakie")
	using BayesianShockEstimation
	using CairoMakie
end

# ╔═╡ 2c9c65ce-375a-4eb9-a896-6b650822de31
begin
	# Generate Boundary
	σ_w = 0.2
	σ_v = 1
	spatial_steps = 400

	boundary = BoundaryParams(σ_w, σ_v)
	
	true_boundary = boundary.σ_w * randn(spatial_steps)
	for i in 2:length(true_boundary)
		true_boundary[i] += true_boundary[i-1]
	end
	true_boundary .-= sum(true_boundary)/spatial_steps

	observations = true_boundary + boundary.σ_v * randn(spatial_steps)
end;

# ╔═╡ e85f2a4c-28f6-488b-b5f5-0890d59edd5e
begin
	# Configure the bigger simulations
	stepsize_x = 0.05
	stepsize_t = 0.005
	stepcount_f = 2000
	ballsize = 0.051
	threshold = 1
	
	shock_params = ShockParams( stepsize_x,
							  	stepsize_t,
								stepcount_f,
								ballsize,
								threshold);



	sample_count = 50
end;

# ╔═╡ 632b23f0-750c-408f-bbf7-ec0a3eca33f8
true_solution = godunov_burgers_1D( true_boundary, 
									stepsize_x, 
									stepsize_t,
									stepcount_f);

# ╔═╡ fde3b39b-f18f-4949-b3ba-81461d895c1c
begin
	# Generate Samples from RTS Sampler
	rts_sample_count = 3
	estimates, kalman_variances, rts_variances = RTS_Smooth(observations, σ_w, σ_v)
	samples = zeros(length(true_boundary),rts_sample_count)
	for i in 1:rts_sample_count
		samples[:,i] = RTS_sample(estimates, kalman_variances, rts_variances, σ_w)
	end
end;

# ╔═╡ 634a0428-685c-4d43-8617-f3759e7b423d
begin
	thread_count = Threads.nthreads()
	density_multi = zeros(length(observations), stepcount_f, thread_count)
	Threads.@threads for i in 1:thread_count
		density_multi[:,:,i] .= ShockDensity(shock_params, boundary, observations, sample_count);
	end
	density = dropdims(sum(density_multi, dims=3), dims=3) / thread_count
end;

# ╔═╡ 340063eb-d478-4645-840c-c1cd7ccfb2e8
begin
	downsample_rate = 4
	
	figa = Figure(resolution = (800,700),fontsize=20, font="Serif")
	ax1a = Axis(figa[2,2], xlabel="Space", ylabel="Time", title="Shock Arrival Rate")
	hma = heatmap!(ax1a, stepsize_x * (1:downsample_rate:spatial_steps-1), stepsize_t * (1:downsample_rate:stepcount_f-1), density[1:downsample_rate:end,1:downsample_rate:end], colorrange=(0,2))
	Colorbar(figa[2,3], hma, label="Arrival Rate")
	xlims!(ax1a, 0, stepsize_x * spatial_steps)
	ylims!(ax1a, 0, stepsize_t * stepcount_f)
	Label(figa[1,2,Bottom()], "(b)", padding = (0,0,10,60))


	downsample_rate = 2
	ax2a = Axis(figa[1,2], xlabel="Space", ylabel="Time", title="True Solution")
	hma2 = heatmap!(ax2a, stepsize_x * (1:downsample_rate:spatial_steps), stepsize_t * (1:downsample_rate:stepcount_f), true_solution[1:downsample_rate:end,1:downsample_rate:end])
	Colorbar(figa[1,3], hma2, label="Flowrate")
	xlims!(ax2a, 0, stepsize_x * spatial_steps)
	ylims!(ax2a, 0, stepsize_t * stepcount_f)
	Label(figa[2,2,Bottom()], "(d)", padding = (0,0,0,60))


	ax3a = Axis(figa[1,1], xlabel="Position", ylabel="Flowrate", title="Initial Condition")
	lines!(ax3a, stepsize_x * (1:length(observations)), observations, label="Observations")
	lines!(ax3a, stepsize_x * (1:length(observations)), true_boundary, label="Truth")
	xlims!(ax3a, 0, stepsize_x * length(observations))
	ylims!(ax3a, minimum(observations) - 2, maximum(observations) + .5)
	axislegend(position=:lb)
	Label(figa[1,1,Bottom()], "(a)", padding = (0,0,10,60))
	
	ax4a = Axis(figa[2,1], xlabel="Position", ylabel="Flowrate", title="$(rts_sample_count) Example RTS Samples")
	for i in 1:rts_sample_count
		lines!(ax4a,stepsize_x * (1:length(observations)), samples[:,i])
	end
	xlims!(ax4a, 0, stepsize_x * length(observations))
	Label(figa[2,1,Bottom()], "(c)", padding = (0,0,0,60))

	figa
end

# ╔═╡ e68830f0-5d6c-4074-a7a8-b036263dacc3
save("sspfig.pdf", figa)

# ╔═╡ Cell order:
# ╠═981bb432-8ba4-11ed-1108-bd491a5d58cc
# ╠═2c9c65ce-375a-4eb9-a896-6b650822de31
# ╠═e85f2a4c-28f6-488b-b5f5-0890d59edd5e
# ╠═632b23f0-750c-408f-bbf7-ec0a3eca33f8
# ╠═fde3b39b-f18f-4949-b3ba-81461d895c1c
# ╠═634a0428-685c-4d43-8617-f3759e7b423d
# ╠═340063eb-d478-4645-840c-c1cd7ccfb2e8
# ╠═e68830f0-5d6c-4074-a7a8-b036263dacc3
