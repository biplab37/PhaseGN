mu_list = [0.0, 0.9, 1.0, 1.1]

srange = -2:0.01:2

param
pot_data = zeros(length(srange), length(mu_list))

for i in eachindex(srange)
	for j in eachindex(mu_list)
		pot_data[i, j] = omega(srange[i], 0.01, mu_list[j], param0)
	end
end

pot_data2 = zeros(length(srange), length(mu_list))

for i in eachindex(srange)
	for j in eachindex(mu_list)
		pot_data2[i, j] = omega(srange[i], 0.01, mu_list[j], param)
	end
end

writedlm("fig_2_a_grand_potential.csv", hcat(srange, pot_data) )
writedlm("fig_2_b_grand_potential.csv", hcat(srange, pot_data2) )