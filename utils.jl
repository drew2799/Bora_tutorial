using Printf
using LaTeXStrings

function get_zeff(zbin)
    if zbin == "z1"
        return 1.0027
    elseif zbin == "z2"
        return 1.2016
    elseif zbin == "z3"
        return 1.4007
    elseif zbin == "z4"
        return 1.65
    else
        throw(ArgumentError("Invalid zbin. It must be 'z1', 'z2', 'z3' or 'z4', but you gave $(zbin)."))
    end
end

function load_emu(path, zbin, ℓ, architecture, r_test)
    weights = npzread(path*"/"*string(ℓ)*"/weights.npy")
    trained_emu = Bora.SimpleChainsEmulator(Architecture=architecture, Weights=weights)
    ξℓ_emu = Bora.ξℓEmulator(TrainedEmulator = trained_emu, rgrid=r_test,
                                InMinMax = npzread(path*"/"*string(ℓ)*"/inminmax.npy"),
                                OutMinMax = npzread(path*"/"*string(ℓ)*"/outminmax.npy"))
    return ξℓ_emu
end

function LoadEmulator(emu_path, zbin, s)
    println("🔄 Loading BORA emulator for redshift "*string(zbin))
    ξℓ0_emu = load_emu(emu_path, zbin, 0, mlpd, s)
    ξℓ2_emu = load_emu(emu_path, zbin, 2, mlpd, s)
    ξℓ4_emu = load_emu(emu_path, zbin, 4, mlpd, s)
    complete_emu = Bora.CompleteEmulator(rgrid=s, ξℓMono=ξℓ0_emu, ξℓQuad=ξℓ2_emu, ξℓHexa=ξℓ4_emu);
    println("✅ Loading completed")
    return complete_emu
end

function load_emu(path, ℓ, architecture, r_test)
    weights = npzread(path*"/"*string(ℓ)*"/weights.npy")
    trained_emu = Bora.SimpleChainsEmulator(Architecture=architecture, Weights=weights)
    ξℓ_emu = Bora.ξℓEmulator(TrainedEmulator = trained_emu, rgrid=r_test,
                                InMinMax = npzread(path*"/"*string(ℓ)*"/inminmax.npy"),
                                OutMinMax = npzread(path*"/"*string(ℓ)*"/outminmax.npy"))
    return ξℓ_emu
end

function LoadEmulator(emu_path, s)
    @info "🔄 Loading BORA emulator"
    ξℓ0_emu = load_emu(emu_path, 0, mlpd, s)
    ξℓ2_emu = load_emu(emu_path, 2, mlpd, s)
    ξℓ4_emu = load_emu(emu_path, 4, mlpd, s)
    complete_emu = Bora.CompleteEmulator(rgrid=s, ξℓMono=ξℓ0_emu, ξℓQuad=ξℓ2_emu, ξℓHexa=ξℓ4_emu);
    @info "✅ Loading completed"
    return complete_emu
end

function Hartlap_factor(N_mocks, N_obs)
    return (N_mocks-N_obs-2)/(N_mocks-1)
end

function Percival_factor(N_mocks, N_obs, N_params)
    B = (N_mocks-N_obs-2)/((N_mocks-N_obs-1)*(N_mocks-N_obs-4))
    f = (N_mocks-1)*(1+B*(N_obs-N_params))/(N_mocks-N_obs+N_params-1)
    return f
end

function leave_1_out_cov(ξℓ_mocks, mock_id)

    mask = trues(size(ξℓ_mocks,1))
    if mock_id != 0
        mask[mock_id] = false
    end
    masked_ξℓs = ξℓ_mocks[mask,:,:]
    
    reshaped_ξℓ = PermutedDimsArray(masked_ξℓs, (1, 3, 2))
    ξℓs = hcat(reshaped_ξℓ[:,:,1], reshaped_ξℓ[:,:,2], reshaped_ξℓ[:,:,3])

    l1o_cov = cov(ξℓs, dims=1, corrected=true)

    return l1o_cov
end

function LoadData(mocks_path, s_min, s_max, N_params; ells=(0, 2, 4), mock_id=0)

    @info "🔄 Loading data 🔄"
    
    s = npzread(mocks_path*"separations_rebin_5_s_0_200.npy")
    ξℓ_mocks = npzread(mocks_path*"multipoles_rebin_5_s_0_200.npy")

    if mock_id==0
        print("Using the mean of 1000 mocks.")
        ξℓ_mean = mean(ξℓ_mocks, dims=1)
        ξ0_mean = ξℓ_mean[1,1,:]
        ξ2_mean = ξℓ_mean[1,2,:]
        ξ4_mean = ξℓ_mean[1,3,:]
        corr_f = 0
    elseif 1<=mock_id<=1000
        @printf("Using mock %04d.", mock_id)
        ξℓ_mean = ξℓ_mocks[mock_id,:,:]
        ξ0_mean = ξℓ_mean[1,:]
        ξ2_mean = ξℓ_mean[2,:]
        ξ4_mean = ξℓ_mean[3,:]
        corr_f = 1
    else
        throw(ArgumentError("Invalid mock_id. It must be 0 to use mean vector or in [1,1000], but you gave $(mock_id)."))
    end

    @info "✂️ Separation range cutting ✂️"

    nℓ = length(ells)
    ξ_list = (ξ0_mean, ξ2_mean, ξ4_mean)
    ξℓ_mean_flat = vcat(ξ_list[1:nℓ]...)

    mask = (s .>= s_min) .& (s .<= s_max)
    full_mask = repeat(mask, nℓ)

    cut_ξℓ_mean_flat = ξℓ_mean_flat[full_mask]

    @info "🧮 Covariance and Percival correction 🧮"
    
    N_mocks = size(ξℓ_mocks)[1]
    N_obs_mocks = size(ξℓ_mocks)[2]*size(ξℓ_mocks)[3]
    N_obs = length(cut_ξℓ_mean_flat)
    
    uncorr_cov = leave_1_out_cov(ξℓ_mocks, mock_id)
    cut_uncorr_cov = uncorr_cov[1:Int(nℓ*N_obs_mocks/3), 1:Int(nℓ*N_obs_mocks/3)]
    cut_uncorr_cov = uncorr_cov[full_mask, full_mask]
    
    α = Percival_factor(N_mocks-corr_f, N_obs, N_params)
    cut_cov = α .* cut_uncorr_cov

    @info "✅ Loading completed"
    
    return s[mask], cut_ξℓ_mean_flat, Hermitian(cut_cov)
end

function LoadData2(mocks_path, s_min, s_max, N_params; ells=(0, 2, 4), mock_id=0, correct=nothing)

    @info "🔄 Loading data 🔄"
    
    s = npzread(mocks_path*"separations_rebin_5_s_0_200.npy")
    ξℓ_mocks = npzread(mocks_path*"multipoles_rebin_5_s_0_200.npy")

    if mock_id==0
        print("Using the mean of 1000 mocks.")
        ξℓ_mean = mean(ξℓ_mocks, dims=1)
        ξ0_mean = ξℓ_mean[1,1,:]
        ξ2_mean = ξℓ_mean[1,2,:]
        ξ4_mean = ξℓ_mean[1,3,:]
    elseif 1<=mock_id<=1000
        @printf("Using mock %04d.", mock_id)
        ξℓ_mean = ξℓ_mocks[mock_id,:,:]
        ξ0_mean = ξℓ_mean[1,:]
        ξ2_mean = ξℓ_mean[2,:]
        ξ4_mean = ξℓ_mean[3,:]
    else
        throw(ArgumentError("Invalid mock_id. It must be 0 to use mean vector or in [1,1000], but you gave $(mock_id)."))
    end

    @info "✂️ Separation range cutting ✂️"

    nℓ = length(ells)
    ξ_list = (ξ0_mean, ξ2_mean, ξ4_mean)
    ξℓ_mean_flat = vcat(ξ_list[1:nℓ]...)

    mask = (s .>= s_min) .& (s .<= s_max)
    full_mask = repeat(mask, nℓ)

    cut_ξℓ_mean_flat = ξℓ_mean_flat[full_mask]

    @info "🧮 Covariance and Percival correction 🧮"
    
    N_mocks = size(ξℓ_mocks)[1]
    N_obs_mocks = size(ξℓ_mocks)[2]*size(ξℓ_mocks)[3]
    N_obs = length(cut_ξℓ_mean_flat)
    
    uncorr_cov = leave_1_out_cov(ξℓ_mocks, 0)
    cut_uncorr_cov = uncorr_cov[1:Int(nℓ*N_obs_mocks/3), 1:Int(nℓ*N_obs_mocks/3)]
    cut_uncorr_cov = uncorr_cov[full_mask, full_mask]

    if correct==nothing
        @info "✅ Loading completed"
        return s[mask], cut_ξℓ_mean_flat, Hermitian(cut_uncorr_cov)
    elseif correct=="hartlap"
        α = Hartlap_factor(N_mocks, N_obs)
        cut_cov = inv(α.*inv(cut_uncorr_cov))
    elseif correct=="percival"
        α = Percival_factor(N_mocks, N_obs, N_params)
        cut_cov = α .* cut_uncorr_cov
    else
        throw(ArgumentError("Invalid correction to covariance. It must be nothing, 'hartlap' or 'percival', but you gave $(correct)."))
    end

    @info "✅ Loading completed"
    return s[mask], cut_ξℓ_mean_flat, Hermitian(cut_cov)
end

function chisq(data, model, cov; nparams=0)
    @assert length(data) == length(model) "data and model size mismatch"
    @assert size(cov,1) == size(cov,2) == length(data) "covariance matrix size mismatch"

    r = data .- model
    chi2 = r' * inv(cov) * r

    dof = length(data) - nparams
    chi2_red = chi2 / dof

    return chi2, chi2_red
end

function sgrid_piecewise(N; edges=(0.0,50.0,80.0,120.0,150.0,200.0),
                            weights=(1.0,2.0,6.0,2.0,1.0))
    @assert length(edges) == length(weights) + 1

    L = diff(collect(edges))
    w = collect(weights)
    alloc = w .* L
    nseg = round.(Int, N .* alloc ./ sum(alloc))

    # correggi eventuale mismatch per arrotondamenti
    Δ = N - sum(nseg)
    if Δ != 0
        order = sortperm(alloc; rev=true)
        for k in 1:abs(Δ)
            nseg[order[(k-1) % length(nseg) + 1]] += sign(Δ)
        end
    end

    s = Float64[]
    for i in 1:length(nseg)
        a, b = edges[i], edges[i+1]
        ni = nseg[i]
        ni <= 0 && continue
        append!(s, range(a, b; length=ni+1)[1:end-1])  # evita duplicati ai bordi
    end
    
    return s
end

function perturb(nt::NamedTuple)
    NamedTuple{keys(nt)}(values(nt) .+ rand(Normal(0, 0.5), length(nt)))
end