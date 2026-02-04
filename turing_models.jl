using Turing


ℓ_to_ℓidx = Dict(0=>1, 2=>2, 4=>3)

@model function polyBB_BORAmodel(D, inv_Γ, emu, s, s_test, zbin, growth; ells=(0, 2, 4), iso_ap=false, cov=nothing)

    ℓidx = getindex.(Ref(ℓ_to_ℓidx), collect(ells))
    
    if iso_ap
        α̃ISO ~ 100*Uniform(0.8, 1.2)
        α̃AP ~ 100*Uniform(0.8, 1.2)
        αpar = (α̃ISO/100) * (α̃AP/100)^(2/3)
        αperp = (α̃ISO/100) * (α̃AP/100)^(-1/3)
    else
        α̃par ~ 100*Uniform(0.8, 1.2)
        α̃perp ~ 100*Uniform(0.8, 1.2)
        αpar = α̃par/100
        αperp = α̃perp/100
    end
    b̃ ~ 10*Uniform(0., 5.)  
    f̃ ~ 10*Uniform(0., 2.)
    Σpar ~ Uniform(0., 20.)
    Σperp ~ Uniform(0., 20.)
    Σs ~ Uniform(0., 10.)

    lim = 20
    
    B00 ~ Uniform(-lim, lim)
    B01 ~ Uniform(-lim, lim)
    B02 ~ Uniform(-lim, lim)
    if 2 in ells
        B20 ~ Uniform(-lim, lim)
        B21 ~ Uniform(-lim, lim)
        B22 ~ Uniform(-lim, lim)
    else
        B20 = 0
        B21 = 0
        B22 = 0
    end
    if 4 in ells
        B40 ~ Uniform(-lim, lim)
        B41 ~ Uniform(-lim, lim)
        B42 ~ Uniform(-lim, lim)
    else
        B40 = 0
        B41 = 0
        B42 = 0
    end

    bao_params = [αpar, αperp, b̃/10, f̃/10, Σpar, Σperp, Σs, get_zeff(zbin)]
    bb_params = [B00, B01, B02, B20, B21, B22, B40, B41, B42]

    theory = Bora.get_ξℓs(bao_params, emu)
    BB = Bora.get_broadband(s, bb_params)
    
    ξ_list = [(growth^2) .* DataInterpolations.AkimaInterpolation(view(theory, idx, :), s_test).(s) .+ view(BB, idx, :) for idx in ℓidx]
    interp_ξ = permutedims(hcat(ξ_list...))
    
    if cov == nothing
        flattheory = reshape((interp_ξ)', (size(inv_Γ)[1],))
        prediction = (inv_Γ * flattheory)
        D ~ MvNormal(prediction, I)
    else
        flattheory = reshape((interp_ξ)', (size(cov)[1],))
        D ~ MvNormal(flattheory, cov)
    end
    
    return nothing
end

@model function splineBB_BORAmodel(D, inv_Γ, emu, s, s_test, zbin, growth; ells=(0, 2), iso_ap=false)
    
    if iso_ap
        α̃ISO ~ 100*Uniform(0.8, 1.2)
        α̃AP ~ 100*Uniform(0.8, 1.2)
        αpar = (α̃ISO/100) * (α̃AP/100)^(2/3)
        αperp = (α̃ISO/100) * (α̃AP/100)^(-1/3)
    else
        α̃par ~ 100*Uniform(0.8, 1.2)
        α̃perp ~ 100*Uniform(0.8, 1.2)
        αpar = α̃par/100
        αperp = α̃perp/100
    end
    b̃ ~ 10*Uniform(0., 5.)  
    f̃ ~ 10*Uniform(0., 2.)
    Σpar ~ Uniform(0., 20.)
    Σperp ~ Uniform(0., 20.)
    Σs ~ Uniform(0., 10.)

    lim = 20
    
    a0 ~ Uniform(-lim, lim)
    a1 ~ Uniform(-lim, lim)
    a2 ~ Uniform(-lim, lim)
    a3 ~ Uniform(-lim, lim)
    a4 ~ Uniform(-lim, lim)
    a5 ~ Uniform(-lim, lim)
    if 4 in ells
        a6 ~ Uniform(-lim, lim)
        a7 ~ Uniform(-lim, lim)
        a8 ~ Uniform(-lim, lim)
    else
        a6 = 0
        a7 = 0
        a8 = 0
    end
    
    bao_params = [αpar, αperp, b̃/10, f̃/10, Σpar, Σperp, Σs, get_zeff(zbin)]
    bb_params = [a0*0.006, a1*0.006, a2*0.006, a3*0.006, a4*100, a5, a6, a7, a8]

    theory = Bora.get_ξℓs(bao_params, emu)

    BB = get_spline_broadband(s, bb_params; Δ=0.08, kmin=0.02, n_hexa=3)
    
    interp_ξ0 = DataInterpolations.AkimaInterpolation((growth^2).*theory[1,:], s_test).(s) .+ BB[1,:]
    interp_ξ2 = DataInterpolations.AkimaInterpolation((growth^2).*theory[2,:], s_test).(s) .+ BB[2,:]
    if 4 in ells
        interp_ξ4 = DataInterpolations.AkimaInterpolation((growth^2).*theory[3,:], s_test).(s) .+ BB[3,:]
        interp_ξ = Matrix(hcat(interp_ξ0, interp_ξ2, interp_ξ4)')
    else
        interp_ξ = Matrix(hcat(interp_ξ0, interp_ξ2)')
    end
    
    flattheory = reshape((interp_ξ)', (size(inv_Γ)[1],))
    
    prediction = (inv_Γ * flattheory)

    D ~ MvNormal(prediction, I)
    
    return nothing
end

function find_best_fit(posterior_model, initial_params, estimator; ntry=10)
    
    lp = -1.0e4
    BF = nothing
    status = nothing
    bf_best = nothing

    for i in 1:ntry
        
        init_p = rand(MvNormal(initial_params, I))
        
        if estimator == "MAP"
            bf = maximum_a_posteriori(posterior_model; initial_params=init_p)
        elseif estimator == "MLE"
            bf = maximum_likelihood(posterior_model; initial_params=init_p)
        else
            throw(ArgumentError("Invalid estimator. It must be 'MAP' or 'MLE', but you gave $(estimator)."))
        end
        
        if bf.lp > lp
            BF = bf.values.array
            lp = bf.lp
            status = bf.optim_result.retcode
            bf_best = bf
        end
    end

    return BF, lp, status, bf_best
end

function loglike_SH(x0::AbstractVector, μ::AbstractVector, inv_S, ns::Int, m::Real)
    r = x0 .- μ
    chi2 = dot(r, inv_S * r)
    return -(m/2) * log1p(chi2 / (ns - 1))
end

function m(ns, nd, nθ)
    B = (ns - nd - 2)/((ns - nd - 1)*(ns - nd - 4))
    return nθ + 2 + (ns - 1 + B*(nd - nθ))/(1 + B*(nd - nθ))
end

@model function polyBB_t_student_BORAmodel(d, inv_S, emu, s, s_test, zbin, growth, ns, nd, nθ; ells=(0, 2, 4), iso_ap=false, cov=nothing)

    ℓidx = getindex.(Ref(ℓ_to_ℓidx), collect(ells))
    
    if iso_ap
        α̃ISO ~ 100*Uniform(0.8, 1.2)
        α̃AP ~ 100*Uniform(0.8, 1.2)
        αpar = (α̃ISO/100) * (α̃AP/100)^(2/3)
        αperp = (α̃ISO/100) * (α̃AP/100)^(-1/3)
    else
        α̃par ~ 100*Uniform(0.8, 1.2)
        α̃perp ~ 100*Uniform(0.8, 1.2)
        αpar = α̃par/100
        αperp = α̃perp/100
    end
    b̃ ~ 10*Uniform(0., 5.)  
    f̃ ~ 10*Uniform(0., 2.)
    Σpar ~ Uniform(0., 20.)
    Σperp ~ Uniform(0., 20.)
    Σs ~ Uniform(0., 10.)

    lim = 20
    
    B00 ~ Uniform(-lim, lim)
    B01 ~ Uniform(-lim, lim)
    B02 ~ Uniform(-lim, lim)
    if 2 in ells
        B20 ~ Uniform(-lim, lim)
        B21 ~ Uniform(-lim, lim)
        B22 ~ Uniform(-lim, lim)
    else
        B20 = 0
        B21 = 0
        B22 = 0
    end
    if 4 in ells
        B40 ~ Uniform(-lim, lim)
        B41 ~ Uniform(-lim, lim)
        B42 ~ Uniform(-lim, lim)
    else
        B40 = 0
        B41 = 0
        B42 = 0
    end

    bao_params = [αpar, αperp, b̃/10, f̃/10, Σpar, Σperp, Σs, get_zeff(zbin)]
    bb_params = [B00, B01, B02, B20, B21, B22, B40, B41, B42]

    theory = Bora.get_ξℓs(bao_params, emu)
    BB = Bora.get_broadband(s, bb_params)
    
    ξ_list = [(growth^2) .* DataInterpolations.AkimaInterpolation(view(theory, idx, :), s_test).(s) .+ view(BB, idx, :) for idx in ℓidx]
    interp_ξ = permutedims(hcat(ξ_list...))

    m̃ = m(ns, nd, nθ)
    
    flattheory = reshape((interp_ξ)', (size(inv_S)[1],))
    Turing.@addlogprob! loglike_SH(d, flattheory, inv_S, ns, m̃)
    
    return nothing
end