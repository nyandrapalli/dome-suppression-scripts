using CSV
using DataFrames
using Ipopt
using JuMP
using LinearAlgebra
using Statistics
using TOML

const MOI = JuMP.MOI
const ROOT = normpath(joinpath(@__DIR__, ".."))
const CONTRACT_PATH = joinpath(ROOT, "fig7_geometry_contract.toml")
const MODE = get(ENV, "FIG7_V5_MODE", "audit") # audit | smoke
const RUN_TAG = get(ENV, "FIG7_V5_RUN_TAG", "unlabeled")
const OUT = get(ENV, "FIG7_V5_OUT", joinpath(@__DIR__, "outputs", "fig7_fluid_domain_v5_$(RUN_TAG)"))
const N_NODES = parse(Int, get(ENV, "FIG7_V5_N_NODES", "81"))
const IPOPT_MAXITER = parse(Int, get(ENV, "FIG7_V5_IPOPT_MAXITER", "2500"))
const IPOPT_TOL = parse(Float64, get(ENV, "FIG7_V5_IPOPT_TOL", "1e-8"))
const KAPPA_LO_KBT = parse(Float64, get(ENV, "FIG7_V5_KAPPA_LO_KBT", "80.0"))
const KAPPA_LD_KBT = parse(Float64, get(ENV, "FIG7_V5_KAPPA_LD_KBT", "25.0"))
const LINE_TENSION_PN = parse(Float64, get(ENV, "FIG7_V5_LINE_TENSION_PN", "0.5"))
const KA_LO_J_M2 = parse(Float64, get(ENV, "FIG7_V5_KA_LO_J_M2", "0.20"))
const KA_LD_J_M2 = parse(Float64, get(ENV, "FIG7_V5_KA_LD_J_M2", "0.20"))
const KDELTA_LO_J_M2 = parse(Float64, get(ENV, "FIG7_V5_KDELTA_LO_J_M2", "0.20"))
const KDELTA_LD_J_M2 = parse(Float64, get(ENV, "FIG7_V5_KDELTA_LD_J_M2", "0.20"))
const NEUTRAL_OFFSET_UM = parse(Float64, get(ENV, "FIG7_V5_NEUTRAL_OFFSET_UM", "0.002"))
const DELTA_KAPPA_G_KBT = parse(Float64, get(ENV, "FIG7_V5_DELTA_KAPPA_G_KBT", "0.0"))
const KBT_J = parse(Float64, get(ENV, "FIG7_V5_KBT_J", "4.114e-21"))
const CTB_LO_AREA_EXPANSION = parse(Float64, get(ENV, "FIG7_V5_CTB_LO_AREA_EXPANSION", "0.0458"))
const CTB_LO_MISMATCH = parse(Float64, get(ENV, "FIG7_V5_CTB_LO_MISMATCH", "0.005"))
const STREP_LD_AREA_EXPANSION = parse(Float64, get(ENV, "FIG7_V5_STREP10_LD_AREA_EXPANSION", "0.001"))
const STREP_LD_MISMATCH = parse(Float64, get(ENV, "FIG7_V5_STREP10_LD_MISMATCH", "0.0005"))
const PRE_SPHERE_TOL_UM = parse(Float64, get(ENV, "FIG7_V5_PRE_SPHERE_TOL_UM", "0.25"))
const HEIGHT_TOL_UM = parse(Float64, get(ENV, "FIG7_V5_HEIGHT_TOL_UM", "0.15"))

struct Contract
    radius::Float64
    area::Float64
    volume::Float64
    lo_area::Float64
    ld_area::Float64
    boundary_radius::Float64
    target_height::Float64
    target_cap_height::Float64
    target_cap_diameter::Float64
end

function load_contract()
    raw = TOML.parsefile(CONTRACT_PATH)
    pre = raw["pre_ctb_10nM_streptavidin"]
    post = raw["post_ctb_10nM_streptavidin"]
    cap = post["cap_measurement"]
    Contract(pre["vesicle_radius_um"], pre["total_area_um2"], pre["total_volume_um3"],
        pre["lo_area_um2"], pre["ld_area_um2"], pre["reference_lo_boundary_radius_um"],
        post["lo_height_increase_um"], cap["lo_cap_height_um"], cap["lo_cap_diameter_um"])
end

const C = load_contract()
const LINE_TENSION_KBT_UM = LINE_TENSION_PN * 1e-18 / KBT_J
const KA_LO = KA_LO_J_M2 * 1e-12 / KBT_J
const KA_LD = KA_LD_J_M2 * 1e-12 / KBT_J
const KDELTA_LO = KDELTA_LO_J_M2 * 1e-12 / KBT_J
const KDELTA_LD = KDELTA_LD_J_M2 * 1e-12 / KBT_J
const THETA_DB = asin(C.boundary_radius / C.radius)
const N_SEGMENTS = N_NODES - 1
const N_LO_SEGMENTS = clamp(round(Int, N_SEGMENTS * THETA_DB / pi), 4, N_SEGMENTS - 4)
const DB_NODE = N_LO_SEGMENTS + 1
const THETA0 = vcat(collect(range(0.0, THETA_DB, length=DB_NODE)),
                    collect(range(THETA_DB, pi, length=N_NODES - DB_NODE + 1))[2:end])

struct ProteinState
    label::String
    lo_area_expansion::Float64
    lo_mismatch::Float64
    ld_area_expansion::Float64
    ld_mismatch::Float64
end

function ensure_empty_output()
    isdir(OUT) && !isempty(readdir(OUT)) && error("Refusing to overwrite nonempty output directory: $(OUT)")
    mkpath(OUT)
end

function spherical_seed()
    rho = C.radius .* sin.(THETA0)
    z = C.radius .* cos.(THETA0)
    psi = copy(THETA0)
    ell = [hypot(rho[i + 1] - rho[i], z[i + 1] - z[i]) for i in 1:N_SEGMENTS]
    (; rho, z, psi, ell)
end

function domed_seed()
    seed = spherical_seed()
    rho, z = copy(seed.rho), copy(seed.z)
    for i in 1:DB_NODE
        u = THETA0[i] / THETA_DB
        blend = 1 - 3u^2 + 2u^3
        z[i] += C.target_height * blend
    end
    shift = (z[1] + z[end]) / 2
    z .-= shift
    ell = [hypot(rho[i + 1] - rho[i], z[i + 1] - z[i]) for i in 1:N_SEGMENTS]
    psi = zeros(N_NODES)
    psi[1], psi[end] = 0.0, pi
    for i in 2:(N_NODES - 1)
        dr = rho[i + 1] - rho[i - 1]
        dz = z[i + 1] - z[i - 1]
        psi[i] = clamp(atan(-dz, dr), 0.0, pi)
    end
    (; rho, z, psi, ell)
end

function set_start_shape!(rho, z, psi, ell, seed)
    for i in 1:N_NODES
        set_start_value(rho[i], seed.rho[i])
        set_start_value(z[i], seed.z[i])
        set_start_value(psi[i], seed.psi[i])
    end
    for i in 1:N_SEGMENTS
        set_start_value(ell[i], seed.ell[i])
    end
end

function reference_leaflet_difference(area)
    4 * NEUTRAL_OFFSET_UM * area / C.radius
end

function solve_state(state::ProteinState, seed; branch::String)
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "tol", IPOPT_TOL)
    set_optimizer_attribute(model, "constr_viol_tol", IPOPT_TOL)
    set_optimizer_attribute(model, "acceptable_tol", 10 * IPOPT_TOL)
    set_optimizer_attribute(model, "max_iter", IPOPT_MAXITER)
    set_optimizer_attribute(model, "print_level", 0)

    @variable(model, 0.0 <= rho[1:N_NODES] <= 20.0)
    @variable(model, -20.0 <= z[1:N_NODES] <= 20.0)
    @variable(model, 0.0 <= psi[1:N_NODES] <= pi)
    @variable(model, 1e-4 <= ell[1:N_SEGMENTS] <= 5.0)
    set_start_shape!(rho, z, psi, ell, seed)

    fix(rho[1], 0.0; force=true)
    fix(rho[end], 0.0; force=true)
    fix(psi[1], 0.0; force=true)
    fix(psi[end], pi; force=true)
    @constraint(model, z[1] + z[end] == 0.0)

    for i in 1:N_SEGMENTS
        @NLconstraint(model, rho[i + 1] - rho[i] == ell[i] * cos((psi[i] + psi[i + 1]) / 2))
        @NLconstraint(model, z[i + 1] - z[i] == -ell[i] * sin((psi[i] + psi[i + 1]) / 2))
    end

    @NLexpression(model, area_seg[i=1:N_SEGMENTS], pi * (rho[i] + rho[i + 1]) * ell[i])
    @NLexpression(model, ks[i=1:N_SEGMENTS], (psi[i + 1] - psi[i]) / ell[i])
    @NLexpression(model, kphi[i=1:N_SEGMENTS],
        sin((psi[i] + psi[i + 1]) / 2) / ((rho[i] + rho[i + 1]) / 2 + 1e-10))
    @NLexpression(model, two_m[i=1:N_SEGMENTS], ks[i] + kphi[i])
    @NLexpression(model, total_volume,
        sum(pi * (rho[i]^2 + rho[i + 1]^2) * (z[i] - z[i + 1]) / 2 for i in 1:N_SEGMENTS))
    @NLconstraint(model, total_volume == C.volume)

    @NLexpression(model, area_lo, sum(area_seg[i] for i in 1:N_LO_SEGMENTS))
    @NLexpression(model, area_ld, sum(area_seg[i] for i in (N_LO_SEGMENTS + 1):N_SEGMENTS))
    @NLexpression(model, delta_area_lo,
        2 * NEUTRAL_OFFSET_UM * sum(two_m[i] * area_seg[i] for i in 1:N_LO_SEGMENTS))
    @NLexpression(model, delta_area_ld,
        2 * NEUTRAL_OFFSET_UM * sum(two_m[i] * area_seg[i] for i in (N_LO_SEGMENTS + 1):N_SEGMENTS))

    area_pref_lo = C.lo_area * (1 + state.lo_area_expansion)
    area_pref_ld = C.ld_area * (1 + state.ld_area_expansion)
    delta_pref_lo = reference_leaflet_difference(C.lo_area) + state.lo_mismatch * C.lo_area
    delta_pref_ld = reference_leaflet_difference(C.ld_area) - state.ld_mismatch * C.ld_area

    @NLexpression(model, bend_lo,
        0.5 * KAPPA_LO_KBT * sum(two_m[i]^2 * area_seg[i] for i in 1:N_LO_SEGMENTS))
    @NLexpression(model, bend_ld,
        0.5 * KAPPA_LD_KBT * sum(two_m[i]^2 * area_seg[i] for i in (N_LO_SEGMENTS + 1):N_SEGMENTS))
    @NLexpression(model, area_energy_lo, KA_LO * (area_lo - area_pref_lo)^2 / (2 * C.lo_area))
    @NLexpression(model, area_energy_ld, KA_LD * (area_ld - area_pref_ld)^2 / (2 * C.ld_area))
    @NLexpression(model, mismatch_energy_lo,
        KDELTA_LO * (delta_area_lo - delta_pref_lo)^2 / (2 * C.lo_area))
    @NLexpression(model, mismatch_energy_ld,
        KDELTA_LD * (delta_area_ld - delta_pref_ld)^2 / (2 * C.ld_area))
    @NLexpression(model, line_energy, 2 * pi * LINE_TENSION_KBT_UM * rho[DB_NODE])
    @NLexpression(model, gaussian_energy, 2 * pi * DELTA_KAPPA_G_KBT * cos(psi[DB_NODE]))
    @NLobjective(model, Min, bend_lo + bend_ld + area_energy_lo + area_energy_ld +
        mismatch_energy_lo + mismatch_energy_ld + line_energy + gaussian_energy)

    optimize!(model)
    term = termination_status(model)
    primal = primal_status(model)
    has_solution = primal in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    has_solution || return (status=string(term), primal_status=string(primal), branch=branch,
        accepted=false, model=model, shape=nothing, metrics=nothing)

    shape = (rho=value.(rho), z=value.(z), psi=value.(psi), ell=value.(ell))
    metrics = compute_metrics(shape, state)
    accepted = term in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL) &&
               metrics.max_geometry_error < 1e-6 && abs(metrics.volume_error) < 1e-8
    (status=string(term), primal_status=string(primal), branch=branch, accepted=accepted,
        objective=objective_value(model), solve_time=solve_time(model), model=model, shape=shape, metrics=metrics)
end

function compute_metrics(shape, state)
    rho, z, psi, ell = shape.rho, shape.z, shape.psi, shape.ell
    area = [pi * (rho[i] + rho[i + 1]) * ell[i] for i in 1:N_SEGMENTS]
    ks = [(psi[i + 1] - psi[i]) / ell[i] for i in 1:N_SEGMENTS]
    kphi = [sin((psi[i] + psi[i + 1]) / 2) / ((rho[i] + rho[i + 1]) / 2 + 1e-10) for i in 1:N_SEGMENTS]
    two_m = ks .+ kphi
    volume = sum(pi * (rho[i]^2 + rho[i + 1]^2) * (z[i] - z[i + 1]) / 2 for i in 1:N_SEGMENTS)
    alo, ald = sum(area[1:N_LO_SEGMENTS]), sum(area[(N_LO_SEGMENTS + 1):end])
    dalo = 2 * NEUTRAL_OFFSET_UM * sum(two_m[1:N_LO_SEGMENTS] .* area[1:N_LO_SEGMENTS])
    dald = 2 * NEUTRAL_OFFSET_UM * sum(two_m[(N_LO_SEGMENTS + 1):end] .* area[(N_LO_SEGMENTS + 1):end])
    alo_pref = C.lo_area * (1 + state.lo_area_expansion)
    ald_pref = C.ld_area * (1 + state.ld_area_expansion)
    dalo_pref = reference_leaflet_difference(C.lo_area) + state.lo_mismatch * C.lo_area
    dald_pref = reference_leaflet_difference(C.ld_area) - state.ld_mismatch * C.ld_area
    bend_lo = 0.5 * KAPPA_LO_KBT * sum(two_m[1:N_LO_SEGMENTS].^2 .* area[1:N_LO_SEGMENTS])
    bend_ld = 0.5 * KAPPA_LD_KBT * sum(two_m[(N_LO_SEGMENTS + 1):end].^2 .* area[(N_LO_SEGMENTS + 1):end])
    ealo = KA_LO * (alo - alo_pref)^2 / (2 * C.lo_area)
    eald = KA_LD * (ald - ald_pref)^2 / (2 * C.ld_area)
    edlo = KDELTA_LO * (dalo - dalo_pref)^2 / (2 * C.lo_area)
    edld = KDELTA_LD * (dald - dald_pref)^2 / (2 * C.ld_area)
    eline = 2pi * LINE_TENSION_KBT_UM * rho[DB_NODE]
    egauss = 2pi * DELTA_KAPPA_G_KBT * cos(psi[DB_NODE])
    geom_errors = Float64[]
    for i in 1:N_SEGMENTS
        push!(geom_errors, rho[i + 1] - rho[i] - ell[i] * cos((psi[i] + psi[i + 1]) / 2))
        push!(geom_errors, z[i + 1] - z[i] + ell[i] * sin((psi[i] + psi[i + 1]) / 2))
    end
    radial = hypot.(rho, z)
    (; volume, volume_error=(volume - C.volume) / C.volume, area_lo=alo, area_ld=ald,
       total_area=alo + ald, delta_area_lo=dalo, delta_area_ld=dald,
       area_pref_lo=alo_pref, area_pref_ld=ald_pref, delta_pref_lo=dalo_pref, delta_pref_ld=dald_pref,
       bend_lo, bend_ld, area_energy_lo=ealo, area_energy_ld=eald,
       mismatch_energy_lo=edlo, mismatch_energy_ld=edld, line_energy=eline, gaussian_energy=egauss,
       total_energy=bend_lo + bend_ld + ealo + eald + edlo + edld + eline + egauss,
       apex_displacement=z[1] - C.radius, cap_height=z[1] - z[DB_NODE],
       cap_diameter=2rho[DB_NODE], boundary_radius=rho[DB_NODE],
       radial_cv=std(radial) / mean(radial), max_geometry_error=maximum(abs.(geom_errors)),
       segment_area=area, ks, kphi, two_m)
end

function write_solution(prefix, result)
    summary = DataFrame([merge((branch=result.branch, termination_status=result.status,
        primal_status=result.primal_status, equilibrium_accepted=result.accepted,
        objective_kbt=get(result, :objective, NaN), solve_time_s=get(result, :solve_time, NaN)),
        result.metrics === nothing ? NamedTuple() : NamedTuple{propertynames(result.metrics)[1:end-4]}(
            Tuple(getproperty(result.metrics, p) for p in propertynames(result.metrics)[1:end-4])))])
    CSV.write(joinpath(OUT, "$(prefix)_summary.csv"), summary)
    result.shape === nothing && return
    s, m = result.shape, result.metrics
    CSV.write(joinpath(OUT, "$(prefix)_nodes.csv"), DataFrame(node=1:N_NODES, rho_um=s.rho,
        z_um=s.z, psi_rad=s.psi, domain=[i <= DB_NODE ? "Lo" : "Ld" for i in 1:N_NODES]))
    CSV.write(joinpath(OUT, "$(prefix)_segments.csv"), DataFrame(segment=1:N_SEGMENTS,
        domain=[i <= N_LO_SEGMENTS ? "Lo" : "Ld" for i in 1:N_SEGMENTS], ell_um=s.ell,
        area_um2=m.segment_area, meridional_curvature_um_inv=m.ks,
        azimuthal_curvature_um_inv=m.kphi, two_mean_curvature_um_inv=m.two_m))
end

function state_rows(results)
    rows = NamedTuple[]
    for (label, result) in results
        m = result.metrics
        push!(rows, m === nothing ? (label=label, branch=result.branch, status=result.status,
            accepted=false) : (label=label, branch=result.branch, status=result.status,
            accepted=result.accepted, objective_kbt=result.objective, volume_error=m.volume_error,
            geometry_error=m.max_geometry_error, apex_displacement_um=m.apex_displacement,
            cap_height_um=m.cap_height, cap_diameter_um=m.cap_diameter,
            lo_area_um2=m.area_lo, ld_area_um2=m.area_ld, total_area_um2=m.total_area,
            lo_area_change_percent=100 * (m.area_lo / C.lo_area - 1),
            bending_lo_kbt=m.bend_lo, bending_ld_kbt=m.bend_ld,
            area_energy_lo_kbt=m.area_energy_lo, area_energy_ld_kbt=m.area_energy_ld,
            mismatch_energy_lo_kbt=m.mismatch_energy_lo, mismatch_energy_ld_kbt=m.mismatch_energy_ld,
            line_energy_kbt=m.line_energy, total_energy_kbt=m.total_energy))
    end
    DataFrame(rows)
end

function write_audit()
    CSV.write(joinpath(OUT, "condition_audit.csv"), DataFrame(
        field=["shape", "constraint solver", "volume", "Lo/Ld area", "leaflet mismatch",
            "CTB Lo area expansion", "CTB Lo mismatch", "Strep Ld area expansion",
            "Strep Ld mismatch", "protein C0", "direct cross-domain coupling"],
        value=["$(N_NODES)-node axisymmetric meridian", "Ipopt nonlinear equality constraints",
            "fixed at $(C.volume) um3", "fluid domain-level preferred-area energy",
            "domain-level outer-minus-inner preferred-area energy", CTB_LO_AREA_EXPANSION,
            CTB_LO_MISMATCH, STREP_LD_AREA_EXPANSION, STREP_LD_MISMATCH, "inactive", "inactive"],
    ))
    CSV.write(joinpath(OUT, "run_config.csv"), DataFrame(parameter=["mode", "n_nodes", "ipopt_tol",
        "ipopt_maxiter", "target_apex_displacement_um", "target_cap_height_um",
        "target_cap_diameter_um", "measured_lo_area_change_percent"], value=[MODE, N_NODES,
        IPOPT_TOL, IPOPT_MAXITER, C.target_height, C.target_cap_height, C.target_cap_diameter,
        100 * (pi * ((C.target_cap_diameter / 2)^2 + C.target_cap_height^2) / C.lo_area - 1)]))
end

function main()
    MODE in ("audit", "smoke") || error("FIG7_V5_MODE must be audit or smoke.")
    N_NODES >= 21 || error("FIG7_V5_N_NODES must be at least 21.")
    ensure_empty_output(); write_audit()
    if MODE == "audit"
        write(joinpath(OUT, "status.txt"), "audit_only: no minimization was run.\n")
        return
    end
    baseline_state = ProteinState("baseline", 0.0, 0.0, 0.0, 0.0)
    pre_state = ProteinState("pre_ctb_10nM_strep", 0.0, 0.0,
        STREP_LD_AREA_EXPANSION, STREP_LD_MISMATCH)
    post_state = ProteinState("post_ctb_10nM_strep", CTB_LO_AREA_EXPANSION,
        CTB_LO_MISMATCH, STREP_LD_AREA_EXPANSION, STREP_LD_MISMATCH)

    baseline = solve_state(baseline_state, spherical_seed(); branch="sphere")
    baseline.shape === nothing && error("Baseline solve failed without a feasible point.")
    pre = solve_state(pre_state, baseline.shape; branch="continuation")
    pre.shape === nothing && error("Pre-CTB solve failed without a feasible point.")
    post_sphere = solve_state(post_state, pre.shape; branch="near_sphere")
    post_dome = solve_state(post_state, domed_seed(); branch="fig7_domed_seed")
    results = [("baseline", baseline), ("pre", pre), ("post_sphere", post_sphere), ("post_dome", post_dome)]
    for (label, result) in results
        write_solution(label, result)
    end
    summary = state_rows(results)
    CSV.write(joinpath(OUT, "branch_summary.csv"), summary)
    valid_post = [r for r in (post_sphere, post_dome) if r.accepted]
    if isempty(valid_post)
        write(joinpath(OUT, "status.txt"), "no_converged_post_ctb_branch\n")
        return
    end
    best = valid_post[argmin([r.objective for r in valid_post])]
    mpre, mpost = pre.metrics, best.metrics
    pre_error = max(abs(mpre.apex_displacement), C.radius * mpre.radial_cv)
    height_error = abs(mpost.apex_displacement - C.target_height)
    calibration_ok = baseline.accepted && pre.accepted && pre_error <= PRE_SPHERE_TOL_UM && height_error <= HEIGHT_TOL_UM
    CSV.write(joinpath(OUT, "selected_branch.csv"), DataFrame(branch=[best.branch],
        objective_kbt=[best.objective], pre_sphere_error_um=[pre_error],
        post_height_error_um=[height_error], post_apex_displacement_um=[mpost.apex_displacement],
        post_cap_height_um=[mpost.cap_height], post_cap_diameter_um=[mpost.cap_diameter],
        post_lo_area_change_percent=[100 * (mpost.area_lo / C.lo_area - 1)],
        calibration_accepted=[calibration_ok]))
    write(joinpath(OUT, "status.txt"), calibration_ok ? "feasible_calibration_found\n" : "converged_branches_but_no_calibration\n")
end

get(ENV, "FIG7_V5_LIBRARY_MODE", "false") == "true" || main()
