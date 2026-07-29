# V7.7.0: combine the corrected free-interface V7.6 mechanics with a local
# inner-Ld leaflet-area source. Streptavidin concentration changes the
# preferred inner-Ld area, while the validated 10 nM spontaneous curvatures
# remain fixed. No direct Lo suppression or fitted cross-domain force is used.
#
# The relative magnitudes come from the triangulated four-control study. A
# single transfer factor preserves the V7.3 net Ld curvature at 10 nM, while
# the CTB component is adjusted to preserve the V7.3 post-CTB net Lo curvature.

const V5_SOURCE = normpath(joinpath(@__DIR__, "..", "v5", "run_fig7_fluid_domain_area_mismatch_v5.jl"))
const V73_LINE_TENSION_PN = parse(Float64, get(ENV, "FIG7_V770_LINE_TENSION_PN", "0.3"))
const V73_GAUSSIAN_RATIO = parse(Float64, get(ENV, "FIG7_V770_GAUSSIAN_RATIO", "0.9"))
const V73_KAPPA_LO_KBT = parse(Float64, get(ENV, "FIG7_V5_KAPPA_LO_KBT", "80.0"))
const V73_KAPPA_LD_KBT = parse(Float64, get(ENV, "FIG7_V5_KAPPA_LD_KBT", "25.0"))
const V73_DELTA_KAPPA_G_KBT =
    V73_GAUSSIAN_RATIO * (V73_KAPPA_LO_KBT - V73_KAPPA_LD_KBT)
ENV["FIG7_V5_LINE_TENSION_PN"] = string(V73_LINE_TENSION_PN)
ENV["FIG7_V5_DELTA_KAPPA_G_KBT"] = string(V73_DELTA_KAPPA_G_KBT)
ENV["FIG7_V5_LIBRARY_MODE"] = "true"
include(V5_SOURCE)
delete!(ENV, "FIG7_V5_LIBRARY_MODE")
delete!(ENV, "FIG7_V5_LINE_TENSION_PN")
delete!(ENV, "FIG7_V5_DELTA_KAPPA_G_KBT")

const V51_MODE = get(ENV, "FIG7_V770_MODE", "audit") # audit | smoke
const V51_TAG = get(ENV, "FIG7_V770_RUN_TAG", "unlabeled")
const V51_OUT = get(ENV, "FIG7_V770_OUT", joinpath(@__DIR__, "outputs", "fig7_local_leaflet_interface_v7_7_0_$(V51_TAG)"))
const V51_IPOPT_MAXITER = parse(Int, get(ENV, "FIG7_V770_IPOPT_MAXITER", "7000"))
const V51_IPOPT_TOL = parse(Float64, get(ENV, "FIG7_V770_IPOPT_TOL", "1e-8"))
const V51_IPOPT_PRINT_LEVEL = parse(Int, get(ENV, "FIG7_V770_IPOPT_PRINT_LEVEL", "0"))
const V760_GEOMETRY_TOL = parse(Float64,
    get(ENV, "FIG7_V770_GEOMETRY_TOL", string(10 * V51_IPOPT_TOL)))
const V760_AXIS_CONTACT_TOL_UM = parse(Float64,
    get(ENV, "FIG7_V770_AXIS_CONTACT_TOL_UM", "1e-5"))
const MAX_SEGMENT_RATIO = parse(Float64, get(ENV, "FIG7_V770_MAX_SEGMENT_RATIO", "2.0"))
const EPS_CTB_OUTER_LO = parse(Float64, get(ENV, "FIG7_V770_EPS_CTB_OUTER_LO", "0.0916"))
const PRIOR_C0_LD_LIPID_MAG_UM_INV =
    parse(Float64, get(ENV, "FIG7_V770_PRIOR_C0_LD_LIPID_MAG_UM_INV", "0.02"))
const PRIOR_C0_LD_STREP10_MAG_UM_INV =
    parse(Float64, get(ENV, "FIG7_V770_PRIOR_C0_LD_STREP10_MAG_UM_INV", "0.0625"))
const PRIOR_C0_LO_LIPID_MAG_UM_INV =
    parse(Float64, get(ENV, "FIG7_V770_PRIOR_C0_LO_LIPID_MAG_UM_INV", "0.035"))
const TARGET_C0_LD_COMBINED10_UM_INV =
    parse(Float64, get(ENV, "FIG7_V770_TARGET_C0_LD_COMBINED10_UM_INV", "-0.015"))
const TARGET_C0_LO_POST_UM_INV =
    parse(Float64, get(ENV, "FIG7_V770_TARGET_C0_LO_POST_UM_INV", string(2 / 6.0)))
const V7510_TRANSFER_SCALE =
    abs(TARGET_C0_LD_COMBINED10_UM_INV) /
    (PRIOR_C0_LD_LIPID_MAG_UM_INV + PRIOR_C0_LD_STREP10_MAG_UM_INV)
const C0_LD_LIPID_UM_INV =
    -V7510_TRANSFER_SCALE * PRIOR_C0_LD_LIPID_MAG_UM_INV
const C0_LO_LIPID_UM_INV =
    V7510_TRANSFER_SCALE * PRIOR_C0_LO_LIPID_MAG_UM_INV
const C0_LD_STREP10_UM_INV =
    -V7510_TRANSFER_SCALE * PRIOR_C0_LD_STREP10_MAG_UM_INV
const C0_LD_COMBINED10_UM_INV = C0_LD_LIPID_UM_INV + C0_LD_STREP10_UM_INV
const C0_LO_CTB_UM_INV = TARGET_C0_LO_POST_UM_INV - C0_LO_LIPID_UM_INV
const C0_LO_POST_UM_INV = C0_LO_LIPID_UM_INV + C0_LO_CTB_UM_INV
const V7510_STREP_SATURATION_NM =
    parse(Float64, get(ENV, "FIG7_V770_STREP_SATURATION_NM", "80.0"))
const INNER_LD_BIOTIN_MOLE_FRACTION =
    parse(Float64, get(ENV, "FIG7_V770_INNER_LD_BIOTIN_MOLE_FRACTION", "0.002"))
const BOUND_SITE_AREA_COUPLING =
    parse(Float64, get(ENV, "FIG7_V770_BOUND_SITE_AREA_COUPLING", "1.0"))
const EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM =
    parse(Float64, get(ENV, "FIG7_V770_DOME_INCREASE_TARGET_UM", "2.95"))
const CORRECTED_DOME_ACCEPTANCE_TOL_UM =
    parse(Float64, get(ENV, "FIG7_V770_DOME_INCREASE_TOL_UM", "0.50"))
const V770_CONCENTRATIONS_NM = parse.(Float64,
    filter(x -> !isempty(x), split(strip(get(ENV, "FIG7_V770_CONCENTRATIONS_NM", "10,40,80")),
        r"[,;\s]+")))
const V770_MONOTONIC_TOL_UM =
    parse(Float64, get(ENV, "FIG7_V770_MONOTONIC_TOL_UM", "0.02"))
const V770_INTERFACE_MOMENT_REL_TOL =
    parse(Float64, get(ENV, "FIG7_V770_INTERFACE_MOMENT_REL_TOL", "0.10"))
const V770_INTERFACE_TRANSVERSE_REL_TOL =
    parse(Float64, get(ENV, "FIG7_V770_INTERFACE_TRANSVERSE_REL_TOL", "0.50"))
const TARGET_PROFILE_PATH = "inactive: obsolete V6 20.3413-um contour excluded"

struct V7State
    label::String
    lo_area_expansion::Float64
    lo_mismatch::Float64
    ld_area_expansion::Float64
    ld_mismatch::Float64
    c0_lo_um_inv::Float64
    c0_ld_um_inv::Float64
end

function ensure_empty_v51_output()
    isdir(V51_OUT) && !isempty(readdir(V51_OUT)) && error("Refusing to overwrite nonempty output directory: $(V51_OUT)")
    mkpath(V51_OUT)
end

function segment_lengths_v51(ell_lo, ell_ld)
    [i <= N_LO_SEGMENTS ? ell_lo : ell_ld for i in 1:N_SEGMENTS]
end

function raw_metrics_v51(shape)
    rho, z, psi = shape.rho, shape.z, shape.psi
    ell = segment_lengths_v51(shape.ell_lo, shape.ell_ld)
    area = [pi * (rho[i] + rho[i + 1]) * ell[i] for i in 1:N_SEGMENTS]
    ks = [(psi[i + 1] - psi[i]) / ell[i] for i in 1:N_SEGMENTS]
    kphi = [sin((psi[i] + psi[i + 1]) / 2) / ((rho[i] + rho[i + 1]) / 2 + 1e-10) for i in 1:N_SEGMENTS]
    two_m = ks .+ kphi
    volume = sum(pi * (z[i] - z[i + 1]) *
        (rho[i]^2 + rho[i] * rho[i + 1] + rho[i + 1]^2) / 3 for i in 1:N_SEGMENTS)
    (; ell, area, ks, kphi, two_m, volume)
end

function reference_v51()
    seed = spherical_seed()
    shape = (; rho=seed.rho, z=seed.z, psi=seed.psi,
        ell_lo=mean(seed.ell[1:N_LO_SEGMENTS]),
        ell_ld=mean(seed.ell[(N_LO_SEGMENTS + 1):end]))
    raw = raw_metrics_v51(shape)
    area_lo = sum(raw.area[1:N_LO_SEGMENTS])
    area_ld = sum(raw.area[(N_LO_SEGMENTS + 1):end])
    delta_lo = 2 * NEUTRAL_OFFSET_UM * sum(raw.two_m[1:N_LO_SEGMENTS] .* raw.area[1:N_LO_SEGMENTS])
    delta_ld = 2 * NEUTRAL_OFFSET_UM * sum(raw.two_m[(N_LO_SEGMENTS + 1):end] .* raw.area[(N_LO_SEGMENTS + 1):end])
    (; shape, volume=raw.volume, area_lo, area_ld, delta_lo, delta_ld)
end

const REF51 = reference_v51()
const CONTRACT_RAW_V71 = TOML.parsefile(CONTRACT_PATH)
const PRE_CAP_HEIGHT_V71 = Float64(
    CONTRACT_RAW_V71["pre_ctb_10nM_streptavidin"]["reference_lo_cap_height_um"])
const POST_CAP_HEIGHT_V71 = Float64(
    CONTRACT_RAW_V71["post_ctb_10nM_streptavidin"]["cap_measurement"]["lo_cap_height_um"])
const LO_DELTA_AREA_SHIFT_V71 = 8pi * NEUTRAL_OFFSET_UM *
    (POST_CAP_HEIGHT_V71 - PRE_CAP_HEIGHT_V71)
const LO_GEOMETRIC_MISMATCH_V71 = LO_DELTA_AREA_SHIFT_V71 / REF51.area_lo
const LD_GEOMETRIC_MISMATCH_V71 = 0.0

function start_shape_v51(seed)
    (; rho=seed.rho, z=seed.z, psi=seed.psi,
       ell_lo=mean(seed.ell[1:N_LO_SEGMENTS]),
       ell_ld=mean(seed.ell[(N_LO_SEGMENTS + 1):end]))
end

function profile_self_intersections_v760(rho, z; tol=V760_GEOMETRY_TOL)
    orient(ax, ay, bx, by, cx, cy) = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    intersections = 0
    for i in 1:N_SEGMENTS
        i + 2 > N_SEGMENTS && continue
        for j in (i + 2):N_SEGMENTS
            maximum((rho[i], rho[i + 1])) + tol < minimum((rho[j], rho[j + 1])) && continue
            maximum((rho[j], rho[j + 1])) + tol < minimum((rho[i], rho[i + 1])) && continue
            maximum((z[i], z[i + 1])) + tol < minimum((z[j], z[j + 1])) && continue
            maximum((z[j], z[j + 1])) + tol < minimum((z[i], z[i + 1])) && continue
            o1 = orient(rho[i], z[i], rho[i + 1], z[i + 1], rho[j], z[j])
            o2 = orient(rho[i], z[i], rho[i + 1], z[i + 1], rho[j + 1], z[j + 1])
            o3 = orient(rho[j], z[j], rho[j + 1], z[j + 1], rho[i], z[i])
            o4 = orient(rho[j], z[j], rho[j + 1], z[j + 1], rho[i + 1], z[i + 1])
            if o1 * o2 < -tol^2 && o3 * o4 < -tol^2
                intersections += 1
            end
        end
    end
    intersections
end

function interface_balance_v760(shape, raw, state, mismatch_stress_lo, mismatch_stress_ld)
    N_LO_SEGMENTS >= 2 || error("At least two Lo segments are required.")
    N_SEGMENTS - N_LO_SEGMENTS >= 2 || error("At least two Ld segments are required.")
    lo_last = N_LO_SEGMENTS
    ld_first = N_LO_SEGMENTS + 1
    two_h_lo_b = 1.5 * raw.two_m[lo_last] - 0.5 * raw.two_m[lo_last - 1]
    two_h_ld_b = 1.5 * raw.two_m[ld_first] - 0.5 * raw.two_m[ld_first + 1]
    kphi_b = sin(shape.psi[DB_NODE]) / shape.rho[DB_NODE]

    # The nonlocal area-difference term adds 2*d*q to the domain bending moment.
    moment_lo = KAPPA_LO_KBT * (two_h_lo_b - state.c0_lo_um_inv) +
        2 * NEUTRAL_OFFSET_UM * mismatch_stress_lo
    moment_ld = KAPPA_LD_KBT * (two_h_ld_b - state.c0_ld_um_inv) +
        2 * NEUTRAL_OFFSET_UM * mismatch_stress_ld
    gaussian_moment = DELTA_KAPPA_G_KBT * kphi_b
    # For the north-to-south profile and Delta kappa_G = kappa_G,Ld-kappa_G,Lo,
    # stationarity of 2*pi*Delta kappa_G*cos(psi_b) gives
    # Delta kappa_G*k_phi = m_Lo-m_Ld.
    moment_residual = gaussian_moment - (moment_lo - moment_ld)
    moment_scale = max(abs(gaussian_moment), abs(moment_lo - moment_ld), 1e-12)

    d_two_h_lo_ds = (raw.two_m[lo_last] - raw.two_m[lo_last - 1]) / raw.ell[lo_last]
    d_two_h_ld_ds = (raw.two_m[ld_first + 1] - raw.two_m[ld_first]) / raw.ell[ld_first]
    transverse_bending = -KAPPA_LO_KBT * d_two_h_lo_ds +
        KAPPA_LD_KBT * d_two_h_ld_ds
    transverse_line = LINE_TENSION_KBT_UM * kphi_b
    transverse_residual = transverse_line - transverse_bending
    transverse_scale = max(abs(transverse_line), abs(transverse_bending), 1e-12)

    (; phase_mapping="north-to-south profile: Lo then Ld; delta_kappa_g=kappaG_Ld-kappaG_Lo",
       two_h_lo_boundary=two_h_lo_b, two_h_ld_boundary=two_h_ld_b,
       kphi_boundary=kphi_b, moment_lo, moment_ld, gaussian_moment,
       moment_residual, moment_relative_residual=abs(moment_residual) / moment_scale,
       d_two_h_lo_ds, d_two_h_ld_ds, transverse_line, transverse_bending,
       transverse_residual,
       transverse_relative_residual=abs(transverse_residual) / transverse_scale)
end

function compute_metrics_v51(shape, state)
    raw = raw_metrics_v51(shape)
    area_lo = sum(raw.area[1:N_LO_SEGMENTS])
    area_ld = sum(raw.area[(N_LO_SEGMENTS + 1):end])
    delta_lo = 2 * NEUTRAL_OFFSET_UM * sum(raw.two_m[1:N_LO_SEGMENTS] .* raw.area[1:N_LO_SEGMENTS])
    delta_ld = 2 * NEUTRAL_OFFSET_UM * sum(raw.two_m[(N_LO_SEGMENTS + 1):end] .* raw.area[(N_LO_SEGMENTS + 1):end])
    area_pref_lo = REF51.area_lo * (1 + state.lo_area_expansion)
    area_pref_ld = REF51.area_ld * (1 + state.ld_area_expansion)
    delta_pref_lo = REF51.delta_lo + state.lo_mismatch * REF51.area_lo
    delta_pref_ld = REF51.delta_ld + state.ld_mismatch * REF51.area_ld
    bend_lo = 0.5 * KAPPA_LO_KBT * sum(
        (raw.two_m[1:N_LO_SEGMENTS] .- state.c0_lo_um_inv).^2 .*
        raw.area[1:N_LO_SEGMENTS])
    bend_ld = 0.5 * KAPPA_LD_KBT * sum(
        (raw.two_m[(N_LO_SEGMENTS + 1):end] .- state.c0_ld_um_inv).^2 .*
        raw.area[(N_LO_SEGMENTS + 1):end])
    ealo = KA_LO * (area_lo - area_pref_lo)^2 / (2 * REF51.area_lo)
    eald = KA_LD * (area_ld - area_pref_ld)^2 / (2 * REF51.area_ld)
    edlo = KDELTA_LO * (delta_lo - delta_pref_lo)^2 / (2 * REF51.area_lo)
    edld = KDELTA_LD * (delta_ld - delta_pref_ld)^2 / (2 * REF51.area_ld)
    mismatch_stress_lo = KDELTA_LO * (delta_lo - delta_pref_lo) / REF51.area_lo
    mismatch_stress_ld = KDELTA_LD * (delta_ld - delta_pref_ld) / REF51.area_ld
    eline = 2pi * LINE_TENSION_KBT_UM * shape.rho[DB_NODE]
    egauss = 2pi * DELTA_KAPPA_G_KBT * cos(shape.psi[DB_NODE])
    geom = Float64[]
    for i in 1:N_SEGMENTS
        push!(geom, shape.rho[i + 1] - shape.rho[i] - raw.ell[i] * cos((shape.psi[i] + shape.psi[i + 1]) / 2))
        push!(geom, shape.z[i + 1] - shape.z[i] + raw.ell[i] * sin((shape.psi[i] + shape.psi[i + 1]) / 2))
    end
    dpsi = diff(shape.psi)
    negative_meridional_segments = count(<(-V760_GEOMETRY_TOL), dpsi)
    axial_order_violations = count(>(V760_GEOMETRY_TOL), diff(shape.z))
    interior_axis_contacts = count(<=(V760_AXIS_CONTACT_TOL_UM), shape.rho[2:(end - 1)])
    profile_self_intersections = profile_self_intersections_v760(shape.rho, shape.z)
    interface = interface_balance_v760(
        shape, raw, state, mismatch_stress_lo, mismatch_stress_ld)
    ratio = maximum(raw.ell) / minimum(raw.ell)
    radial = hypot.(shape.rho, shape.z)
    total_area = area_lo + area_ld
    sphericity = pi^(1 / 3) * (6 * abs(raw.volume))^(2 / 3) / total_area
    (; volume=raw.volume, volume_error=(raw.volume - REF51.volume) / REF51.volume,
       area_lo, area_ld, total_area, sphericity, delta_area_lo=delta_lo,
       delta_area_ld=delta_ld, area_pref_lo, area_pref_ld, delta_pref_lo, delta_pref_ld,
       c0_lo_um_inv=state.c0_lo_um_inv, c0_ld_um_inv=state.c0_ld_um_inv,
       bend_lo, bend_ld, area_energy_lo=ealo, area_energy_ld=eald,
       mismatch_energy_lo=edlo, mismatch_energy_ld=edld, line_energy=eline,
       gaussian_energy=egauss, total_energy=bend_lo + bend_ld + ealo + eald + edlo + edld + eline + egauss,
       apex_displacement=shape.z[1] - C.radius, cap_height=shape.z[1] - shape.z[DB_NODE],
       cap_diameter=2shape.rho[DB_NODE], boundary_radius=shape.rho[DB_NODE],
       radial_cv=std(radial) / mean(radial), max_geometry_error=maximum(abs.(geom)),
       negative_meridional_segments, min_dpsi=minimum(dpsi),
       axial_order_violations, interior_axis_contacts, profile_self_intersections,
       mismatch_stress_lo, mismatch_stress_ld, interface,
       segment_ratio=ratio, min_segment=minimum(raw.ell),
       max_segment=maximum(raw.ell), segment_area=raw.area, ks=raw.ks, kphi=raw.kphi,
       two_m=raw.two_m, segment_length=raw.ell)
end

function normalized_profile_v6(rho, z)
    zspan = maximum(z) - minimum(z)
    rscale = maximum(rho)
    zspan > 0 || error("Profile has zero axial extent.")
    rscale > 0 || error("Profile has zero radial extent.")
    zn = (z .- minimum(z)) ./ zspan
    rn = rho ./ rscale
    order = sortperm(zn)
    zn[order], rn[order]
end

function interpolate_linear_v6(x, y, query)
    query <= x[1] && return y[1]
    query >= x[end] && return y[end]
    upper = searchsortedfirst(x, query)
    lower = upper - 1
    weight = (query - x[lower]) / (x[upper] - x[lower])
    (1 - weight) * y[lower] + weight * y[upper]
end

contour_shape_rms_v6(shape) = NaN

function ld_sphere_fit_v6(shape)
    indices = (DB_NODE + 2):N_NODES
    rho = shape.rho[indices]
    z = shape.z[indices]
    response = rho.^2 .+ z.^2
    design = hcat(2z, ones(length(z)))
    coefficients = design \ response
    center_z = coefficients[1]
    radius = sqrt(coefficients[2] + center_z^2)
    residual = sqrt.(rho.^2 .+ (z .- center_z).^2) .- radius
    protrusion = shape.z[1] - (center_z + radius)
    (; center_z, radius, rms=sqrt(mean(residual.^2)), protrusion)
end

function validation_metrics_v6(result)
    result.metrics === nothing && return nothing
    shape, metrics = result.shape, result.metrics
    ld_fit = ld_sphere_fit_v6(shape)
    north_radius = 2 / abs(metrics.two_m[1])
    south_radius = 2 / abs(metrics.two_m[end])
    boundary_osculating_radius = 1 / abs(metrics.ks[N_LO_SEGMENTS])
    (; contour_shape_rms=contour_shape_rms_v6(shape),
       apex_displacement=metrics.apex_displacement,
       apex_error=metrics.apex_displacement - C.target_height,
       cap_height=metrics.cap_height,
       cap_height_error=metrics.cap_height - C.target_cap_height,
       cap_diameter=metrics.cap_diameter,
       cap_diameter_error=metrics.cap_diameter - C.target_cap_diameter,
       boundary_radius=metrics.boundary_radius,
       boundary_radius_error=metrics.boundary_radius - C.boundary_radius,
       boundary_normal_tilt_deg=rad2deg(shape.psi[DB_NODE]),
       north_osculating_radius=north_radius,
       south_osculating_radius=south_radius,
       boundary_osculating_radius,
       lo_area_change_percent=100 * (metrics.area_lo / REF51.area_lo - 1),
       ld_fitted_radius=ld_fit.radius, ld_fit_rms=ld_fit.rms,
       lo_protrusion_from_ld_fit=ld_fit.protrusion)
end

function solve_state_v51(state, seed; branch)
    model = Model(Ipopt.Optimizer)
    if V51_IPOPT_PRINT_LEVEL == 0
        set_silent(model)
    else
        set_optimizer_attribute(model, "print_level", V51_IPOPT_PRINT_LEVEL)
    end
    set_optimizer_attribute(model, "tol", V51_IPOPT_TOL)
    set_optimizer_attribute(model, "constr_viol_tol", V51_IPOPT_TOL)
    set_optimizer_attribute(model, "acceptable_tol", 10V51_IPOPT_TOL)
    set_optimizer_attribute(model, "max_iter", V51_IPOPT_MAXITER)

    @variable(model, 0.0 <= rho[1:N_NODES] <= 20.0)
    @variable(model, -20.0 <= z[1:N_NODES] <= 20.0)
    @variable(model, 0.0 <= psi[1:N_NODES] <= pi)
    @variable(model, 0.25 * REF51.shape.ell_lo <= ell_lo <= 4.0 * REF51.shape.ell_lo)
    @variable(model, 0.25 * REF51.shape.ell_ld <= ell_ld <= 4.0 * REF51.shape.ell_ld)

    start = seed isa NamedTuple && hasproperty(seed, :ell_lo) ? seed : start_shape_v51(seed)
    for i in 1:N_NODES
        set_start_value(rho[i], start.rho[i]); set_start_value(z[i], start.z[i]); set_start_value(psi[i], start.psi[i])
    end
    set_start_value(ell_lo, start.ell_lo); set_start_value(ell_ld, start.ell_ld)

    fix(rho[1], 0.0; force=true); fix(rho[end], 0.0; force=true)
    fix(psi[1], 0.0; force=true); fix(psi[end], pi; force=true)
    @constraint(model, z[1] + z[end] == 0.0)
    for i in 1:N_SEGMENTS
        local_length = i <= N_LO_SEGMENTS ? ell_lo : ell_ld
        @NLconstraint(model, rho[i + 1] - rho[i] == local_length * cos((psi[i] + psi[i + 1]) / 2))
        @NLconstraint(model, z[i + 1] - z[i] == -local_length * sin((psi[i] + psi[i + 1]) / 2))
    end

    @NLexpression(model, area_lo,
        sum(pi * (rho[i] + rho[i + 1]) * ell_lo for i in 1:N_LO_SEGMENTS))
    @NLexpression(model, area_ld,
        sum(pi * (rho[i] + rho[i + 1]) * ell_ld for i in (N_LO_SEGMENTS + 1):N_SEGMENTS))
    @NLexpression(model, volume,
        sum(pi * (z[i] - z[i + 1]) * (rho[i]^2 + rho[i] * rho[i + 1] + rho[i + 1]^2) / 3 for i in 1:N_SEGMENTS))
    @NLconstraint(model, volume == REF51.volume)

    @NLexpression(model, ks_lo[i=1:N_LO_SEGMENTS], (psi[i + 1] - psi[i]) / ell_lo)
    @NLexpression(model, kp_lo[i=1:N_LO_SEGMENTS],
        sin((psi[i] + psi[i + 1]) / 2) / ((rho[i] + rho[i + 1]) / 2 + 1e-10))
    @NLexpression(model, tm_lo[i=1:N_LO_SEGMENTS], ks_lo[i] + kp_lo[i])
    @NLexpression(model, ks_ld[i=(N_LO_SEGMENTS + 1):N_SEGMENTS], (psi[i + 1] - psi[i]) / ell_ld)
    @NLexpression(model, kp_ld[i=(N_LO_SEGMENTS + 1):N_SEGMENTS],
        sin((psi[i] + psi[i + 1]) / 2) / ((rho[i] + rho[i + 1]) / 2 + 1e-10))
    @NLexpression(model, tm_ld[i=(N_LO_SEGMENTS + 1):N_SEGMENTS], ks_ld[i] + kp_ld[i])
    @NLexpression(model, delta_lo,
        2 * NEUTRAL_OFFSET_UM * sum(tm_lo[i] * pi * (rho[i] + rho[i + 1]) * ell_lo for i in 1:N_LO_SEGMENTS))
    @NLexpression(model, delta_ld,
        2 * NEUTRAL_OFFSET_UM * sum(tm_ld[i] * pi * (rho[i] + rho[i + 1]) * ell_ld for i in (N_LO_SEGMENTS + 1):N_SEGMENTS))

    area_pref_lo = REF51.area_lo * (1 + state.lo_area_expansion)
    area_pref_ld = REF51.area_ld * (1 + state.ld_area_expansion)
    delta_pref_lo = REF51.delta_lo + state.lo_mismatch * REF51.area_lo
    delta_pref_ld = REF51.delta_ld + state.ld_mismatch * REF51.area_ld
    @NLexpression(model, bend_lo, 0.5 * KAPPA_LO_KBT *
        sum((tm_lo[i] - state.c0_lo_um_inv)^2 *
            pi * (rho[i] + rho[i + 1]) * ell_lo for i in 1:N_LO_SEGMENTS))
    @NLexpression(model, bend_ld, 0.5 * KAPPA_LD_KBT *
        sum((tm_ld[i] - state.c0_ld_um_inv)^2 *
            pi * (rho[i] + rho[i + 1]) * ell_ld for i in (N_LO_SEGMENTS + 1):N_SEGMENTS))
    @NLexpression(model, ealo, KA_LO * (area_lo - area_pref_lo)^2 / (2 * REF51.area_lo))
    @NLexpression(model, eald, KA_LD * (area_ld - area_pref_ld)^2 / (2 * REF51.area_ld))
    @NLexpression(model, edlo, KDELTA_LO * (delta_lo - delta_pref_lo)^2 / (2 * REF51.area_lo))
    @NLexpression(model, edld, KDELTA_LD * (delta_ld - delta_pref_ld)^2 / (2 * REF51.area_ld))
    @NLexpression(model, eline, 2pi * LINE_TENSION_KBT_UM * rho[DB_NODE])
    @NLexpression(model, egauss, 2pi * DELTA_KAPPA_G_KBT * cos(psi[DB_NODE]))
    @NLobjective(model, Min, bend_lo + bend_ld + ealo + eald + edlo + edld + eline + egauss)

    optimize!(model)
    term, primal = termination_status(model), primal_status(model)
    has_solution = primal in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    has_solution || return (status=string(term), primal_status=string(primal), branch=branch,
        accepted=false, objective=NaN, solve_time=solve_time(model), shape=nothing, metrics=nothing)
    shape = (; rho=value.(rho), z=value.(z), psi=value.(psi),
        ell_lo=value(ell_lo), ell_ld=value(ell_ld))
    metrics = compute_metrics_v51(shape, state)
    solver_ok = term in (MOI.LOCALLY_SOLVED, MOI.OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL)
    mesh_ok = metrics.max_geometry_error < 1e-7 && abs(metrics.volume_error) < 1e-8 &&
              metrics.axial_order_violations == 0 &&
              metrics.interior_axis_contacts == 0 &&
              metrics.profile_self_intersections == 0 &&
              metrics.segment_ratio <= MAX_SEGMENT_RATIO
    (status=string(term), primal_status=string(primal), branch=branch,
        accepted=solver_ok && mesh_ok, objective=objective_value(model), solve_time=solve_time(model),
        shape, metrics)
end

function write_result_v51(prefix, result)
    if result.metrics === nothing
        CSV.write(joinpath(V51_OUT, "$(prefix)_summary.csv"), DataFrame(branch=[result.branch],
            status=[result.status], accepted=[false]))
        return
    end
    m, s = result.metrics, result.shape
    CSV.write(joinpath(V51_OUT, "$(prefix)_summary.csv"), DataFrame(branch=[result.branch],
        status=[result.status], primal_status=[result.primal_status], accepted=[result.accepted],
        objective_kbt=[result.objective], solve_time_s=[result.solve_time], volume_error=[m.volume_error],
        geometry_error=[m.max_geometry_error],
        negative_meridional_segments=[m.negative_meridional_segments],
        min_dpsi_rad=[m.min_dpsi],
        axial_order_violations=[m.axial_order_violations],
        interior_axis_contacts=[m.interior_axis_contacts],
        profile_self_intersections=[m.profile_self_intersections],
        segment_ratio=[m.segment_ratio],
        min_segment_um=[m.min_segment], max_segment_um=[m.max_segment], apex_displacement_um=[m.apex_displacement],
        cap_height_um=[m.cap_height], cap_diameter_um=[m.cap_diameter], lo_area_um2=[m.area_lo],
        ld_area_um2=[m.area_ld], lo_area_change_percent=[100 * (m.area_lo / REF51.area_lo - 1)],
        effective_c0_lo_um_inv=[m.c0_lo_um_inv], effective_c0_ld_um_inv=[m.c0_ld_um_inv],
        lo_bending_energy_kbt=[m.bend_lo], ld_bending_energy_kbt=[m.bend_ld],
        lo_area_energy_kbt=[m.area_energy_lo], ld_area_energy_kbt=[m.area_energy_ld],
        lo_mismatch_energy_kbt=[m.mismatch_energy_lo], ld_mismatch_energy_kbt=[m.mismatch_energy_ld],
        line_energy_kbt=[m.line_energy], gaussian_energy_kbt=[m.gaussian_energy],
        interface_moment_residual_kbt_um=[m.interface.moment_residual],
        interface_moment_relative_residual=[m.interface.moment_relative_residual],
        interface_transverse_residual_kbt_um2=[m.interface.transverse_residual],
        interface_transverse_relative_residual=[m.interface.transverse_relative_residual],
        total_energy_kbt=[m.total_energy]))
    CSV.write(joinpath(V51_OUT, "$(prefix)_nodes.csv"), DataFrame(node=1:N_NODES, rho_um=s.rho,
        z_um=s.z, psi_rad=s.psi, domain=[i <= DB_NODE ? "Lo" : "Ld" for i in 1:N_NODES]))
    CSV.write(joinpath(V51_OUT, "$(prefix)_segments.csv"), DataFrame(segment=1:N_SEGMENTS,
        domain=[i <= N_LO_SEGMENTS ? "Lo" : "Ld" for i in 1:N_SEGMENTS],
        length_um=m.segment_length, area_um2=m.segment_area,
        meridional_curvature_um_inv=m.ks, azimuthal_curvature_um_inv=m.kphi,
        two_mean_curvature_um_inv=m.two_m))
end

function summary_row_v51(label, result)
    m = result.metrics
    m === nothing && return (label=label, branch=result.branch, status=result.status, accepted=false,
        objective_kbt=NaN, volume_error=NaN, geometry_error=NaN,
        negative_meridional_segments=-1, min_dpsi_rad=NaN,
        axial_order_violations=-1, interior_axis_contacts=-1,
        profile_self_intersections=-1, segment_ratio=NaN,
        apex_displacement_um=NaN, cap_height_um=NaN,
        cap_diameter_um=NaN, lo_area_change_percent=NaN,
        effective_c0_lo_um_inv=NaN, effective_c0_ld_um_inv=NaN, lo_bending_energy_kbt=NaN,
        ld_bending_energy_kbt=NaN, lo_area_energy_kbt=NaN, ld_area_energy_kbt=NaN,
        lo_mismatch_energy_kbt=NaN, ld_mismatch_energy_kbt=NaN, line_energy_kbt=NaN,
        gaussian_energy_kbt=NaN, interface_moment_residual_kbt_um=NaN,
        interface_moment_relative_residual=NaN,
        interface_transverse_residual_kbt_um2=NaN,
        interface_transverse_relative_residual=NaN, total_energy_kbt=NaN)
    (label=label, branch=result.branch, status=result.status, accepted=result.accepted,
        objective_kbt=result.objective, volume_error=m.volume_error,
        geometry_error=m.max_geometry_error,
        negative_meridional_segments=m.negative_meridional_segments,
        min_dpsi_rad=m.min_dpsi,
        axial_order_violations=m.axial_order_violations,
        interior_axis_contacts=m.interior_axis_contacts,
        profile_self_intersections=m.profile_self_intersections,
        segment_ratio=m.segment_ratio,
        apex_displacement_um=m.apex_displacement,
        cap_height_um=m.cap_height, cap_diameter_um=m.cap_diameter,
        lo_area_change_percent=100 * (m.area_lo / REF51.area_lo - 1),
        effective_c0_lo_um_inv=m.c0_lo_um_inv,
        effective_c0_ld_um_inv=m.c0_ld_um_inv,
        lo_bending_energy_kbt=m.bend_lo, ld_bending_energy_kbt=m.bend_ld,
        lo_area_energy_kbt=m.area_energy_lo, ld_area_energy_kbt=m.area_energy_ld,
        lo_mismatch_energy_kbt=m.mismatch_energy_lo,
        ld_mismatch_energy_kbt=m.mismatch_energy_ld, line_energy_kbt=m.line_energy,
        gaussian_energy_kbt=m.gaussian_energy,
        interface_moment_residual_kbt_um=m.interface.moment_residual,
        interface_moment_relative_residual=m.interface.moment_relative_residual,
        interface_transverse_residual_kbt_um2=m.interface.transverse_residual,
        interface_transverse_relative_residual=m.interface.transverse_relative_residual,
        total_energy_kbt=m.total_energy)
end

function write_validation_v6(result)
    validation = validation_metrics_v6(result)
    validation === nothing && return
    CSV.write(joinpath(V51_OUT, "post_validation.csv"), DataFrame(
        branch=[result.branch],
        normalized_contour_rms=[validation.contour_shape_rms],
        apex_displacement_um=[validation.apex_displacement],
        apex_target_um=[C.target_height],
        apex_error_um=[validation.apex_error],
        cap_height_um=[validation.cap_height],
        cap_height_target_um=[C.target_cap_height],
        cap_height_error_um=[validation.cap_height_error],
        cap_diameter_um=[validation.cap_diameter],
        cap_diameter_target_um=[C.target_cap_diameter],
        cap_diameter_error_um=[validation.cap_diameter_error],
        boundary_radius_um=[validation.boundary_radius],
        boundary_radius_target_um=[C.boundary_radius],
        boundary_radius_error_um=[validation.boundary_radius_error],
        boundary_normal_tilt_deg=[validation.boundary_normal_tilt_deg],
        north_osculating_radius_um=[validation.north_osculating_radius],
        south_osculating_radius_um=[validation.south_osculating_radius],
        boundary_osculating_radius_um=[validation.boundary_osculating_radius],
        lo_area_change_percent=[validation.lo_area_change_percent],
        ld_fitted_radius_um=[validation.ld_fitted_radius],
        ld_fit_rms_um=[validation.ld_fit_rms],
        lo_protrusion_from_ld_fit_um=[validation.lo_protrusion_from_ld_fit]))
end

function strep_coverage_v770(concentration_nm)
    concentration_nm >= 0 || error("Streptavidin concentration must be nonnegative.")
    clamp(concentration_nm / V7510_STREP_SATURATION_NM, 0.0, 1.0)
end

function incremental_inner_ld_area_fraction_v770(concentration_nm)
    coverage = strep_coverage_v770(concentration_nm)
    reference_coverage = strep_coverage_v770(10.0)
    coverage + 1e-12 >= reference_coverage ||
        error("The suppression sweep cannot go below the 10 nM reference.")
    INNER_LD_BIOTIN_MOLE_FRACTION * BOUND_SITE_AREA_COUPLING *
        max(coverage - reference_coverage, 0.0)
end

function states_v7510()
    (
        baseline=V7State("baseline", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        lipid_only=V7State("lipid_asymmetry_only", 0.0, 0.0, 0.0, 0.0,
            C0_LO_LIPID_UM_INV, C0_LD_LIPID_UM_INV),
        strep_only=V7State("strep10_only", 0.0, 0.0, 0.0, 0.0,
            0.0, C0_LD_STREP10_UM_INV),
        combined_pre=V7State("lipid_plus_strep10_pre_ctb", 0.0, 0.0, 0.0, 0.0,
            C0_LO_LIPID_UM_INV, C0_LD_COMBINED10_UM_INV),
        combined_post=V7State("lipid_plus_strep10_post_ctb", EPS_CTB_OUTER_LO / 2,
            LO_GEOMETRIC_MISMATCH_V71, 0.0, 0.0,
            C0_LO_POST_UM_INV, C0_LD_COMBINED10_UM_INV),
    )
end

function concentration_slug_v770(concentration_nm)
    replace(string(concentration_nm), "." => "p")
end

function post_state_v770(concentration_nm)
    inner_ld_area_fraction =
        incremental_inner_ld_area_fraction_v770(concentration_nm)
    V7State(
        "post_ctb_strep_$(concentration_slug_v770(concentration_nm))nM",
        EPS_CTB_OUTER_LO / 2,
        LO_GEOMETRIC_MISMATCH_V71,
        inner_ld_area_fraction / 2,
        -inner_ld_area_fraction,
        C0_LO_POST_UM_INV,
        C0_LD_COMBINED10_UM_INV,
    )
end

function write_audit_v51()
    states = states_v7510()
    labels = [states.combined_pre.label, states.combined_post.label]
    ctb = [false, true]
    strep_nm = [10.0, 10.0]
    c0_ld_lipid = fill(C0_LD_LIPID_UM_INV, 2)
    c0_ld_strep = fill(C0_LD_STREP10_UM_INV, 2)
    c0_lo_lipid = fill(C0_LO_LIPID_UM_INV, 2)
    c0_lo_ctb = [0.0, C0_LO_CTB_UM_INV]
    effective_ld = c0_ld_lipid .+ c0_ld_strep
    effective_lo = c0_lo_lipid .+ c0_lo_ctb
    lo_area_expansion = [0.0, EPS_CTB_OUTER_LO / 2]
    lo_mismatch = [0.0, LO_GEOMETRIC_MISMATCH_V71]

    CSV.write(joinpath(V51_OUT, "condition_audit.csv"), DataFrame(
        field=[
            "purpose", "shape mesh", "tangent-angle admissibility",
            "profile validity", "volume constraint", "area mechanics",
            "curvature energy", "curvature sign convention",
            "interface mechanics", "interface residual estimator",
            "parameter transfer status", "direct coupling", "protein work",
            "target-shape force", "validation observable", "Ipopt Hessian",
        ],
        value=[
            "test concentration-dependent local inner-Ld leaflet strain with the corrected natural two-phase interface",
            "$(N_NODES)-node axisymmetric meridian",
            "0 <= psi <= pi; local decreases are allowed and represent negative meridional curvature",
            "axial ordering, interior-axis contacts, and nonadjacent profile intersections checked after solve",
            "fixed discrete reference volume using exact conical frusta",
            "fixed CTB Lo area/mismatch plus concentration-dependent preferred inner-Ld area",
            "domain-local Helfrich energy (2H-C0_domain)^2",
            "outward normal: sphere 2H>0; outer-leaflet sources C0>0; inner-leaflet sources C0<0",
            "not added as duplicate constraints; follows variationally from the full two-phase energy",
            "one-sided linear extrapolation of 2H plus moment and transverse-shear residuals",
            "relative magnitudes from four-control study; scale fixed by V7.3 effective Ld curvature",
            "inactive; domains communicate through shared shape, volume, and boundary",
            "inactive", "inactive",
            "post-minus-pre Lo protrusion relative to the fitted Ld sphere",
            "exact",
        ],
    ))
    CSV.write(joinpath(V51_OUT, "run_config.csv"), DataFrame(
        parameter=[
            "mode", "n_nodes", "ipopt_tol", "ipopt_maxiter",
            "psi_monotonicity_constraint", "geometry_tolerance",
            "interior_axis_contact_tolerance_um",
            "kappa_lo_kbt", "kappa_ld_kbt", "line_tension_pn",
            "gaussian_ratio", "delta_kappa_g_kbt",
            "prior_c0_ld_lipid_magnitude_um_inv",
            "prior_c0_ld_strep10_magnitude_um_inv",
            "prior_c0_lo_lipid_magnitude_um_inv",
            "ratio_preserving_transfer_scale",
            "target_c0_ld_combined10_um_inv",
            "target_c0_lo_post_um_inv",
            "c0_ld_lipid_um_inv", "c0_lo_lipid_um_inv",
            "c0_ld_strep10_um_inv", "c0_ld_combined10_um_inv",
            "c0_lo_ctb_um_inv", "c0_lo_post_total_um_inv",
            "strep_saturation_nm",
            "inner_ld_biotin_mole_fraction",
            "bound_site_area_coupling",
            "ctb_lo_area_expansion", "ctb_lo_leaflet_mismatch",
            "pre_sphere_tolerance_um",
            "corrected_dome_increase_target_um",
            "corrected_dome_acceptance_tolerance_um",
            "parameter_provenance", "target_profile_path",
        ],
        value=[
            V51_MODE, N_NODES, V51_IPOPT_TOL, V51_IPOPT_MAXITER,
            false, V760_GEOMETRY_TOL, V760_AXIS_CONTACT_TOL_UM,
            KAPPA_LO_KBT, KAPPA_LD_KBT, LINE_TENSION_PN,
            V73_GAUSSIAN_RATIO, DELTA_KAPPA_G_KBT,
            PRIOR_C0_LD_LIPID_MAG_UM_INV,
            PRIOR_C0_LD_STREP10_MAG_UM_INV,
            PRIOR_C0_LO_LIPID_MAG_UM_INV,
            V7510_TRANSFER_SCALE,
            TARGET_C0_LD_COMBINED10_UM_INV,
            TARGET_C0_LO_POST_UM_INV,
            C0_LD_LIPID_UM_INV, C0_LO_LIPID_UM_INV,
            C0_LD_STREP10_UM_INV, C0_LD_COMBINED10_UM_INV,
            C0_LO_CTB_UM_INV, C0_LO_POST_UM_INV,
            V7510_STREP_SATURATION_NM,
            INNER_LD_BIOTIN_MOLE_FRACTION,
            BOUND_SITE_AREA_COUPLING,
            EPS_CTB_OUTER_LO / 2, LO_GEOMETRIC_MISMATCH_V71,
            PRE_SPHERE_TOL_UM,
            EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM,
            CORRECTED_DOME_ACCEPTANCE_TOL_UM,
            "V7.3/V7.4 10 nM curvatures retained; additional occupancy acts only through preferred inner-Ld area",
            TARGET_PROFILE_PATH,
        ],
    ))
    CSV.write(joinpath(V51_OUT, "state_table.csv"), DataFrame(
        state=labels,
        ctb_present=ctb,
        strep_concentration_nm=strep_nm,
        c0_ld_lipid_um_inv=c0_ld_lipid,
        c0_ld_strep_um_inv=c0_ld_strep,
        c0_lo_lipid_um_inv=c0_lo_lipid,
        c0_lo_ctb_um_inv=c0_lo_ctb,
        effective_c0_ld_um_inv=effective_ld,
        effective_c0_lo_um_inv=effective_lo,
        lo_area_expansion=lo_area_expansion,
        lo_leaflet_mismatch=lo_mismatch,
        ld_area_expansion=zeros(2),
        ld_leaflet_mismatch=zeros(2),
        minimized=trues(2),
    ))
    CSV.write(joinpath(V51_OUT, "energy_term_audit.csv"), DataFrame(
        term=[
            "Lo Helfrich bending", "Ld Helfrich bending",
            "Lo area elasticity", "Ld area elasticity",
            "CTB shift of preferred Lo area",
            "Lo area-difference elasticity", "Ld area-difference elasticity",
            "CTB shift of preferred Lo area difference",
            "Lo/Ld line energy", "Gaussian boundary energy",
            "direct Lo-Ld curvature coupling", "protein work", "target-shape force",
        ],
        combined_pre=[true, true, true, true, false, true, true, false,
            true, true, false, false, false],
        combined_post=[true, true, true, true, true, true, true, true,
            true, true, false, false, false],
    ))
    CSV.write(joinpath(V51_OUT, "parameter_role_audit.csv"), DataFrame(
        parameter=[
            "C0_Ld_lipid", "C0_Lo_lipid", "C0_Ld_Strep10", "C0_Lo_CTB",
            "CTB Lo preferred-area shift", "CTB Lo area-difference shift",
            "strep concentration label", "inner-Ld biotin mole fraction",
            "direct Lo-Ld coupling coefficient", "protein work coefficient",
        ],
        value=[
            C0_LD_LIPID_UM_INV, C0_LO_LIPID_UM_INV,
            C0_LD_STREP10_UM_INV, C0_LO_CTB_UM_INV,
            EPS_CTB_OUTER_LO / 2, LO_GEOMETRIC_MISMATCH_V71,
            10.0, INNER_LD_BIOTIN_MOLE_FRACTION, 0.0, 0.0,
        ],
        energy_active=[
            true, true, true, true, true, true,
            false, false, false, false,
        ],
        metadata_only=[
            false, false, false, false, false, false,
            true, true, false, false,
        ],
        role=[
            "Ld Helfrich target curvature",
            "Lo Helfrich target curvature",
            "Ld Helfrich target curvature",
            "Lo Helfrich target curvature",
            "changes preferred Lo area after CTB",
            "changes preferred Lo area difference after CTB",
            "condition identifier; concentration mapping is audited separately",
            "composition record; not converted independently into energy",
            "explicitly inactive",
            "explicitly inactive",
        ],
    ))
    concentration_nm = [10.0, 20.0, 40.0, 80.0]
    c0_strep = fill(C0_LD_STREP10_UM_INV, length(concentration_nm))
    incremental_inner_ld_area =
        incremental_inner_ld_area_fraction_v770.(concentration_nm)
    CSV.write(joinpath(V51_OUT, "concentration_mapping.csv"), DataFrame(
        strep_concentration_nm=concentration_nm,
        relative_coverage=strep_coverage_v770.(concentration_nm),
        incremental_inner_ld_area_fraction=incremental_inner_ld_area,
        ld_midplane_preferred_area_fraction_shift=
            incremental_inner_ld_area ./ 2,
        ld_outer_minus_inner_preferred_area_fraction_shift=
            -incremental_inner_ld_area,
        c0_ld_lipid_um_inv=fill(C0_LD_LIPID_UM_INV, length(concentration_nm)),
        c0_ld_strep_um_inv=c0_strep,
        effective_c0_ld_um_inv=C0_LD_LIPID_UM_INV .+ c0_strep,
        c0_lo_lipid_um_inv=fill(C0_LO_LIPID_UM_INV, length(concentration_nm)),
        c0_lo_ctb_um_inv=fill(C0_LO_CTB_UM_INV, length(concentration_nm)),
        effective_c0_lo_post_um_inv=fill(C0_LO_POST_UM_INV, length(concentration_nm)),
        energy_active=trues(length(concentration_nm)),
        role=fill(
            "concentration changes only the preferred inner-Ld leaflet area; spontaneous curvatures remain fixed",
            length(concentration_nm)),
    ))
    CSV.write(joinpath(V51_OUT, "expected_outputs.csv"), DataFrame(
        path=[
            "state_table.csv", "energy_term_audit.csv",
            "parameter_role_audit.csv", "concentration_mapping.csv",
            "branch_summary.csv",
            "selected_states.csv", "*_nodes.csv", "*_segments.csv",
            "interface_balance.csv", "reference_validation.csv",
            "status.txt",
        ],
        created_by=[
            "audit", "audit", "audit", "audit", "smoke", "smoke", "smoke", "smoke",
            "smoke when a post-CTB branch is feasible",
            "smoke when pre- and post-CTB branches are feasible",
            "audit or smoke",
        ],
    ))
end

function write_audit_v770()
    write_audit_v51()

    concentrations = V770_CONCENTRATIONS_NM
    coverage = strep_coverage_v770.(concentrations)
    incremental_inner_ld_area =
        incremental_inner_ld_area_fraction_v770.(concentrations)
    c0_strep = fill(C0_LD_STREP10_UM_INV, length(concentrations))
    CSV.write(joinpath(V51_OUT, "concentration_mapping.csv"), DataFrame(
        strep_concentration_nm=concentrations,
        bound_fraction_of_inner_ld_biotin=coverage,
        inner_ld_biotin_mole_fraction=
            fill(INNER_LD_BIOTIN_MOLE_FRACTION, length(concentrations)),
        bound_site_area_coupling=
            fill(BOUND_SITE_AREA_COUPLING, length(concentrations)),
        incremental_inner_ld_area_fraction=incremental_inner_ld_area,
        ld_midplane_preferred_area_fraction_shift=
            incremental_inner_ld_area ./ 2,
        ld_outer_minus_inner_preferred_area_fraction_shift=
            -incremental_inner_ld_area,
        c0_ld_lipid_um_inv=fill(C0_LD_LIPID_UM_INV, length(concentrations)),
        c0_ld_strep_um_inv=c0_strep,
        effective_c0_ld_um_inv=C0_LD_LIPID_UM_INV .+ c0_strep,
        c0_lo_lipid_um_inv=fill(C0_LO_LIPID_UM_INV, length(concentrations)),
        c0_lo_ctb_um_inv=fill(C0_LO_CTB_UM_INV, length(concentrations)),
        effective_c0_lo_post_um_inv=fill(C0_LO_POST_UM_INV, length(concentrations)),
        mapping_enters_energy=trues(length(concentrations)),
        interpretation=fill(
            "concentration shifts the preferred inner-Ld leaflet area; the 10 nM spontaneous curvatures remain fixed",
            length(concentrations)),
    ))

    state_names = ["pre_ctb_strep_10nM";
        ["post_ctb_strep_$(concentration_slug_v770(c))nM" for c in concentrations]]
    ctb_present = [false; trues(length(concentrations))]
    strep_concentration = [10.0; concentrations]
    c0_ld_strep = [C0_LD_STREP10_UM_INV; c0_strep]
    c0_lo_ctb = [0.0; fill(C0_LO_CTB_UM_INV, length(concentrations))]
    lo_area_expansion = [0.0; fill(EPS_CTB_OUTER_LO / 2, length(concentrations))]
    lo_mismatch = [0.0; fill(LO_GEOMETRIC_MISMATCH_V71, length(concentrations))]
    ld_area_expansion = [0.0; incremental_inner_ld_area ./ 2]
    ld_mismatch = [0.0; -incremental_inner_ld_area]
    CSV.write(joinpath(V51_OUT, "state_table.csv"), DataFrame(
        state=state_names,
        ctb_present=ctb_present,
        strep_concentration_nm=strep_concentration,
        c0_ld_lipid_um_inv=fill(C0_LD_LIPID_UM_INV, length(state_names)),
        c0_ld_strep_um_inv=c0_ld_strep,
        effective_c0_ld_um_inv=C0_LD_LIPID_UM_INV .+ c0_ld_strep,
        c0_lo_lipid_um_inv=fill(C0_LO_LIPID_UM_INV, length(state_names)),
        c0_lo_ctb_um_inv=c0_lo_ctb,
        effective_c0_lo_um_inv=C0_LO_LIPID_UM_INV .+ c0_lo_ctb,
        lo_area_expansion=lo_area_expansion,
        lo_leaflet_mismatch=lo_mismatch,
        ld_area_expansion=ld_area_expansion,
        ld_leaflet_mismatch=ld_mismatch,
        minimized=[V51_MODE == "smoke"; fill(V51_MODE == "smoke", length(concentrations))],
    ))

    CSV.write(joinpath(V51_OUT, "suppression_model_audit.csv"), DataFrame(
        item=[
            "independent variable",
            "energy-active concentration mapping",
            "fixed CTB Lo curvature",
            "fixed CTB Lo preferred-area shift",
            "fixed CTB Lo area-difference shift",
            "fixed lipid spontaneous curvatures",
            "fixed mechanical moduli",
            "shared geometric coupling",
            "inactive occupancy energy",
            "inactive protein work",
            "inactive direct Lo-Ld coupling",
            "inactive target-shape force",
            "branch selection",
            "mesh",
        ],
        value=[
            "streptavidin concentration (nM)",
            "epsilon_inner,Ld(c)=f_biotin*g_area*[theta(c)-theta(10 nM)]",
            string(C0_LO_CTB_UM_INV, " um^-1"),
            string(EPS_CTB_OUTER_LO / 2),
            string(LO_GEOMETRIC_MISMATCH_V71),
            "C0_Ld_lipid=$(C0_LD_LIPID_UM_INV) um^-1; C0_Lo_lipid=$(C0_LO_LIPID_UM_INV) um^-1",
            "kappa_Lo=$(KAPPA_LO_KBT) kBT; kappa_Ld=$(KAPPA_LD_KBT) kBT; line tension=$(LINE_TENSION_PN) pN",
            "shared boundary geometry; force and moment balances follow from variation of domain, line, and Gaussian energies",
            "zero",
            "zero",
            "zero",
            "zero",
            "lowest objective among accepted continuation, pre-reference, and domed seeds",
            "$(N_NODES)-node axisymmetric meridian",
        ],
    ))

    CSV.write(joinpath(V51_OUT, "parameter_role_audit.csv"), DataFrame(
        parameter=[
            "streptavidin concentration",
            "inner-Ld biotin mole fraction",
            "C0_Ld_lipid",
            "C0_Ld_Strep(10 nM reference)",
            "C0_Lo_lipid",
            "C0_Lo_CTB",
            "CTB Lo preferred-area shift",
            "CTB Lo area-difference shift",
            "incremental inner-Ld preferred-area shift",
            "occupancy free energy",
            "protein work",
            "direct Lo-Ld coupling",
            "target-shape force",
        ],
        energy_active=[
            true, true, true, true, true, true, true, true, true,
            false, false, false, false,
        ],
        metadata_only=[
            false, false, false, false, false, false, false, false, false,
            false, false, false, false,
        ],
        role=[
            "mapped to bound-site fraction and then to preferred inner-Ld area",
            "sets the maximum number of area-active inner-Ld binding sites",
            "fixed Ld Helfrich target curvature",
            "fixed Ld Helfrich target curvature at its validated 10 nM value",
            "fixed Lo Helfrich target curvature",
            "fixed post-CTB Lo Helfrich target curvature",
            "fixed post-CTB preferred Lo area",
            "fixed post-CTB preferred Lo area difference",
            "changes Ld preferred midplane area by +epsilon/2 and outer-minus-inner area difference by -epsilon",
            "explicitly inactive",
            "explicitly inactive",
            "explicitly inactive",
            "explicitly inactive",
        ],
    ))

    CSV.write(joinpath(V51_OUT, "energy_term_audit.csv"), DataFrame(
        term=[
            "Lo Helfrich bending",
            "Ld Helfrich bending",
            "Lo midplane-area elasticity",
            "Ld midplane-area elasticity",
            "Lo leaflet-area-difference elasticity",
            "Ld leaflet-area-difference elasticity",
            "CTB shift of preferred Lo area",
            "CTB shift of preferred Lo area difference",
            "Strep shift of preferred inner-Ld area",
            "Lo/Ld line energy",
            "Gaussian boundary energy",
            "exact enclosed-volume constraint",
            "direct Lo-Ld curvature coupling",
            "protein work",
            "target-shape force",
            "global ADE",
        ],
        pre_ctb_10nM=[
            true, true, true, true, true, true,
            false, false, false, true, true, true,
            false, false, false, false,
        ],
        post_ctb_10nM=[
            true, true, true, true, true, true,
            true, true, false, true, true, true,
            false, false, false, false,
        ],
        post_ctb_above_10nM=[
            true, true, true, true, true, true,
            true, true, true, true, true, true,
            false, false, false, false,
        ],
        implementation=[
            "0.5*kappa_Lo*integral(2H-C0_Lo)^2 dA",
            "0.5*kappa_Ld*integral(2H-C0_Ld)^2 dA",
            "KA_Lo*(A_Lo-A0_Lo)^2/(2*A0_Lo)",
            "KA_Ld*(A_Ld-A0_Ld)^2/(2*A0_Ld)",
            "KDelta_Lo*(DeltaA_Lo-DeltaA0_Lo)^2/(2*A0_Lo)",
            "KDelta_Ld*(DeltaA_Ld-DeltaA0_Ld)^2/(2*A0_Ld)",
            "A0_Lo -> A0_Lo*(1+epsilon_CTB/2)",
            "DeltaA0_Lo -> DeltaA0_Lo+epsilon_geometry*A0_Lo",
            "A0_Ld -> A0_Ld*(1+epsilon_inner/2); DeltaA0_Ld -> DeltaA0_Ld-epsilon_inner*A0_Ld",
            "lambda*L_LoLd",
            "2*pi*DeltaKappaG*cos(psi_boundary)",
            "V=V_reference",
            "inactive",
            "inactive",
            "inactive",
            "inactive",
        ],
    ))

    CSV.write(joinpath(V51_OUT, "branch_plan.csv"), DataFrame(
        branch=["continuation", "pre_reference", "fig7_domed_seed"],
        initial_state=[
            "previous selected concentration; pre-CTB selected state at 10 nM",
            "selected pre-CTB 10 nM state",
            "independent Fig. 7 domed seed",
        ],
        selection_rule=fill(
            "lowest objective among solver- and geometry-accepted branches",
            3),
    ))

    CSV.write(joinpath(V51_OUT, "expected_outputs.csv"), DataFrame(
        path=[
            "state_table.csv",
            "concentration_mapping.csv",
            "suppression_model_audit.csv",
            "parameter_role_audit.csv",
            "energy_term_audit.csv",
            "branch_plan.csv",
            "candidate_branch_summary.csv",
            "selected_reference.csv",
            "selected_concentration_trace.csv",
            "selected_interface_balance.csv",
            "smoke_validation.csv",
            "*_summary.csv",
            "*_nodes.csv",
            "*_segments.csv",
            "status.txt",
        ],
        created_by=[
            "audit", "audit", "audit", "audit", "audit", "audit",
            "smoke", "smoke", "smoke", "smoke", "smoke",
            "smoke for every candidate", "smoke for every candidate",
            "smoke for every candidate", "audit or smoke",
        ],
    ))
end

function solve_two_seeds_v7510(label, state, sphere_seed)
    sphere = solve_state_v51(state, sphere_seed; branch="near_sphere")
    dome = solve_state_v51(state, domed_seed(); branch="fig7_domed_seed")
    ((label * "_sphere", sphere), (label * "_dome", dome))
end

function select_lowest_v7510(candidates)
    valid = [result for (_, result) in candidates if result.accepted]
    isempty(valid) && return nothing
    valid[argmin([result.objective for result in valid])]
end

function selected_row_v7510(label, result)
    result === nothing && return (
        state=label, branch="none", accepted=false, objective_kbt=NaN,
        apex_displacement_um=NaN, radial_cv=NaN, sphericity=NaN,
        cap_height_um=NaN, cap_diameter_um=NaN,
        effective_c0_ld_um_inv=NaN, effective_c0_lo_um_inv=NaN,
        sphere_error_um=NaN, lo_protrusion_from_ld_sphere_um=NaN,
        negative_meridional_segments=-1,
        interface_moment_residual_kbt_um=NaN,
        interface_moment_relative_residual=NaN,
        interface_transverse_residual_kbt_um2=NaN,
        interface_transverse_relative_residual=NaN,
    )
    metrics = result.metrics
    (
        state=label,
        branch=result.branch,
        accepted=result.accepted,
        objective_kbt=result.objective,
        apex_displacement_um=metrics.apex_displacement,
        radial_cv=metrics.radial_cv,
        sphericity=metrics.sphericity,
        cap_height_um=metrics.cap_height,
        cap_diameter_um=metrics.cap_diameter,
        effective_c0_ld_um_inv=metrics.c0_ld_um_inv,
        effective_c0_lo_um_inv=metrics.c0_lo_um_inv,
        sphere_error_um=max(abs(metrics.apex_displacement), C.radius * metrics.radial_cv),
        lo_protrusion_from_ld_sphere_um=ld_sphere_fit_v6(result.shape).protrusion,
        negative_meridional_segments=metrics.negative_meridional_segments,
        interface_moment_residual_kbt_um=metrics.interface.moment_residual,
        interface_moment_relative_residual=metrics.interface.moment_relative_residual,
        interface_transverse_residual_kbt_um2=metrics.interface.transverse_residual,
        interface_transverse_relative_residual=metrics.interface.transverse_relative_residual,
    )
end

function concentration_candidates_v770(concentration_nm, state, previous_shape, pre_shape)
    slug = concentration_slug_v770(concentration_nm)
    seeds = [
        ("continuation", previous_shape),
        ("pre_reference", pre_shape),
        ("fig7_domed_seed", domed_seed()),
    ]
    [
        ("c$(slug)_$(seed_label)",
            solve_state_v51(state, seed;
                branch="c$(slug)_$(seed_label)"))
        for (seed_label, seed) in seeds
    ]
end

function candidate_summary_row_v770(concentration_nm, label, result)
    merge(
        (strep_concentration_nm=concentration_nm,),
        summary_row_v51(label, result),
    )
end

function selected_reference_row_v770(pre)
    fit = ld_sphere_fit_v6(pre.shape)
    m = pre.metrics
    (
        state="pre_ctb_strep_10nM",
        branch=pre.branch,
        accepted=pre.accepted,
        objective_kbt=pre.objective,
        effective_c0_ld_um_inv=m.c0_ld_um_inv,
        effective_c0_lo_um_inv=m.c0_lo_um_inv,
        sphericity=m.sphericity,
        radial_cv=m.radial_cv,
        lo_protrusion_from_ld_sphere_um=fit.protrusion,
        cap_height_um=m.cap_height,
        cap_diameter_um=m.cap_diameter,
        lo_bending_energy_kbt=m.bend_lo,
        ld_bending_energy_kbt=m.bend_ld,
        lo_area_energy_kbt=m.area_energy_lo,
        ld_area_energy_kbt=m.area_energy_ld,
        lo_mismatch_energy_kbt=m.mismatch_energy_lo,
        ld_mismatch_energy_kbt=m.mismatch_energy_ld,
        line_energy_kbt=m.line_energy,
        gaussian_energy_kbt=m.gaussian_energy,
        total_energy_kbt=m.total_energy,
        interface_moment_relative_residual=m.interface.moment_relative_residual,
        interface_transverse_relative_residual=m.interface.transverse_relative_residual,
    )
end

function selected_trace_row_v770(
    concentration_nm,
    result,
    pre_protrusion,
    ten_nm_protrusion,
)
    m = result.metrics
    fit = ld_sphere_fit_v6(result.shape)
    area_tension_lo = KA_LO * (m.area_lo - m.area_pref_lo) / REF51.area_lo
    area_tension_ld = KA_LD * (m.area_ld - m.area_pref_ld) / REF51.area_ld
    kbt_um2_to_pn_um = KBT_J * 1e18
    corrected_from_pre = fit.protrusion - pre_protrusion
    suppression_from_10 = ten_nm_protrusion - fit.protrusion
    incremental_inner_ld_area =
        incremental_inner_ld_area_fraction_v770(concentration_nm)
    (
        strep_concentration_nm=concentration_nm,
        bound_fraction_of_inner_ld_biotin=
            strep_coverage_v770(concentration_nm),
        inner_ld_biotin_mole_fraction=INNER_LD_BIOTIN_MOLE_FRACTION,
        incremental_inner_ld_area_fraction=incremental_inner_ld_area,
        branch=result.branch,
        accepted=result.accepted,
        objective_kbt=result.objective,
        c0_ld_lipid_um_inv=C0_LD_LIPID_UM_INV,
        c0_ld_strep_um_inv=C0_LD_STREP10_UM_INV,
        effective_c0_ld_um_inv=m.c0_ld_um_inv,
        c0_lo_lipid_um_inv=C0_LO_LIPID_UM_INV,
        c0_lo_ctb_um_inv=C0_LO_CTB_UM_INV,
        effective_c0_lo_um_inv=m.c0_lo_um_inv,
        lo_area_expansion=EPS_CTB_OUTER_LO / 2,
        lo_leaflet_mismatch=LO_GEOMETRIC_MISMATCH_V71,
        ld_area_expansion=incremental_inner_ld_area / 2,
        ld_leaflet_mismatch=-incremental_inner_ld_area,
        sphericity=m.sphericity,
        radial_cv=m.radial_cv,
        apex_displacement_um=m.apex_displacement,
        raw_cap_height_um=m.cap_height,
        cap_diameter_um=m.cap_diameter,
        boundary_radius_um=m.boundary_radius,
        boundary_length_um=2pi * m.boundary_radius,
        ld_fitted_sphere_radius_um=fit.radius,
        ld_sphere_fit_rms_um=fit.rms,
        lo_protrusion_from_ld_sphere_um=fit.protrusion,
        corrected_lo_protrusion_from_pre_um=corrected_from_pre,
        lo_dome_suppression_from_10nM_um=suppression_from_10,
        restored_to_pre_within_0p25um=abs(corrected_from_pre) <= PRE_SPHERE_TOL_UM,
        lo_area_um2=m.area_lo,
        ld_area_um2=m.area_ld,
        lo_area_change_percent=100 * (m.area_lo / REF51.area_lo - 1),
        ld_area_change_percent=100 * (m.area_ld / REF51.area_ld - 1),
        lo_delta_area_um2=m.delta_area_lo,
        ld_delta_area_um2=m.delta_area_ld,
        lo_bending_energy_kbt=m.bend_lo,
        ld_bending_energy_kbt=m.bend_ld,
        lo_area_energy_kbt=m.area_energy_lo,
        ld_area_energy_kbt=m.area_energy_ld,
        lo_mismatch_energy_kbt=m.mismatch_energy_lo,
        ld_mismatch_energy_kbt=m.mismatch_energy_ld,
        line_energy_kbt=m.line_energy,
        gaussian_energy_kbt=m.gaussian_energy,
        total_energy_kbt=m.total_energy,
        lo_area_tension_kbt_um2=area_tension_lo,
        ld_area_tension_kbt_um2=area_tension_ld,
        lo_area_tension_pn_um=area_tension_lo * kbt_um2_to_pn_um,
        ld_area_tension_pn_um=area_tension_ld * kbt_um2_to_pn_um,
        lo_mismatch_stress_kbt_um2=m.mismatch_stress_lo,
        ld_mismatch_stress_kbt_um2=m.mismatch_stress_ld,
        lo_mismatch_stress_pn_um=m.mismatch_stress_lo * kbt_um2_to_pn_um,
        ld_mismatch_stress_pn_um=m.mismatch_stress_ld * kbt_um2_to_pn_um,
        line_tension_pn=LINE_TENSION_PN,
        negative_meridional_segments=m.negative_meridional_segments,
        min_dpsi_rad=m.min_dpsi,
        volume_relative_error=m.volume_error,
        geometry_error=m.max_geometry_error,
        axial_order_violations=m.axial_order_violations,
        interior_axis_contacts=m.interior_axis_contacts,
        profile_self_intersections=m.profile_self_intersections,
        segment_ratio=m.segment_ratio,
        interface_moment_residual_kbt_um=m.interface.moment_residual,
        interface_moment_relative_residual=m.interface.moment_relative_residual,
        interface_transverse_residual_kbt_um2=m.interface.transverse_residual,
        interface_transverse_relative_residual=m.interface.transverse_relative_residual,
    )
end

function selected_interface_row_v770(label, concentration_nm, result)
    metrics = result.metrics
    i = metrics.interface
    area_tension_lo =
        KA_LO * (metrics.area_lo - metrics.area_pref_lo) / REF51.area_lo
    area_tension_ld =
        KA_LD * (metrics.area_ld - metrics.area_pref_ld) / REF51.area_ld
    (
        state=label,
        strep_concentration_nm=concentration_nm,
        branch=result.branch,
        two_h_lo_boundary_um_inv=i.two_h_lo_boundary,
        two_h_ld_boundary_um_inv=i.two_h_ld_boundary,
        lo_effective_moment_kbt_um=i.moment_lo,
        ld_effective_moment_kbt_um=i.moment_ld,
        gaussian_moment_kbt_um=i.gaussian_moment,
        lo_area_tension_kbt_um2=area_tension_lo,
        ld_area_tension_kbt_um2=area_tension_ld,
        lo_area_difference_stress_kbt_um2=metrics.mismatch_stress_lo,
        ld_area_difference_stress_kbt_um2=metrics.mismatch_stress_ld,
        lo_leaflet_moment_contribution_kbt_um=
            2 * NEUTRAL_OFFSET_UM * metrics.mismatch_stress_lo,
        ld_leaflet_moment_contribution_kbt_um=
            2 * NEUTRAL_OFFSET_UM * metrics.mismatch_stress_ld,
        moment_residual_kbt_um=i.moment_residual,
        moment_relative_residual=i.moment_relative_residual,
        transverse_residual_kbt_um2=i.transverse_residual,
        transverse_relative_residual=i.transverse_relative_residual,
    )
end

function main_v760()
    V51_MODE in ("audit", "smoke") ||
        error("FIG7_V770_MODE must be audit or smoke.")
    PRIOR_C0_LD_LIPID_MAG_UM_INV > 0 ||
        error("Prior Ld lipid-curvature magnitude must be positive.")
    PRIOR_C0_LD_STREP10_MAG_UM_INV > 0 ||
        error("Prior Ld Strep-curvature magnitude must be positive.")
    PRIOR_C0_LO_LIPID_MAG_UM_INV >= 0 ||
        error("Prior Lo lipid-curvature magnitude must be nonnegative.")
    TARGET_C0_LD_COMBINED10_UM_INV < 0 ||
        error("The validated combined Ld target must be negative.")
    TARGET_C0_LO_POST_UM_INV > 0 ||
        error("The validated post-CTB Lo target must be positive.")
    EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM > 0 ||
        error("The corrected dome increase target must be positive.")
    CORRECTED_DOME_ACCEPTANCE_TOL_UM > 0 ||
        error("The corrected dome acceptance tolerance must be positive.")
    V7510_STREP_SATURATION_NM >= 80 ||
        error("The saturation concentration must include the 80 nM endpoint.")
    C0_LD_LIPID_UM_INV <= 0 ||
        error("Inner-Ld lipid curvature must be nonpositive in this convention.")
    C0_LD_STREP10_UM_INV <= 0 ||
        error("Inner-Ld Strep curvature must be nonpositive in this convention.")
    C0_LO_LIPID_UM_INV >= 0 ||
        error("Outer-Lo lipid curvature must be nonnegative in this convention.")
    C0_LO_CTB_UM_INV >= 0 ||
        error("Outer-Lo CTB curvature must be nonnegative in this convention.")
    isapprox(C0_LD_COMBINED10_UM_INV,
        TARGET_C0_LD_COMBINED10_UM_INV; atol=1e-12) ||
        error("The decomposed 10 nM Ld curvature does not preserve its target.")
    isapprox(C0_LO_POST_UM_INV,
        TARGET_C0_LO_POST_UM_INV; atol=1e-12) ||
        error("The decomposed post-CTB Lo curvature does not preserve its target.")
    ensure_empty_v51_output()
    write_audit_v51()
    if V51_MODE == "audit"
        write(joinpath(V51_OUT, "status.txt"),
            "audit_only: no minimization was run.\n")
        return
    end

    states = states_v7510()
    pre_candidates = solve_two_seeds_v7510(
        "combined_pre", states.combined_pre, REF51.shape)
    pre = select_lowest_v7510(pre_candidates)
    candidate_results = collect(pre_candidates)
    if pre === nothing
        for (label, result) in candidate_results
            write_result_v51(label, result)
        end
        CSV.write(joinpath(V51_OUT, "branch_summary.csv"),
            DataFrame([summary_row_v51(label, result)
                for (label, result) in candidate_results]))
        CSV.write(joinpath(V51_OUT, "selected_states.csv"), DataFrame([
            selected_row_v7510("combined_pre", pre),
            selected_row_v7510("combined_post", nothing),
        ]))
        write(joinpath(V51_OUT, "status.txt"),
            "combined_pre_mesh_or_solver_validation_failed\n")
        return
    end

    post_candidates = (
        ("combined_post_continuation",
            solve_state_v51(states.combined_post, pre.shape;
                branch="combined_pre_continuation")),
        ("combined_post_dome",
            solve_state_v51(states.combined_post, domed_seed();
                branch="fig7_domed_seed")),
    )
    post = select_lowest_v7510(post_candidates)
    append!(candidate_results, collect(post_candidates))
    for (label, result) in candidate_results
        write_result_v51(label, result)
    end
    CSV.write(joinpath(V51_OUT, "branch_summary.csv"),
        DataFrame([summary_row_v51(label, result)
            for (label, result) in candidate_results]))
    CSV.write(joinpath(V51_OUT, "selected_states.csv"), DataFrame([
        selected_row_v7510("combined_pre", pre),
        selected_row_v7510("combined_post", post),
    ]))
    if post === nothing
        write(joinpath(V51_OUT, "status.txt"),
            "combined_post_mesh_or_solver_validation_failed\n")
        return
    end

    write_validation_v6(post)
    pre_validation = validation_metrics_v6(pre)
    post_validation = validation_metrics_v6(post)
    pre_error = max(abs(pre.metrics.apex_displacement),
        C.radius * pre.metrics.radial_cv)
    corrected_dome_increase =
        post_validation.lo_protrusion_from_ld_fit -
        pre_validation.lo_protrusion_from_ld_fit
    corrected_dome_signed_error =
        corrected_dome_increase - EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM
    corrected_dome_abs_error = abs(corrected_dome_signed_error)
    reference_ok =
        pre_error <= PRE_SPHERE_TOL_UM &&
        corrected_dome_abs_error <= CORRECTED_DOME_ACCEPTANCE_TOL_UM
    CSV.write(joinpath(V51_OUT, "reference_validation.csv"), DataFrame(
        pre_sphere_error_um=[pre_error],
        pre_sphere_tolerance_um=[PRE_SPHERE_TOL_UM],
        pre_lo_protrusion_from_ld_sphere_um=[
            pre_validation.lo_protrusion_from_ld_fit],
        post_lo_protrusion_from_ld_sphere_um=[
            post_validation.lo_protrusion_from_ld_fit],
        corrected_dome_increase_um=[corrected_dome_increase],
        experimental_corrected_dome_increase_um=[
            EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM],
        corrected_dome_signed_error_um=[corrected_dome_signed_error],
        corrected_dome_absolute_error_um=[corrected_dome_abs_error],
        accepted_comparison_tolerance_um=[
            CORRECTED_DOME_ACCEPTANCE_TOL_UM],
        global_apex_displacement_um=[post.metrics.apex_displacement],
        raw_cap_height_um=[post.metrics.cap_height],
        cap_diameter_um=[post.metrics.cap_diameter],
        reference_accepted=[reference_ok],
        interpretation=[reference_ok ?
            "10 nM reference accepted using the domain-relative Lo protrusion change" :
            "10 nM reference outside the accepted domain-relative protrusion tolerance"],
    ))
    CSV.write(joinpath(V51_OUT, "interface_balance.csv"), DataFrame(
        state=["combined_pre", "combined_post"],
        phase_mapping=[
            pre.metrics.interface.phase_mapping,
            post.metrics.interface.phase_mapping,
        ],
        two_h_lo_boundary_um_inv=[
            pre.metrics.interface.two_h_lo_boundary,
            post.metrics.interface.two_h_lo_boundary,
        ],
        two_h_ld_boundary_um_inv=[
            pre.metrics.interface.two_h_ld_boundary,
            post.metrics.interface.two_h_ld_boundary,
        ],
        lo_effective_moment_kbt_um=[
            pre.metrics.interface.moment_lo,
            post.metrics.interface.moment_lo,
        ],
        ld_effective_moment_kbt_um=[
            pre.metrics.interface.moment_ld,
            post.metrics.interface.moment_ld,
        ],
        gaussian_moment_kbt_um=[
            pre.metrics.interface.gaussian_moment,
            post.metrics.interface.gaussian_moment,
        ],
        moment_residual_kbt_um=[
            pre.metrics.interface.moment_residual,
            post.metrics.interface.moment_residual,
        ],
        moment_relative_residual=[
            pre.metrics.interface.moment_relative_residual,
            post.metrics.interface.moment_relative_residual,
        ],
        transverse_residual_kbt_um2=[
            pre.metrics.interface.transverse_residual,
            post.metrics.interface.transverse_residual,
        ],
        transverse_relative_residual=[
            pre.metrics.interface.transverse_relative_residual,
            post.metrics.interface.transverse_relative_residual,
        ],
        interpretation=fill(
            "diagnostic only; convergence must be assessed by mesh refinement",
            2),
    ))
    write(joinpath(V51_OUT, "status.txt"),
        reference_ok ? "reference_validation_passed\n" :
        "valid_mesh_but_reference_validation_failed\n")
end

function main_v770()
    V51_MODE in ("audit", "smoke") ||
        error("FIG7_V770_MODE must be audit or smoke.")
    isempty(V770_CONCENTRATIONS_NM) &&
        error("At least one streptavidin concentration is required.")
    issorted(V770_CONCENTRATIONS_NM) ||
        error("Streptavidin concentrations must be sorted in ascending order.")
    length(unique(V770_CONCENTRATIONS_NM)) == length(V770_CONCENTRATIONS_NM) ||
        error("Streptavidin concentrations must be unique.")
    first(V770_CONCENTRATIONS_NM) == 10.0 ||
        error("The sweep must start at the validated 10 nM post-CTB reference.")
    all(c -> 10.0 <= c <= V7510_STREP_SATURATION_NM,
        V770_CONCENTRATIONS_NM) ||
        error("Concentrations must lie between 10 nM and the audited saturation endpoint.")
    N_NODES >= 81 ||
        error("V7.7.0 requires at least 81 meridional nodes.")
    0.0 < INNER_LD_BIOTIN_MOLE_FRACTION < 1.0 ||
        error("The inner-Ld biotin mole fraction must lie between zero and one.")
    BOUND_SITE_AREA_COUPLING >= 0 ||
        error("The bound-site area coupling must be nonnegative.")
    V770_MONOTONIC_TOL_UM >= 0 ||
        error("The monotonicity diagnostic tolerance must be nonnegative.")
    V770_INTERFACE_MOMENT_REL_TOL > 0 ||
        error("The interface moment-residual tolerance must be positive.")
    V770_INTERFACE_TRANSVERSE_REL_TOL > 0 ||
        error("The interface transverse-residual tolerance must be positive.")
    isapprox(C0_LD_COMBINED10_UM_INV,
        TARGET_C0_LD_COMBINED10_UM_INV; atol=1e-12) ||
        error("The decomposed 10 nM Ld curvature does not preserve its target.")
    isapprox(C0_LO_POST_UM_INV,
        TARGET_C0_LO_POST_UM_INV; atol=1e-12) ||
        error("The fixed post-CTB Lo curvature does not preserve its target.")

    ensure_empty_v51_output()
    write_audit_v770()
    if V51_MODE == "audit"
        write(joinpath(V51_OUT, "status.txt"),
            "audit_only: no minimization was run.\n")
        return
    end

    states = states_v7510()
    pre_candidates = collect(solve_two_seeds_v7510(
        "pre_ctb_10nM", states.combined_pre, REF51.shape))
    pre = select_lowest_v7510(pre_candidates)
    candidate_rows = NamedTuple[]
    for (label, result) in pre_candidates
        write_result_v51(label, result)
        push!(candidate_rows, merge(
            (strep_concentration_nm=10.0, ctb_present=false),
            summary_row_v51(label, result)))
    end
    if pre === nothing
        CSV.write(joinpath(V51_OUT, "candidate_branch_summary.csv"),
            DataFrame(candidate_rows))
        write(joinpath(V51_OUT, "status.txt"),
            "pre_ctb_reference_failed_solver_or_geometry_validation\n")
        return
    end

    CSV.write(joinpath(V51_OUT, "selected_reference.csv"),
        DataFrame([selected_reference_row_v770(pre)]))
    pre_fit = ld_sphere_fit_v6(pre.shape)
    selected_results = Dict{Float64,Any}()
    previous_shape = pre.shape

    for concentration_nm in V770_CONCENTRATIONS_NM
        state = post_state_v770(concentration_nm)
        candidates = concentration_candidates_v770(
            concentration_nm, state, previous_shape, pre.shape)
        for (label, result) in candidates
            write_result_v51(label, result)
            push!(candidate_rows, merge(
                (strep_concentration_nm=concentration_nm, ctb_present=true),
                summary_row_v51(label, result)))
        end
        selected = select_lowest_v7510(candidates)
        selected_results[concentration_nm] = selected
        selected === nothing || (previous_shape = selected.shape)
    end

    CSV.write(joinpath(V51_OUT, "candidate_branch_summary.csv"),
        DataFrame(candidate_rows))
    all_selected = all(selected_results[c] !== nothing
        for c in V770_CONCENTRATIONS_NM)
    if !all_selected
        failed = [c for c in V770_CONCENTRATIONS_NM
            if selected_results[c] === nothing]
        CSV.write(joinpath(V51_OUT, "failed_concentrations.csv"),
            DataFrame(strep_concentration_nm=failed))
        write(joinpath(V51_OUT, "status.txt"),
            "one_or_more_concentrations_failed_solver_or_geometry_validation\n")
        return
    end

    ten_result = selected_results[10.0]
    ten_fit = ld_sphere_fit_v6(ten_result.shape)
    selected_rows = [
        selected_trace_row_v770(
            c, selected_results[c], pre_fit.protrusion, ten_fit.protrusion)
        for c in V770_CONCENTRATIONS_NM
    ]
    CSV.write(joinpath(V51_OUT, "selected_concentration_trace.csv"),
        DataFrame(selected_rows))

    interface_rows = NamedTuple[
        selected_interface_row_v770("pre_ctb_strep_10nM", 10.0, pre),
    ]
    append!(interface_rows, [
        selected_interface_row_v770(
            "post_ctb_strep_$(concentration_slug_v770(c))nM",
            c,
            selected_results[c])
        for c in V770_CONCENTRATIONS_NM
    ])
    CSV.write(joinpath(V51_OUT, "selected_interface_balance.csv"),
        DataFrame(interface_rows))

    write_validation_v6(ten_result)
    pre_sphere_error = max(
        abs(pre.metrics.apex_displacement),
        C.radius * pre.metrics.radial_cv)
    corrected_10nm_increase = ten_fit.protrusion - pre_fit.protrusion
    corrected_10nm_error =
        corrected_10nm_increase - EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM
    reference_ok =
        pre_sphere_error <= PRE_SPHERE_TOL_UM &&
        abs(corrected_10nm_error) <= CORRECTED_DOME_ACCEPTANCE_TOL_UM

    protrusions = [
        ld_sphere_fit_v6(selected_results[c].shape).protrusion
        for c in V770_CONCENTRATIONS_NM
    ]
    successive_changes = diff(protrusions)
    monotonic_suppression =
        all(successive_changes .<= V770_MONOTONIC_TOL_UM)
    mapping_ok = all(begin
        epsilon = incremental_inner_ld_area_fraction_v770(c)
        metrics = selected_results[c].metrics
        isapprox(metrics.c0_ld_um_inv, C0_LD_COMBINED10_UM_INV;
            atol=1e-12) &&
        isapprox(metrics.area_pref_ld,
            REF51.area_ld * (1 + epsilon / 2); atol=1e-10) &&
        isapprox(metrics.delta_pref_ld,
            REF51.delta_ld - epsilon * REF51.area_ld; atol=1e-10)
    end for c in V770_CONCENTRATIONS_NM)
    fixed_lo_ok = all(
        isapprox(
            selected_results[c].metrics.c0_lo_um_inv,
            C0_LO_POST_UM_INV;
            atol=1e-12)
        for c in V770_CONCENTRATIONS_NM)
    endpoint_corrected_protrusion =
        last(protrusions) - pre_fit.protrusion
    endpoint_restored =
        abs(endpoint_corrected_protrusion) <= PRE_SPHERE_TOL_UM
    maximum_interface_moment_relative_residual = maximum(
        selected_results[c].metrics.interface.moment_relative_residual
        for c in V770_CONCENTRATIONS_NM)
    maximum_interface_transverse_relative_residual = maximum(
        selected_results[c].metrics.interface.transverse_relative_residual
        for c in V770_CONCENTRATIONS_NM)
    interface_balance_ok =
        maximum_interface_moment_relative_residual <=
            V770_INTERFACE_MOMENT_REL_TOL &&
        maximum_interface_transverse_relative_residual <=
            V770_INTERFACE_TRANSVERSE_REL_TOL
    smoke_passed =
        reference_ok && monotonic_suppression && mapping_ok && fixed_lo_ok &&
        endpoint_restored && interface_balance_ok

    CSV.write(joinpath(V51_OUT, "reference_validation.csv"), DataFrame(
        pre_sphere_error_um=[pre_sphere_error],
        pre_sphere_tolerance_um=[PRE_SPHERE_TOL_UM],
        pre_lo_protrusion_from_ld_sphere_um=[pre_fit.protrusion],
        post_10nm_lo_protrusion_from_ld_sphere_um=[ten_fit.protrusion],
        corrected_10nm_dome_increase_um=[corrected_10nm_increase],
        experimental_corrected_dome_increase_um=[
            EXPERIMENTAL_CORRECTED_DOME_INCREASE_UM],
        corrected_10nm_signed_error_um=[corrected_10nm_error],
        corrected_10nm_absolute_error_um=[abs(corrected_10nm_error)],
        accepted_comparison_tolerance_um=[
            CORRECTED_DOME_ACCEPTANCE_TOL_UM],
        reference_accepted=[reference_ok],
    ))
    CSV.write(joinpath(V51_OUT, "smoke_validation.csv"), DataFrame(
        n_nodes=[N_NODES],
        concentrations_nm=[join(V770_CONCENTRATIONS_NM, ";")],
        all_concentrations_accepted=[all_selected],
        concentration_mapping_verified=[mapping_ok],
        fixed_lo_ctb_curvature_verified=[fixed_lo_ok],
        ten_nm_reference_accepted=[reference_ok],
        monotonic_tolerance_um=[V770_MONOTONIC_TOL_UM],
        maximum_successive_protrusion_increase_um=[
            isempty(successive_changes) ? 0.0 : maximum(successive_changes)],
        monotonic_suppression=[monotonic_suppression],
        endpoint_corrected_protrusion_from_pre_um=[
            endpoint_corrected_protrusion],
        endpoint_restored_within_pre_tolerance=[endpoint_restored],
        maximum_interface_moment_relative_residual=[
            maximum_interface_moment_relative_residual],
        interface_moment_relative_tolerance=[
            V770_INTERFACE_MOMENT_REL_TOL],
        maximum_interface_transverse_relative_residual=[
            maximum_interface_transverse_relative_residual],
        interface_transverse_relative_tolerance=[
            V770_INTERFACE_TRANSVERSE_REL_TOL],
        interface_balance_within_tolerance=[interface_balance_ok],
        smoke_passed=[smoke_passed],
        interpretation=[smoke_passed ?
            "mechanically valid smoke; production workflow may be prepared" :
            "do not proceed to production: restoration or interface-balance validation failed"],
    ))
    write(joinpath(V51_OUT, "status.txt"),
        smoke_passed ?
        "smoke_passed_monotonic_suppression\n" :
        "smoke_completed_but_not_production_ready\n")
end

main_v770()
