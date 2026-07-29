using CSV
using DataFrames

const V710_BASE_SOURCE = normpath(joinpath(
    @__DIR__, "..", "v7_7_0", "run_fig7_local_leaflet_interface_v7_7_0.jl"))

function load_v770_library()
    source = read(V710_BASE_SOURCE, String)
    library = replace(source, r"\nmain_v770\(\)\s*$" => "\n")
    library == source &&
        error("Could not remove the V7.7.0 main entry point.")
    Base.include_string(Main, library, V710_BASE_SOURCE)
end

load_v770_library()

const V710_CALLBACK_STRIDE = parse(
    Int, get(ENV, "FIG7_V710_CALLBACK_STRIDE", "10"))
const V710_MAX_CALLBACK_FRAMES = parse(
    Int, get(ENV, "FIG7_V710_MAX_CALLBACK_FRAMES", "12"))
const V710_CAPTURE_CALLBACK = lowercase(
    get(ENV, "FIG7_V710_CAPTURE_CALLBACK", "true")) == "true"

function optimizer_output_directory_v710(branch)
    safe = replace(branch, r"[^A-Za-z0-9_.-]" => "_")
    joinpath(V51_OUT, "optimizer_callbacks", safe)
end

function write_optimizer_capture_v710(
        branch, state, optimizer_trace, optimizer_frames, term, primal)
    output = optimizer_output_directory_v710(branch)
    mkpath(output)
    if !isempty(optimizer_trace)
        CSV.write(
            joinpath(output, "optimizer_trace.csv"),
            DataFrame(optimizer_trace),
        )
    end

    frames = joinpath(output, "frames")
    mkpath(frames)
    manifest = NamedTuple[]
    for (serial, frame) in enumerate(optimizer_frames)
        basename = "frame_$(lpad(string(serial - 1), 5, '0'))"
        relative_path = "frames/$(basename)_nodes.csv"
        CSV.write(joinpath(output, relative_path), DataFrame(
            node=1:N_NODES,
            rho_um=frame.rho,
            z_um=frame.z,
            psi_rad=frame.psi,
            domain=[i <= DB_NODE ? "Lo" : "Ld" for i in 1:N_NODES],
        ))
        push!(manifest, (
            serial_frame=serial - 1,
            optimizer_iteration=frame.iteration,
            nodes_file=relative_path,
        ))
    end
    CSV.write(
        joinpath(output, "optimizer_frame_manifest.csv"),
        DataFrame(manifest),
    )
    CSV.write(joinpath(output, "capture_metadata.csv"), DataFrame(
        state=[state.label],
        branch=[branch],
        n_nodes=[N_NODES],
        callback_stride=[V710_CALLBACK_STRIDE],
        maximum_callback_frames=[V710_MAX_CALLBACK_FRAMES],
        captured_frames=[length(optimizer_frames)],
        captured_iterations=[length(optimizer_trace)],
        termination_status=[string(term)],
        primal_status=[string(primal)],
        interpretation=[
            "Numerical Ipopt iterates; these frames are not physical time.",
        ],
    ))
end

# This method intentionally replaces the V7.7 implementation with the same
# objective and constraints plus read-only Ipopt callback capture.
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
    @variable(
        model,
        0.25 * REF51.shape.ell_lo <= ell_lo <=
            4.0 * REF51.shape.ell_lo,
    )
    @variable(
        model,
        0.25 * REF51.shape.ell_ld <= ell_ld <=
            4.0 * REF51.shape.ell_ld,
    )

    start = seed isa NamedTuple && hasproperty(seed, :ell_lo) ?
        seed : start_shape_v51(seed)
    for i in 1:N_NODES
        set_start_value(rho[i], start.rho[i])
        set_start_value(z[i], start.z[i])
        set_start_value(psi[i], start.psi[i])
    end
    set_start_value(ell_lo, start.ell_lo)
    set_start_value(ell_ld, start.ell_ld)

    fix(rho[1], 0.0; force=true)
    fix(rho[end], 0.0; force=true)
    fix(psi[1], 0.0; force=true)
    fix(psi[end], pi; force=true)
    @constraint(model, z[1] + z[end] == 0.0)
    for i in 1:N_SEGMENTS
        local_length = i <= N_LO_SEGMENTS ? ell_lo : ell_ld
        @NLconstraint(
            model,
            rho[i + 1] - rho[i] ==
                local_length * cos((psi[i] + psi[i + 1]) / 2),
        )
        @NLconstraint(
            model,
            z[i + 1] - z[i] ==
                -local_length * sin((psi[i] + psi[i + 1]) / 2),
        )
    end

    @NLexpression(
        model,
        area_lo,
        sum(
            pi * (rho[i] + rho[i + 1]) * ell_lo
            for i in 1:N_LO_SEGMENTS
        ),
    )
    @NLexpression(
        model,
        area_ld,
        sum(
            pi * (rho[i] + rho[i + 1]) * ell_ld
            for i in (N_LO_SEGMENTS + 1):N_SEGMENTS
        ),
    )
    @NLexpression(
        model,
        volume,
        sum(
            pi * (z[i] - z[i + 1]) *
            (rho[i]^2 + rho[i] * rho[i + 1] + rho[i + 1]^2) / 3
            for i in 1:N_SEGMENTS
        ),
    )
    @NLconstraint(model, volume == REF51.volume)

    @NLexpression(
        model,
        ks_lo[i=1:N_LO_SEGMENTS],
        (psi[i + 1] - psi[i]) / ell_lo,
    )
    @NLexpression(
        model,
        kp_lo[i=1:N_LO_SEGMENTS],
        sin((psi[i] + psi[i + 1]) / 2) /
            ((rho[i] + rho[i + 1]) / 2 + 1e-10),
    )
    @NLexpression(
        model,
        tm_lo[i=1:N_LO_SEGMENTS],
        ks_lo[i] + kp_lo[i],
    )
    @NLexpression(
        model,
        ks_ld[i=(N_LO_SEGMENTS + 1):N_SEGMENTS],
        (psi[i + 1] - psi[i]) / ell_ld,
    )
    @NLexpression(
        model,
        kp_ld[i=(N_LO_SEGMENTS + 1):N_SEGMENTS],
        sin((psi[i] + psi[i + 1]) / 2) /
            ((rho[i] + rho[i + 1]) / 2 + 1e-10),
    )
    @NLexpression(
        model,
        tm_ld[i=(N_LO_SEGMENTS + 1):N_SEGMENTS],
        ks_ld[i] + kp_ld[i],
    )
    @NLexpression(
        model,
        delta_lo,
        2 * NEUTRAL_OFFSET_UM * sum(
            tm_lo[i] * pi * (rho[i] + rho[i + 1]) * ell_lo
            for i in 1:N_LO_SEGMENTS
        ),
    )
    @NLexpression(
        model,
        delta_ld,
        2 * NEUTRAL_OFFSET_UM * sum(
            tm_ld[i] * pi * (rho[i] + rho[i + 1]) * ell_ld
            for i in (N_LO_SEGMENTS + 1):N_SEGMENTS
        ),
    )

    area_pref_lo = REF51.area_lo * (1 + state.lo_area_expansion)
    area_pref_ld = REF51.area_ld * (1 + state.ld_area_expansion)
    delta_pref_lo =
        REF51.delta_lo + state.lo_mismatch * REF51.area_lo
    delta_pref_ld =
        REF51.delta_ld + state.ld_mismatch * REF51.area_ld
    @NLexpression(
        model,
        bend_lo,
        0.5 * KAPPA_LO_KBT * sum(
            (tm_lo[i] - state.c0_lo_um_inv)^2 *
            pi * (rho[i] + rho[i + 1]) * ell_lo
            for i in 1:N_LO_SEGMENTS
        ),
    )
    @NLexpression(
        model,
        bend_ld,
        0.5 * KAPPA_LD_KBT * sum(
            (tm_ld[i] - state.c0_ld_um_inv)^2 *
            pi * (rho[i] + rho[i + 1]) * ell_ld
            for i in (N_LO_SEGMENTS + 1):N_SEGMENTS
        ),
    )
    @NLexpression(
        model,
        ealo,
        KA_LO * (area_lo - area_pref_lo)^2 / (2 * REF51.area_lo),
    )
    @NLexpression(
        model,
        eald,
        KA_LD * (area_ld - area_pref_ld)^2 / (2 * REF51.area_ld),
    )
    @NLexpression(
        model,
        edlo,
        KDELTA_LO * (delta_lo - delta_pref_lo)^2 /
            (2 * REF51.area_lo),
    )
    @NLexpression(
        model,
        edld,
        KDELTA_LD * (delta_ld - delta_pref_ld)^2 /
            (2 * REF51.area_ld),
    )
    @NLexpression(
        model,
        eline,
        2pi * LINE_TENSION_KBT_UM * rho[DB_NODE],
    )
    @NLexpression(
        model,
        egauss,
        2pi * DELTA_KAPPA_G_KBT * cos(psi[DB_NODE]),
    )
    @NLobjective(
        model,
        Min,
        bend_lo + bend_ld + ealo + eald + edlo + edld + eline + egauss,
    )

    optimizer_trace = NamedTuple[]
    optimizer_frames = NamedTuple[]
    if V710_CAPTURE_CALLBACK
        function v710_callback(
                alg_mod::Cint,
                iter_count::Cint,
                obj_value::Float64,
                inf_pr::Float64,
                inf_du::Float64,
                mu::Float64,
                d_norm::Float64,
                regularization_size::Float64,
                alpha_du::Float64,
                alpha_pr::Float64,
                ls_trials::Cint)
            push!(optimizer_trace, (
                state=state.label,
                branch=branch,
                iteration=Int(iter_count),
                algorithm_mode=Int(alg_mod),
                objective_kbt=obj_value,
                primal_infeasibility=inf_pr,
                dual_infeasibility=inf_du,
                barrier_parameter=mu,
                step_norm=d_norm,
                regularization_size=regularization_size,
                dual_step=alpha_du,
                primal_step=alpha_pr,
                line_search_trials=Int(ls_trials),
            ))
            if iter_count % V710_CALLBACK_STRIDE == 0 &&
                    length(optimizer_frames) < V710_MAX_CALLBACK_FRAMES
                push!(optimizer_frames, (
                    iteration=Int(iter_count),
                    rho=[
                        callback_value(model, rho[i])
                        for i in 1:N_NODES
                    ],
                    z=[
                        callback_value(model, z[i])
                        for i in 1:N_NODES
                    ],
                    psi=[
                        callback_value(model, psi[i])
                        for i in 1:N_NODES
                    ],
                ))
            end
            true
        end
        MOI.set(model, Ipopt.CallbackFunction(), v710_callback)
    end

    optimize!(model)
    term = termination_status(model)
    primal = primal_status(model)
    write_optimizer_capture_v710(
        branch,
        state,
        optimizer_trace,
        optimizer_frames,
        term,
        primal,
    )

    has_solution =
        primal in (MOI.FEASIBLE_POINT, MOI.NEARLY_FEASIBLE_POINT)
    has_solution || return (
        status=string(term),
        primal_status=string(primal),
        branch=branch,
        accepted=false,
        objective=NaN,
        solve_time=solve_time(model),
        shape=nothing,
        metrics=nothing,
    )

    shape = (
        rho=value.(rho),
        z=value.(z),
        psi=value.(psi),
        ell_lo=value(ell_lo),
        ell_ld=value(ell_ld),
    )
    metrics = compute_metrics_v51(shape, state)
    solver_ok = term in (
        MOI.LOCALLY_SOLVED,
        MOI.OPTIMAL,
        MOI.ALMOST_LOCALLY_SOLVED,
        MOI.ALMOST_OPTIMAL,
    )
    mesh_ok =
        metrics.max_geometry_error < 1e-7 &&
        abs(metrics.volume_error) < 1e-8 &&
        metrics.axial_order_violations == 0 &&
        metrics.interior_axis_contacts == 0 &&
        metrics.profile_self_intersections == 0 &&
        metrics.segment_ratio <= MAX_SEGMENT_RATIO
    (
        status=string(term),
        primal_status=string(primal),
        branch=branch,
        accepted=solver_ok && mesh_ok,
        objective=objective_value(model),
        solve_time=solve_time(model),
        shape=shape,
        metrics=metrics,
    )
end

V710_CALLBACK_STRIDE >= 1 ||
    error("FIG7_V710_CALLBACK_STRIDE must be at least one.")
V710_MAX_CALLBACK_FRAMES >= 1 ||
    error("FIG7_V710_MAX_CALLBACK_FRAMES must be at least one.")

Base.invokelatest(getfield(Main, :main_v770))
