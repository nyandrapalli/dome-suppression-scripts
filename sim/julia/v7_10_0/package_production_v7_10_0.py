from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource, to_rgba
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


LD_COLOR = "#168447"
LO_COLOR = "#d9362b"
BACKGROUND = "#f5f6f8"
FULL_AZIMUTH = 144
RENDER_AZIMUTH = 72
VIDEO_FPS = 12
EQUILIBRIUM_HOLD_FRAMES = 12


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_profile(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = read_csv(path)
    if len(rows) < 3:
        raise ValueError(f"Profile contains fewer than three nodes: {path}")
    rho = np.array([float(row["rho_um"]) for row in rows], dtype=float)
    z = np.array([float(row["z_um"]) for row in rows], dtype=float)
    domains = np.array([row["domain"] for row in rows])
    if not np.isfinite(rho).all() or not np.isfinite(z).all():
        raise ValueError(f"Profile contains nonfinite coordinates: {path}")
    return rho, z, domains


def mesh_from_profile(
    rho: np.ndarray,
    z: np.ndarray,
    node_domains: np.ndarray,
    n_azimuth: int,
) -> tuple[np.ndarray, list[list[int]], list[str]]:
    if len(rho) != len(z) or len(rho) != len(node_domains):
        raise ValueError("Profile arrays must have identical lengths.")
    if n_azimuth < 12:
        raise ValueError("At least 12 azimuthal subdivisions are required.")

    phi = np.linspace(0.0, 2.0 * np.pi, n_azimuth, endpoint=False)
    vertices: list[list[float]] = [[0.0, 0.0, float(z[0])]]
    for index in range(1, len(rho) - 1):
        vertices.extend(
            [
                [
                    float(rho[index] * math.cos(angle)),
                    float(rho[index] * math.sin(angle)),
                    float(z[index]),
                ]
                for angle in phi
            ]
        )
    south_index = len(vertices)
    vertices.append([0.0, 0.0, float(z[-1])])

    def ring_start(profile_index: int) -> int:
        return 1 + (profile_index - 1) * n_azimuth

    faces: list[list[int]] = []
    face_domains: list[str] = []
    first_ring = ring_start(1)
    for azimuth in range(n_azimuth):
        nxt = (azimuth + 1) % n_azimuth
        faces.append([0, first_ring + azimuth, first_ring + nxt])
        face_domains.append(str(node_domains[1]))

    for profile_index in range(1, len(rho) - 2):
        upper = ring_start(profile_index)
        lower = ring_start(profile_index + 1)
        for azimuth in range(n_azimuth):
            nxt = (azimuth + 1) % n_azimuth
            faces.append(
                [
                    upper + azimuth,
                    upper + nxt,
                    lower + nxt,
                    lower + azimuth,
                ]
            )
            face_domains.append(str(node_domains[profile_index + 1]))

    last_ring = ring_start(len(rho) - 2)
    for azimuth in range(n_azimuth):
        nxt = (azimuth + 1) % n_azimuth
        faces.append([last_ring + azimuth, last_ring + nxt, south_index])
        face_domains.append(str(node_domains[-1]))
    return np.asarray(vertices, dtype=float), faces, face_domains


def downsample_profile(
    rho: np.ndarray,
    z: np.ndarray,
    domains: np.ndarray,
    maximum_nodes: int = 181,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if len(rho) <= maximum_nodes:
        return rho, z, domains
    indices = set(
        np.linspace(0, len(rho) - 1, maximum_nodes, dtype=int).tolist()
    )
    transitions = np.flatnonzero(domains[1:] != domains[:-1])
    for transition in transitions:
        for offset in (-1, 0, 1, 2):
            candidate = int(transition + offset)
            if 0 <= candidate < len(rho):
                indices.add(candidate)
    selected = np.array(sorted(indices), dtype=int)
    return rho[selected], z[selected], domains[selected]


def write_mtl(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "newmtl Lo",
                "Kd 0.847 0.212 0.169",
                "Ka 0.120 0.030 0.025",
                "Ks 0.150 0.150 0.150",
                "Ns 40",
                "",
                "newmtl Ld",
                "Kd 0.086 0.518 0.278",
                "Ka 0.020 0.100 0.050",
                "Ks 0.150 0.150 0.150",
                "Ns 40",
                "",
            ]
        ),
        encoding="ascii",
    )


def write_obj(
    path: Path,
    vertices: np.ndarray,
    faces: list[list[int]],
    domains: list[str],
) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("# V7.10 axisymmetric equilibrium mesh; coordinates in um\n")
        handle.write("mtllib domains.mtl\n")
        for vertex in vertices:
            handle.write(
                f"v {vertex[0]:.10g} {vertex[1]:.10g} "
                f"{vertex[2]:.10g}\n"
            )
        current = None
        for face, domain in zip(faces, domains):
            if domain != current:
                handle.write(f"g {domain}\nusemtl {domain}\n")
                current = domain
            indices = " ".join(str(index + 1) for index in face)
            handle.write(f"f {indices}\n")


def render_profile(
    path: Path,
    rho: np.ndarray,
    z: np.ndarray,
    node_domains: np.ndarray,
    limit: float,
    title: str,
) -> None:
    render_rho, render_z, render_domains = downsample_profile(
        rho, z, node_domains
    )
    vertices, faces, domains = mesh_from_profile(
        render_rho, render_z, render_domains, RENDER_AZIMUTH
    )
    polygons = [[vertices[index] for index in face] for face in faces]
    colors = [
        to_rgba(LO_COLOR if domain == "Lo" else LD_COLOR)
        for domain in domains
    ]
    figure = plt.figure(figsize=(5, 5), dpi=180, facecolor="white")
    axis = figure.add_subplot(111, projection="3d")
    axis.add_collection3d(
        Poly3DCollection(
            polygons,
            facecolors=colors,
            edgecolors=colors,
            linewidths=0.0,
            antialiased=False,
            shade=True,
            lightsource=LightSource(azdeg=315, altdeg=45),
        )
    )
    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_zlim(-limit, limit)
    axis.set_box_aspect((1, 1, 1))
    axis.view_init(elev=18, azim=-58, roll=-90)
    axis.set_axis_off()
    axis.set_facecolor(BACKGROUND)
    for coordinate_axis in (axis.xaxis, axis.yaxis, axis.zaxis):
        coordinate_axis.pane.set_facecolor(BACKGROUND)
        coordinate_axis.pane.set_edgecolor(BACKGROUND)
        coordinate_axis.pane.set_alpha(1.0)
    axis.set_title(title, fontsize=11, pad=0)
    figure.subplots_adjust(left=0, right=1, bottom=0, top=0.94)
    figure.savefig(path, facecolor="white")
    plt.close(figure)


def encode_video(frame_directory: Path, output: Path) -> None:
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        raise RuntimeError("ffmpeg is required but was not found on PATH.")
    subprocess.run(
        [
            ffmpeg,
            "-hide_banner",
            "-loglevel",
            "error",
            "-y",
            "-framerate",
            str(VIDEO_FPS),
            "-i",
            str(frame_directory / "frame_%05d.png"),
            "-c:v",
            "libx264",
            "-crf",
            "18",
            "-pix_fmt",
            "yuv420p",
            str(output),
        ],
        check=True,
    )


def selected_pre_profile(run: Path) -> tuple[Path, dict[str, str]]:
    selected = read_csv(run / "selected_reference.csv")
    if len(selected) != 1:
        raise ValueError("Expected one selected pre-CTB reference.")
    row = selected[0]
    names = {
        "near_sphere": "pre_ctb_10nM_sphere_nodes.csv",
        "fig7_domed_seed": "pre_ctb_10nM_dome_nodes.csv",
    }
    if row["branch"] not in names:
        raise ValueError(f"Unknown selected pre-CTB branch: {row['branch']}")
    return run / names[row["branch"]], row


def equilibrium_states(
    run: Path,
) -> list[dict[str, object]]:
    pre_path, pre_metrics = selected_pre_profile(run)
    states: list[dict[str, object]] = [
        {
            "tag": "pre_ctb_10nM",
            "title": "Pre-CTB, 10 nM streptavidin",
            "source": pre_path,
            "metrics": pre_metrics,
            "concentration_nm": 10.0,
            "ctb_present": False,
            "branch": pre_metrics["branch"],
        }
    ]
    for row in read_csv(run / "selected_concentration_trace.csv"):
        concentration = float(row["strep_concentration_nm"])
        concentration_tag = f"{concentration:g}".replace(".", "p")
        states.append(
            {
                "tag": f"post_ctb_strep_{concentration_tag}nM",
                "title": f"Post-CTB, streptavidin = {concentration:g} nM",
                "source": run / f"{row['branch']}_nodes.csv",
                "metrics": row,
                "concentration_nm": concentration,
                "ctb_present": True,
                "branch": row["branch"],
            }
        )
    return states


def make_equilibrium_video(
    png_paths: list[Path],
    frame_directory: Path,
    output: Path,
) -> int:
    frame_directory.mkdir(parents=True)
    serial = 0
    for state_index, png in enumerate(png_paths):
        repeats = EQUILIBRIUM_HOLD_FRAMES
        if state_index in (0, len(png_paths) - 1):
            repeats *= 2
        for _ in range(repeats):
            shutil.copy2(png, frame_directory / f"frame_{serial:05d}.png")
            serial += 1
    encode_video(frame_directory, output)
    return serial


def optimizer_frames(
    run: Path,
    concentration_states: list[dict[str, object]],
    output: Path,
    limit: float,
) -> tuple[int, list[dict[str, str]]]:
    frame_directory = output / "numerical_equilibration_frames"
    frame_directory.mkdir(parents=True)
    aggregate_trace: list[dict[str, str]] = []
    serial = 0
    for state in concentration_states:
        branch = str(state["branch"])
        concentration = float(state["concentration_nm"])
        callback = run / "optimizer_callbacks" / branch
        manifest_path = callback / "optimizer_frame_manifest.csv"
        trace_path = callback / "optimizer_trace.csv"
        if not manifest_path.exists() or not trace_path.exists():
            raise FileNotFoundError(
                f"Missing optimizer callback output for selected branch {branch}"
            )
        for row in read_csv(trace_path):
            aggregate_trace.append(
                {
                    "strep_concentration_nm": f"{concentration:g}",
                    **row,
                }
            )
        manifest = read_csv(manifest_path)
        if not manifest:
            raise ValueError(f"No optimizer frames captured for {branch}")
        for row in manifest:
            rho, z, domains = read_profile(callback / row["nodes_file"])
            render_profile(
                frame_directory / f"frame_{serial:05d}.png",
                rho,
                z,
                domains,
                limit,
                (
                    f"Numerical minimization: {concentration:g} nM, "
                    f"iteration {int(row['optimizer_iteration'])}"
                ),
            )
            serial += 1
    if not aggregate_trace:
        raise ValueError("Selected optimizer traces are empty.")
    with (output / "selected_optimizer_trace.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(aggregate_trace[0])
        )
        writer.writeheader()
        writer.writerows(aggregate_trace)
    encode_video(frame_directory, output / "numerical_equilibration.mp4")
    return serial, aggregate_trace


def make_graphs(
    trace: list[dict[str, str]],
    output: Path,
) -> None:
    concentration = np.array(
        [float(row["strep_concentration_nm"]) for row in trace]
    )
    corrected = np.array(
        [float(row["corrected_lo_protrusion_from_pre_um"]) for row in trace]
    )
    raw = np.array(
        [float(row["lo_protrusion_from_ld_sphere_um"]) for row in trace]
    )

    figure, axis = plt.subplots(figsize=(5.2, 3.8), dpi=220)
    axis.plot(concentration, raw, "o-", color="#d9362b", label="Lo protrusion")
    axis.plot(
        concentration,
        corrected,
        "s-",
        color="#1767a0",
        label="Corrected from pre-CTB",
    )
    axis.axhline(0.25, color="#555555", linestyle="--", linewidth=1)
    axis.set_xlabel("Streptavidin concentration (nM)")
    axis.set_ylabel("Lo protrusion (µm)")
    axis.legend(frameon=False)
    axis.grid(color="#d9dde1", linewidth=0.6)
    figure.tight_layout()
    figure.savefig(output / "concentration_response.png")
    plt.close(figure)

    components = {
        "Lo bending": "lo_bending_energy_kbt",
        "Ld bending": "ld_bending_energy_kbt",
        "Lo mismatch": "lo_mismatch_energy_kbt",
        "Ld mismatch": "ld_mismatch_energy_kbt",
        "Total": "total_energy_kbt",
    }
    figure, axis = plt.subplots(figsize=(5.2, 3.8), dpi=220)
    for label, field in components.items():
        values = np.array([float(row[field]) for row in trace])
        axis.plot(concentration, values, "o-", label=label)
    axis.set_yscale("log")
    axis.set_xlabel("Streptavidin concentration (nM)")
    axis.set_ylabel("Energy ($k_\\mathrm{B}T$)")
    axis.legend(frameon=False, fontsize=8)
    axis.grid(color="#d9dde1", linewidth=0.6, which="both")
    figure.tight_layout()
    figure.savefig(output / "energy_partition.png")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run = args.run.resolve()
    output = args.output.resolve()

    if not run.is_dir():
        raise FileNotFoundError(run)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"Refusing to overwrite nonempty media directory: {output}"
        )
    output.mkdir(parents=True, exist_ok=True)

    states = equilibrium_states(run)
    profiles = []
    for state in states:
        rho, z, domains = read_profile(Path(state["source"]))
        profiles.append((rho, z, domains))
    limit = 1.08 * max(
        max(float(np.max(np.abs(rho))), float(np.max(np.abs(z))))
        for rho, z, _ in profiles
    )

    obj_directory = output / "equilibrium_obj"
    png_directory = output / "equilibrium_png"
    obj_directory.mkdir()
    png_directory.mkdir()
    write_mtl(obj_directory / "domains.mtl")

    _, faces, face_domains = mesh_from_profile(
        *profiles[0], FULL_AZIMUTH
    )
    json_frames = []
    png_paths = []
    for state, profile in zip(states, profiles):
        rho, z, domains = profile
        vertices, current_faces, current_domains = mesh_from_profile(
            rho, z, domains, FULL_AZIMUTH
        )
        if current_faces != faces or current_domains != face_domains:
            raise ValueError("Mesh topology or domain assignment changed.")
        tag = str(state["tag"])
        obj_path = obj_directory / f"{tag}.obj"
        png_path = png_directory / f"{tag}.png"
        write_obj(obj_path, vertices, faces, face_domains)
        render_profile(
            png_path,
            rho,
            z,
            domains,
            limit,
            str(state["title"]),
        )
        png_paths.append(png_path)
        json_frames.append(
            {
                "tag": tag,
                "title": state["title"],
                "strep_concentration_nm": state["concentration_nm"],
                "ctb_present": state["ctb_present"],
                "branch": state["branch"],
                "metrics": state["metrics"],
                "vertices_um": vertices.round(9).tolist(),
                "source_nodes_csv": str(Path(state["source"]).resolve()),
            }
        )

    trajectory = {
        "schema": "fig7-v7.10-axisymmetric-production-trajectory-v1",
        "units": {"length": "um", "energy": "kBT"},
        "color_rule": {"Ld": LD_COLOR, "Lo": LO_COLOR},
        "warnings": [
            "Concentration frames are separately minimized equilibria, not physical-time dynamics.",
            "Optimizer frames are numerical Ipopt iterates, not physical-time dynamics.",
            "The 3D meshes are surfaces of revolution of the axisymmetric solver profile.",
            "The concentration-to-occupancy mapping is a declared model assumption.",
        ],
        "topology": {
            "azimuthal_subdivisions": FULL_AZIMUTH,
            "faces_zero_based": faces,
            "face_domain": face_domains,
        },
        "frames": json_frames,
    }
    with (output / "trajectory.json").open("w", encoding="utf-8") as handle:
        json.dump(trajectory, handle, separators=(",", ":"))

    equilibrium_count = make_equilibrium_video(
        png_paths,
        output / "concentration_equilibrium_frames",
        output / "concentration_equilibria.mp4",
    )
    concentration_states = [
        state for state in states if bool(state["ctb_present"])
    ]
    optimizer_count, optimizer_trace = optimizer_frames(
        run, concentration_states, output, limit
    )

    selected_trace = read_csv(run / "selected_concentration_trace.csv")
    shutil.copy2(
        run / "selected_concentration_trace.csv",
        output / "production_metrics.csv",
    )
    make_graphs(selected_trace, output)

    manifest = {
        "schema": "fig7-v7.10-production-manifest-v1",
        "raw_run": str(run),
        "media_output": str(output),
        "equilibrium_state_count": len(states),
        "equilibrium_obj_count": len(list(obj_directory.glob("*.obj"))),
        "equilibrium_png_count": len(list(png_directory.glob("*.png"))),
        "equilibrium_video_frame_count": equilibrium_count,
        "optimizer_video_frame_count": optimizer_count,
        "optimizer_trace_row_count": len(optimizer_trace),
        "full_azimuthal_subdivisions": FULL_AZIMUTH,
        "render_azimuthal_subdivisions": RENDER_AZIMUTH,
        "colors": {"Ld": LD_COLOR, "Lo": LO_COLOR},
        "videos": {
            "concentration_equilibria": "concentration_equilibria.mp4",
            "numerical_equilibration": "numerical_equilibration.mp4",
        },
    }
    (output / "media_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    (output / "status.txt").write_text("complete\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
