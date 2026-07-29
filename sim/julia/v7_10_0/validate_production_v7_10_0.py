from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_config(path: Path) -> dict[str, str]:
    return {
        row["parameter"]: row["value"]
        for row in read_csv(path)
    }


def close(actual: str, expected: float, tolerance: float = 1.0e-12) -> bool:
    return math.isclose(
        float(actual), expected, rel_tol=0.0, abs_tol=tolerance
    )


def video_metadata(path: Path) -> dict[str, object]:
    ffprobe = shutil.which("ffprobe")
    if not ffprobe:
        raise RuntimeError("ffprobe is required but was not found on PATH.")
    result = subprocess.run(
        [
            ffprobe,
            "-v",
            "error",
            "-select_streams",
            "v:0",
            "-show_entries",
            "stream=width,height,duration,nb_frames",
            "-of",
            "json",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    streams = json.loads(result.stdout).get("streams", [])
    if len(streams) != 1:
        raise ValueError(f"Expected one video stream in {path}")
    return streams[0]


def count_obj(path: Path) -> tuple[int, int, set[str]]:
    vertices = 0
    faces = 0
    groups: set[str] = set()
    with path.open(encoding="ascii") as handle:
        for line in handle:
            if line.startswith("v "):
                vertices += 1
            elif line.startswith("f "):
                faces += 1
            elif line.startswith("g "):
                groups.add(line.split(maxsplit=1)[1].strip())
    return vertices, faces, groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--media", type=Path, required=True)
    parser.add_argument("--nodes", type=int, required=True)
    parser.add_argument("--concentrations", required=True)
    parser.add_argument("--require-restoration", action="store_true")
    args = parser.parse_args()
    run = args.run.resolve()
    media = args.media.resolve()
    expected_concentrations = [
        float(value) for value in args.concentrations.split(",")
    ]

    config = read_config(run / "run_config.csv")
    trace = read_csv(run / "selected_concentration_trace.csv")
    reference = read_csv(run / "reference_validation.csv")
    if len(reference) != 1:
        raise ValueError("Expected one reference-validation row.")
    concentrations = [
        float(row["strep_concentration_nm"]) for row in trace
    ]
    checks: list[tuple[str, bool, object, object]] = [
        (
            "nodes",
            config.get("n_nodes") == str(args.nodes),
            config.get("n_nodes"),
            str(args.nodes),
        ),
        (
            "gamma_A",
            close(config["bound_site_area_coupling"], 4.0),
            config["bound_site_area_coupling"],
            "4.0",
        ),
        (
            "concentrations",
            concentrations == expected_concentrations,
            ";".join(f"{value:g}" for value in concentrations),
            ";".join(f"{value:g}" for value in expected_concentrations),
        ),
        (
            "all_selected_states_accepted",
            all(row["accepted"].lower() == "true" for row in trace),
            all(row["accepted"].lower() == "true" for row in trace),
            True,
        ),
        (
            "10nM_reference_accepted",
            reference[0]["reference_accepted"].lower() == "true",
            reference[0]["reference_accepted"],
            "true",
        ),
    ]

    protrusions = [
        float(row["lo_protrusion_from_ld_sphere_um"]) for row in trace
    ]
    checks.append(
        (
            "monotonic_suppression",
            all(
                next_value <= previous + 0.02
                for previous, next_value in zip(
                    protrusions[:-1], protrusions[1:]
                )
            ),
            max(
                (
                    next_value - previous
                    for previous, next_value in zip(
                        protrusions[:-1], protrusions[1:]
                    )
                ),
                default=0.0,
            ),
            "<=0.02 um successive increase",
        )
    )

    for concentration, row in zip(concentrations, trace):
        coverage = min(concentration / 80.0, 1.0)
        expected_area = 0.002 * 4.0 * max(coverage - 0.125, 0.0)
        checks.extend(
            [
                (
                    f"{concentration:g}nM_incremental_area",
                    close(
                        row["incremental_inner_ld_area_fraction"],
                        expected_area,
                    ),
                    row["incremental_inner_ld_area_fraction"],
                    expected_area,
                ),
                (
                    f"{concentration:g}nM_C0_Ld_fixed",
                    close(row["effective_c0_ld_um_inv"], -0.015),
                    row["effective_c0_ld_um_inv"],
                    -0.015,
                ),
                (
                    f"{concentration:g}nM_C0_Lo_fixed",
                    close(
                        row["effective_c0_lo_um_inv"], 1.0 / 3.0
                    ),
                    row["effective_c0_lo_um_inv"],
                    1.0 / 3.0,
                ),
                (
                    f"{concentration:g}nM_volume",
                    abs(float(row["volume_relative_error"])) <= 1.0e-8,
                    row["volume_relative_error"],
                    "<=1e-8 absolute",
                ),
                (
                    f"{concentration:g}nM_geometry",
                    float(row["geometry_error"]) <= 1.0e-7,
                    row["geometry_error"],
                    "<=1e-7 um",
                ),
                (
                    f"{concentration:g}nM_segment_ratio",
                    float(row["segment_ratio"]) <= 2.0,
                    row["segment_ratio"],
                    "<=2.0",
                ),
            ]
        )

    endpoint = trace[-1]
    if args.require_restoration:
        checks.append(
            (
                "80nM_restored_within_0p25um",
                float(
                    endpoint["corrected_lo_protrusion_from_pre_um"]
                )
                <= 0.25,
                endpoint["corrected_lo_protrusion_from_pre_um"],
                "<=0.25 um",
            )
        )

    trajectory_path = media / "trajectory.json"
    manifest_path = media / "media_manifest.json"
    trajectory = json.loads(trajectory_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_states = len(expected_concentrations) + 1
    checks.extend(
        [
            (
                "trajectory_schema",
                trajectory.get("schema")
                == "fig7-v7.10-axisymmetric-production-trajectory-v1",
                trajectory.get("schema"),
                "fig7-v7.10-axisymmetric-production-trajectory-v1",
            ),
            (
                "color_rule",
                trajectory.get("color_rule")
                == {"Ld": "#168447", "Lo": "#d9362b"},
                trajectory.get("color_rule"),
                {"Ld": "#168447", "Lo": "#d9362b"},
            ),
            (
                "trajectory_state_count",
                len(trajectory["frames"]) == expected_states,
                len(trajectory["frames"]),
                expected_states,
            ),
            (
                "manifest_equilibrium_count",
                int(manifest["equilibrium_state_count"])
                == expected_states,
                manifest["equilibrium_state_count"],
                expected_states,
            ),
            (
                "optimizer_video_frames",
                int(manifest["optimizer_video_frame_count"]) > 0,
                manifest["optimizer_video_frame_count"],
                ">0",
            ),
            (
                "optimizer_trace_rows",
                int(manifest["optimizer_trace_row_count"]) > 0,
                manifest["optimizer_trace_row_count"],
                ">0",
            ),
        ]
    )

    obj_paths = sorted((media / "equilibrium_obj").glob("*.obj"))
    png_paths = sorted((media / "equilibrium_png").glob("*.png"))
    checks.extend(
        [
            (
                "OBJ_count",
                len(obj_paths) == expected_states,
                len(obj_paths),
                expected_states,
            ),
            (
                "PNG_count",
                len(png_paths) == expected_states,
                len(png_paths),
                expected_states,
            ),
        ]
    )
    for path in obj_paths:
        vertices, faces, groups = count_obj(path)
        checks.extend(
            [
                (
                    f"{path.stem}_OBJ_vertices",
                    vertices > 1000,
                    vertices,
                    ">1000",
                ),
                (
                    f"{path.stem}_OBJ_faces",
                    faces > 1000,
                    faces,
                    ">1000",
                ),
                (
                    f"{path.stem}_OBJ_domains",
                    groups == {"Lo", "Ld"},
                    ";".join(sorted(groups)),
                    "Ld;Lo",
                ),
            ]
        )

    for name in ("concentration_equilibria.mp4", "numerical_equilibration.mp4"):
        path = media / name
        metadata = video_metadata(path)
        duration = float(metadata.get("duration", 0.0))
        checks.extend(
            [
                (
                    f"{name}_nonempty",
                    path.stat().st_size > 10_000,
                    path.stat().st_size,
                    ">10000 bytes",
                ),
                (
                    f"{name}_duration",
                    duration > 0.5,
                    duration,
                    ">0.5 s",
                ),
                (
                    f"{name}_dimensions",
                    int(metadata["width"]) >= 600
                    and int(metadata["height"]) >= 600,
                    f"{metadata['width']}x{metadata['height']}",
                    ">=600x600",
                ),
            ]
        )

    with (media / "validation_v7_10_0.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["check", "passed", "actual", "expected"])
        writer.writerows(checks)
    failures = [name for name, passed, _, _ in checks if not passed]
    if failures:
        (media / "validation_status.txt").write_text(
            "failed\n", encoding="utf-8"
        )
        raise SystemExit(
            "V7.10.0 validation failed: " + ", ".join(failures)
        )
    (media / "validation_status.txt").write_text(
        "passed\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
