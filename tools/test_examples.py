from pathlib import Path
import os
import subprocess
import sys

EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"

SUCCESS_EXAMPLES = {
    EXAMPLES_DIR / "bam.py",
    EXAMPLES_DIR / "daebuilder.py",
    EXAMPLES_DIR / "docs.py",
    EXAMPLES_DIR / "doube-diff.py",
    EXAMPLES_DIR / "if_else.py",
    EXAMPLES_DIR / "issues" / "dae-optimization.py",
    EXAMPLES_DIR / "minimal_example" / "CSTR.py",
    EXAMPLES_DIR / "minimal_example" / "engine.py",
    EXAMPLES_DIR / "minimal_example" / "test.py",
    EXAMPLES_DIR / "nle" / "bates_examples.py",
    EXAMPLES_DIR / "nle" / "bates_inference.py",
    EXAMPLES_DIR / "nle" / "cstr_collocation_sim.py",
    EXAMPLES_DIR / "nle" / "simple_mixer.py",
    EXAMPLES_DIR / "oed" / "quaglio.py",
    EXAMPLES_DIR / "oed" / "tools.py",
    EXAMPLES_DIR / "pe_gauss_newton.py",
    EXAMPLES_DIR / "performance" / "sparcity.py",
    EXAMPLES_DIR / "performance" / "testing_concepts.py",
    EXAMPLES_DIR / "pid.py",
    EXAMPLES_DIR / "sim_plot.py",
    EXAMPLES_DIR / "ternary_flash.py",
    EXAMPLES_DIR / "test_scaling.py",
    EXAMPLES_DIR / "time_derivative_dae.py",
}

IGNORED_EXAMPLES = {
    EXAMPLES_DIR / "nle" / "pe_direct.py",
    EXAMPLES_DIR / "dae_ode" / "hyfo_dae.py",
    EXAMPLES_DIR / "issues" / "zero_obj.py",
    EXAMPLES_DIR / "parameter_identifiabilty.py",
    EXAMPLES_DIR / "covariance_ellipse.py",
    EXAMPLES_DIR / "dario.py",
    EXAMPLES_DIR / "discrepancy.py",
    EXAMPLES_DIR / "inference_dynamic.py",
    EXAMPLES_DIR / "margules.py",
    EXAMPLES_DIR / "nle" / "direct_optimization.py",
    EXAMPLES_DIR / "nle" / "prediction_error.py",
    EXAMPLES_DIR / "nle" / "puromycin.py",
    EXAMPLES_DIR / "nle" / "scaling.py",
    EXAMPLES_DIR / "nle" / "variance_influence.py",
    EXAMPLES_DIR / "scipyelipse.py",
    EXAMPLES_DIR / "tt.py",
    EXAMPLES_DIR / "black_box_dynamic.py",
    EXAMPLES_DIR / "black_box_raoult.py",
    EXAMPLES_DIR / "cstr.py",
    EXAMPLES_DIR / "dae_ode" / "cstr.py",
    EXAMPLES_DIR / "dae_ode" / "free_fall.py",
    EXAMPLES_DIR / "dae_ode" / "ra_dae.py",
    EXAMPLES_DIR / "felix.py",
    EXAMPLES_DIR / "issues" / "covariance_manipulation.py",
    EXAMPLES_DIR / "issues" / "oed_jacobian.py",
    EXAMPLES_DIR / "mpc.py",
    EXAMPLES_DIR / "multistart.py",
    EXAMPLES_DIR / "nle" / "berty.py",
    EXAMPLES_DIR / "nle" / "monte_carlo_covariance.py",
    EXAMPLES_DIR / "nle" / "oed_nle.py",
    EXAMPLES_DIR / "nle" / "vle.py",
    EXAMPLES_DIR / "nle" / "vle_pe.py",
    EXAMPLES_DIR / "oed" / "article.py",
    EXAMPLES_DIR / "oed" / "deluca2016.py",
    EXAMPLES_DIR / "oed" / "hoang.py",
    EXAMPLES_DIR / "optimize.py",
    EXAMPLES_DIR / "playground.py",
    EXAMPLES_DIR / "saskia_model_analysis.py",
    EXAMPLES_DIR / "sci.py",
    EXAMPLES_DIR / "sensitivity_analysis.py",
    EXAMPLES_DIR / "store_results.py",
    EXAMPLES_DIR / "tikhonov_tuning.py",
    EXAMPLES_DIR / "unknown_uncertainty_pe.py",
}

IGNORED_EXAMPLES.update(SUCCESS_EXAMPLES)

ENV = os.environ.copy()
ENV["MPLBACKEND"] = "Agg"

TIMEOUT_SECONDS = 15


def format_example_path(example: Path) -> str:
    relative = example.relative_to(EXAMPLES_DIR)
    parts = " / ".join(f'"{part}"' for part in relative.parts)
    return f"    EXAMPLES_DIR / {parts},"


def test_examples_run():
    failed = []
    timed_out = []
    worked = []

    for example in sorted(EXAMPLES_DIR.rglob("*.py")):
        if example in IGNORED_EXAMPLES:
            continue

        print(f"RUNNING: {example}")

        try:
            result = subprocess.run(
                [sys.executable, str(example)],
                cwd=EXAMPLES_DIR.parent,
                capture_output=True,
                text=True,
                env=ENV,
                timeout=TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired:
            timed_out.append(example)
            print(f"TIMEOUT: {example}")
            continue

        if result.returncode == 0:
            worked.append(example)
            print(f"OK: {example}")
        else:
            failed.append(example)
            print(f"FAILED: {example}")

    print("\n\nSUCCESSFUL RUNS - copy into IGNORED_EXAMPLES if you want to skip them:")
    print("IGNORED_EXAMPLES = {")
    for example in worked:
        print(format_example_path(example))
    print("}")

    if timed_out:
        print("\n\nTIMED OUT:")
        for example in timed_out:
            print(format_example_path(example))

    print("Failes runs = {")
    for example in failed:
        print(format_example_path(example))
    print("}")


if __name__ == "__main__":
    test_examples_run()
