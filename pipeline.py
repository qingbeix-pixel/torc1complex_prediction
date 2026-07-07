#!/usr/bin/env python3
"""
======================================================================
 Plant TORC1 Structure Prediction Pipeline (Boltz-2)
======================================================================

De novo AI-based structure prediction of the Arabidopsis thaliana
TORC1 complex (TOR + RAPTOR1B + LST8-1) using Boltz-2.

当前部署: Boltz-2
未来支持: AlphaFold3（配置模板见 af3/config/config.yaml）

Protocol summary:
  Stage 1 — Fetch & validate sequences (TAIR/UniProt)
  Stage 2 — Build multi-chain input (YAML/FASTA, no templates)
  Stage 3 — Run Boltz-2 prediction
  Stage 4 — Analyse results (pLDDT, PAE, interfaces, H-bonds)
  Stage 5 — Downstream analysis (ConSurf, FoldX, MD)

Usage:
  python pipeline.py                  # Full pipeline
  python pipeline.py --stage 1        # Only input preparation
  python pipeline.py --stage 2        # Only input building
  python pipeline.py --stage 3        # Only prediction
  python pipeline.py --stage 4        # Only analysis
  python pipeline.py --stage 5        # Only downstream
  python pipeline.py --config my.yaml # Custom config

Author: torc1predict pipeline
======================================================================
"""

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

import yaml

from modules.input_prep import prepare_sequences, SequenceReport
from modules.input_builder import build_inputs
from modules.prediction import run_prediction, parse_rankings
from modules.analysis import analyze_structure, parse_plddt_from_b_factor
from modules.downstream import run_downstream_analyses


# ---------------------------------------------------------------------------
# Pipeline runner
# ---------------------------------------------------------------------------

class Torc1Pipeline:
    """Orchestrates the complete TORC1 structure prediction workflow."""

    def __init__(self, config_path: str = "boltz/config/config.yaml"):
        with open(config_path) as f:
            self.config = yaml.safe_load(f)

        self._setup_logging()
        self.results = {}  # Accumulated results from all stages
        self.timings = {}  # Stage timing

    def _setup_logging(self):
        """Configure logging to both console and file."""
        log_dir = self.config["output"].get("subdirs", {}).get("logs", "logs")
        os.makedirs(log_dir, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(log_dir, f"pipeline_{timestamp}.log")

        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            handlers=[
                logging.StreamHandler(sys.stdout),
                logging.FileHandler(log_file),
            ],
        )
        self.logger = logging.getLogger("torc1predict")
        self.logger.info(f"Log file: {log_file}")

    # ------------------------------------------------------------------
    # Stage 1: Sequence preparation
    # ------------------------------------------------------------------

    def stage1_prepare_sequences(self) -> Dict[str, SequenceReport]:
        """Fetch and validate all target sequences."""
        self.logger.info("\n" + "█"*70)
        self.logger.info(" STAGE 1: Sequence Input Preparation")
        self.logger.info("█"*70)

        t0 = time.time()

        reports = prepare_sequences(
            targets=self.config["targets"],
            source=self.config["input"]["source"],
            checks=self.config.get("sequence_checks"),
            fasta_dir=self.config["input"]["fasta_dir"],
        )

        # Print summary
        self.logger.info("\n--- Stage 1 Summary ---")
        for key, rep in reports.items():
            status = "✓" if rep.valid else "✗"
            self.logger.info(
                f"  {status} {rep.name} ({rep.gene_id}): "
                f"{rep.length} aa, {len(rep.domains)} domains"
            )
            for w in rep.warnings:
                self.logger.warning(f"    ⚠ {w}")

        elapsed = time.time() - t0
        self.timings["stage1"] = elapsed
        self.results["sequences"] = reports

        return reports

    # ------------------------------------------------------------------
    # Stage 2: Input building
    # ------------------------------------------------------------------

    def stage2_build_inputs(self) -> Dict:
        """Generate multi-chain input files for the prediction engine."""
        self.logger.info("\n" + "█"*70)
        self.logger.info(" STAGE 2: Multi-chain Input Building")
        self.logger.info("█"*70)

        t0 = time.time()

        reports = self.results.get("sequences")
        if not reports:
            raise RuntimeError("Stage 1 must be run before Stage 2")

        # 当前只部署了 Boltz-2
        engine = "boltz2"
        self.logger.info(f"Building input for: {engine}")

        inputs = build_inputs(
            reports=reports,
            stoichiometry=self.config["complex"]["chains"],
            engine=engine,
            output_dir=self.config["input"]["fasta_dir"],
            num_samples=self.config["prediction"].get("diffusion_samples", 1),
        )

        elapsed = time.time() - t0
        self.timings["stage2"] = elapsed
        self.results["inputs"] = inputs

        self.logger.info(f"\n--- Stage 2 Summary ---")
        self.logger.info(f"  Engine: {inputs['engine']}")
        self.logger.info(f"  Input:  {inputs['input_file']}")
        self.logger.info(f"  FASTA:  {inputs['fasta_file']}")

        return inputs

    # ------------------------------------------------------------------
    # Stage 3: Structure prediction
    # ------------------------------------------------------------------

    def stage3_run_prediction(self) -> Dict:
        """Execute the AI structure prediction (Boltz-2)."""
        self.logger.info("\n" + "█"*70)
        self.logger.info(" STAGE 3: De Novo Structure Prediction")
        self.logger.info("█"*70)

        t0 = time.time()

        inputs = self.results.get("inputs")
        if not inputs:
            raise RuntimeError("Stage 2 must be run before Stage 3")

        # 当前只用 Boltz-2
        engine = "boltz2"
        boltz_config = self.config.get("prediction", {})

        # 检查是否配置了 AF3（未部署）
        if self.config.get("alphafold3"):
            self.logger.warning(
                "检测到 alphafold3 配置，但 AF3 尚未部署。将使用 Boltz-2 运行。"
            )

        self.logger.info(f"Running {engine.upper()} prediction...")
        self.logger.info(f"This may take hours for a ~4,000 residue complex.")

        result = run_prediction(
            input_file=inputs["input_file"],
            output_dir=self.config["output"]["base_dir"],
            engine=engine,
            engine_config=boltz_config,
            msa_databases=boltz_config.get("msa_databases", {}),
        )

        # Parse rankings
        rankings = parse_rankings(result.get("rankings", ""))
        if rankings:
            self.logger.info("\nModel rankings:")
            for r in rankings:
                self.logger.info(
                    f"  #{r['rank']} {r.get('seed', r.get('model', '?'))}: "
                    f"ipTM={r.get('iptm', 'N/A')}, pTM={r.get('ptm', 'N/A')}"
                )

        elapsed = time.time() - t0
        self.timings["stage3"] = elapsed
        self.results["prediction"] = result
        self.results["rankings"] = rankings

        # Summarize affinity results
        affinity = result.get("affinity", {})
        if affinity:
            if "affinity_probability_binary" in affinity:
                self.logger.info(
                    f"  Binder prob: {affinity['affinity_probability_binary']:.3f}"
                )
            if "affinity_pred_value" in affinity:
                self.logger.info(
                    f"  Affinity:    log10(IC50) = {affinity['affinity_pred_value']:.3f} (in μM)"
                )

        self.logger.info(f"\n--- Stage 3 Summary ---")
        self.logger.info(f"  Engine:     {result['engine']}")
        self.logger.info(f"  Wall time:  {elapsed/3600:.1f} hours")
        self.logger.info(f"  Models:     {len(result['models'])} structures")

        return result

    # ------------------------------------------------------------------
    # Stage 4: Post-prediction analysis
    # ------------------------------------------------------------------

    def stage4_analyze(self) -> Dict:
        """Analyze predicted structures."""
        self.logger.info("\n" + "█"*70)
        self.logger.info(" STAGE 4: Post-Prediction Analysis")
        self.logger.info("█"*70)

        t0 = time.time()

        prediction = self.results.get("prediction")
        if not prediction:
            raise RuntimeError("Stage 3 must be run before Stage 4")

        analysis_config = self.config.get("analysis", {})
        out_config = self.config["output"]
        plots_dir = os.path.join(out_config["base_dir"], out_config["subdirs"]["plots"])

        # Analyze the top-ranked model
        model_files = prediction["models"]
        if not model_files:
            self.logger.error("No model files found!")
            return {}

        # Sort by expected ranking (confidence_score)
        rankings = self.results.get("rankings", [])
        if rankings:
            top_ranked = rankings[0]
            top_model = model_files[0]  # Default first
            # rankings from _parse_boltz_confidence have 'model_file' key
            if top_ranked.get("model_file"):
                top_model = top_ranked["model_file"]
            else:
                for mf in model_files:
                    if top_ranked.get("model", "") in mf:
                        top_model = mf
                        break
        else:
            top_model = model_files[0]

        self.logger.info(f"Analyzing top-ranked model: {top_model}")

        # Find official Boltz auxiliary files for the selected model.
        pred_dir = os.path.join(out_config["base_dir"], "predictions")
        model_stem = Path(top_model).stem

        def _pick_matching_file(files, stem):
            for path in files:
                if stem in Path(path).stem:
                    return str(path)
            return str(files[0]) if files else None

        pae_files = (
            sorted(Path(pred_dir).rglob("pae_*.npz")) +
            sorted(Path(out_config["base_dir"], out_config["subdirs"]["pae"]).glob("pae_*.npz")) +
            sorted(Path(out_config["base_dir"], out_config["subdirs"]["pae"]).glob("pae_*.npy")) +
            sorted(Path(out_config["base_dir"], out_config["subdirs"]["pae"]).glob("pae_*.json"))
        )
        plddt_files = (
            sorted(Path(pred_dir).rglob("plddt_*.npz")) +
            sorted(Path(out_config["base_dir"]).rglob("plddt_*.npz"))
        )
        pae_path = _pick_matching_file(pae_files, model_stem)
        plddt_path = _pick_matching_file(plddt_files, model_stem)

        if pae_path:
            self.logger.info(f"Using PAE file: {pae_path}")
        if plddt_path:
            self.logger.info(f"Using pLDDT file: {plddt_path}")

        # Infer chain ranges
        reports = self.results.get("sequences", {})
        chain_ranges = {}
        offset = 0
        for chain_label, target_key in self.config["complex"]["chains"].items():
            rep = reports.get(target_key)
            if rep and rep.length:
                chain_ranges[chain_label] = (offset, offset + rep.length)
                offset += rep.length

        analysis_results = analyze_structure(
            model_path=top_model,
            pae_path=pae_path,
            plddt_path=plddt_path,
            chain_ranges=chain_ranges,
            config=analysis_config,
            output_dir=plots_dir,
        )

        elapsed = time.time() - t0
        self.timings["stage4"] = elapsed
        self.results["analysis"] = analysis_results

        # Print summary
        self._print_analysis_summary(analysis_results)

        return analysis_results

    def _print_analysis_summary(self, results: Dict):
        """Print a concise analysis summary."""
        self.logger.info(f"\n--- Stage 4 Summary ---")

        plddt = results.get("plddt", {}).get("per_chain", {})
        for ch, stats in plddt.items():
            self.logger.info(
                f"  Chain {ch}: pLDDT={stats['mean']:.1f}±{stats['std']:.1f}, "
                f"{stats['fraction_confident']*100:.0f}% confident (≥70)"
            )

        pae = results.get("pae", {})
        for pair, s in pae.get("inter_chain", {}).items():
            self.logger.info(
                f"  Interface {pair}: PAE={s['mean']:.2f} Å, "
                f"{s['fraction_below_threshold']*100:.0f}% ≤ 5 Å"
            )

        ifaces = results.get("interfaces", {})
        for pair, data in ifaces.items():
            self.logger.info(
                f"  Interface {pair}: {data['n_residue_pairs']} contact residue pairs"
            )

        self.logger.info(
            f"  H-bonds: {results.get('hydrogen_bonds', {}).get('count', 0)} inter-chain"
        )
        self.logger.info(
            f"  Hydrophobic: {results.get('hydrophobic_contacts', {}).get('count', 0)} residue pairs"
        )

        # Affinity interpretation
        affinity = self.results.get("prediction", {}).get("affinity", {})
        if affinity:
            aff_config = self.config.get("analysis", {}).get("affinity", {})
            self._print_affinity_summary(affinity, aff_config)

    def _print_affinity_summary(self, affinity: Dict, aff_config: Dict):
        """Print binding affinity prediction summary (Boltz-2 only)."""
        self.logger.info(f"\n--- Affinity Prediction ---")

        prob = affinity.get("affinity_probability_binary")
        pred_val = affinity.get("affinity_pred_value")

        if prob is not None:
            threshold = aff_config.get("binder_threshold", 0.5)
            is_binder = prob >= threshold
            self.logger.info(
                f"  Binder probability: {prob:.3f} "
                f"({'✓ Predicted binder' if is_binder else '✗ Predicted non-binder'}, "
                f"threshold={threshold})"
            )

        if pred_val is not None:
            strong = aff_config.get("strong_binder_threshold", -2.0)
            weak = aff_config.get("weak_binder_threshold", 0.0)
            self.logger.info(
                f"  Affinity (log10 IC50): {pred_val:.3f} "
                f"(IC50 = {10**pred_val:.2f} μM)"
            )

    # ------------------------------------------------------------------
    # Stage 5: Downstream analysis
    # ------------------------------------------------------------------

    def stage5_downstream(self) -> Dict:
        """Run downstream analyses."""
        self.logger.info("\n" + "█"*70)
        self.logger.info(" STAGE 5: Downstream Analysis")
        self.logger.info("█"*70)

        t0 = time.time()

        analysis = self.results.get("analysis")
        prediction = self.results.get("prediction")
        if not analysis or not prediction:
            raise RuntimeError("Stage 3 & 4 must be run before Stage 5")

        # Get chain sequences
        reports = self.results.get("sequences", {})
        chain_sequences = {}
        for chain_label, target_key in self.config["complex"]["chains"].items():
            rep = reports.get(target_key)
            if rep:
                chain_sequences[chain_label] = rep.sequence

        # Find top model
        top_model = analysis.get("model", prediction["models"][0])

        # Gather interface residues
        interface_residues = {
            pair: data.get("residues", [])
            for pair, data in analysis.get("interfaces", {}).items()
        }

        downstream_dir = os.path.join(
            self.config["output"]["base_dir"], "downstream",
        )

        downstream_results = run_downstream_analyses(
            model_path=top_model,
            interface_residues=interface_residues,
            chain_sequences=chain_sequences,
            config=self.config.get("analysis", {}),
            output_dir=downstream_dir,
        )

        elapsed = time.time() - t0
        self.timings["stage5"] = elapsed
        self.results["downstream"] = downstream_results

        return downstream_results

    # ------------------------------------------------------------------
    # Run pipeline
    # ------------------------------------------------------------------

    def run(self, stages: Optional[list] = None):
        """
        Execute the full pipeline or selected stages.

        Parameters
        ----------
        stages : list of int or None
            Stages to run (1-5). None runs all.
        """
        if stages is None:
            stages = [1, 2, 3, 4, 5]

        self.logger.info("="*70)
        self.logger.info(" Plant TORC1 Structure Prediction Pipeline")
        self.logger.info(f" Complex: {self.config['complex']['name']}")
        self.logger.info(f" Engine:  Boltz-2")
        self.logger.info(f" Stages:  {stages}")
        self.logger.info("="*70)

        pipeline_start = time.time()

        stage_map = {
            1: ("Sequence Preparation", self.stage1_prepare_sequences),
            2: ("Input Building", self.stage2_build_inputs),
            3: ("Structure Prediction", self.stage3_run_prediction),
            4: ("Post-Prediction Analysis", self.stage4_analyze),
            5: ("Downstream Analysis", self.stage5_downstream),
        }

        for stage_num in stages:
            if stage_num not in stage_map:
                self.logger.error(f"Unknown stage: {stage_num}")
                continue

            name, func = stage_map[stage_num]
            self.logger.info(f"\n>>> Running Stage {stage_num}: {name}")

            try:
                func()
            except Exception as e:
                self.logger.error(f"Stage {stage_num} failed: {e}", exc_info=True)
                if stage_num >= 3:  # Don't continue after prediction failure
                    raise
                self.logger.warning("Continuing despite error...")

        # Final summary
        total_time = time.time() - pipeline_start
        self.logger.info("\n" + "="*70)
        self.logger.info(" PIPELINE COMPLETE")
        self.logger.info("="*70)
        self.logger.info(f" Total wall time: {total_time/3600:.1f} hours")
        for stage_num, elapsed in sorted(self.timings.items()):
            self.logger.info(f"  Stage {stage_num}: {elapsed/60:.1f} min")

        # Save consolidated summary
        self._save_summary(total_time)

    def _save_summary(self, total_time: float):
        """Write a summary JSON of the entire pipeline run."""
        summary = {
            "pipeline": "torc1predict",
            "complex": self.config["complex"]["name"],
            "engine": "boltz2",
            "timestamp": datetime.now().isoformat(),
            "total_wall_time_s": total_time,
            "stages": {},
        }

        for i in range(1, 6):
            if f"stage{i}" in self.timings:
                summary["stages"][f"stage{i}"] = {
                    "elapsed_s": self.timings[f"stage{i}"],
                }

        # Add key metrics
        rankings = self.results.get("rankings", [])
        if rankings:
            top = rankings[0]
            summary["top_model"] = {
                "rank": 1,
                "iptm": top.get("iptm"),
                "ptm": top.get("ptm"),
                "seed": top.get("seed"),
            }

        # Affinity summary (Boltz-2)
        affinity = self.results.get("prediction", {}).get("affinity", {})
        if affinity:
            summary["affinity"] = {
                k: v for k, v in affinity.items()
                if k in ("affinity_pred_value", "affinity_probability_binary")
            }

        analysis = self.results.get("analysis", {})
        plddt = analysis.get("plddt", {}).get("per_chain", {})
        summary["plddt"] = {
            ch: {"mean": s["mean"], "fraction_confident": s["fraction_confident"]}
            for ch, s in plddt.items()
        }

        summary_path = os.path.join(
            self.config["output"]["base_dir"], "pipeline_summary.json",
        )
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2, default=str)

        self.logger.info(f"\nSummary saved: {summary_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plant TORC1 Structure Prediction Pipeline (Boltz-2)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pipeline.py                        # Full pipeline
  python pipeline.py --stage 1               # Only fetch sequences
  python pipeline.py --stage 1 2             # Fetch + build inputs
  python pipeline.py --stage 4               # Only analyse existing results
  python pipeline.py --config custom.yaml    # Use custom config

Note: 当前仅部署了 Boltz-2 引擎。AlphaFold3 暂未部署，
如需使用请参考 af3/config/config.yaml 中 alphafold3 的配置说明。
        """,
    )
    parser.add_argument(
        "--stage", type=int, nargs="+", default=None,
        help="Stage(s) to run: 1=sequences, 2=inputs, 3=predict, 4=analyse, 5=downstream",
    )
    parser.add_argument(
        "--config", type=str, default="boltz/config/config.yaml",
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate config and inputs without running predictions",
    )

    args = parser.parse_args()

    # Validate config path
    if not os.path.exists(args.config):
        print(f"Error: config file not found: {args.config}")
        sys.exit(1)

    pipeline = Torc1Pipeline(args.config)

    if args.dry_run:
        pipeline.logger.info("DRY RUN — validating configuration...")
        pipeline.stage1_prepare_sequences()
        pipeline.stage2_build_inputs()
        pipeline.logger.info("Dry run complete. Configuration valid.")
        return

    pipeline.run(stages=args.stage)


if __name__ == "__main__":
    main()
