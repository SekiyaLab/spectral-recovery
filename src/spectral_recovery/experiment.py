from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from .matrices import constant_shrinkage, mp_clip, offdiagonal_error, population_correlation, sample_correlation


def run(config_path: Path, output_path: Path) -> list[dict]:
    config = json.loads(config_path.read_text()); rng = np.random.default_rng(config["seed"]); n = config["dimensions"]; rows=[]
    for q in config["q_values"]:
        observations = round(n / q)
        for rho in config["factor_correlations"]:
            truth=population_correlation(n,rho); chol=np.linalg.cholesky(truth)
            for replication in range(config["replications"]):
                samples=rng.normal(size=(observations,n)) @ chol.T; raw=sample_correlation(samples)
                estimates={"raw":raw,"mp_clip":mp_clip(raw,q),"constant_shrinkage":constant_shrinkage(raw,q)}
                for method,estimate in estimates.items():
                    rows.append({"q":q,"observations":observations,"rho":rho,"replication":replication,"method":method,"error":offdiagonal_error(estimate,truth),"top_eigenvalue":float(np.linalg.eigvalsh(estimate)[-1])})
    output_path.parent.mkdir(parents=True,exist_ok=True); output_path.write_text(json.dumps({"config":config,"rows":rows},indent=2)+"\n"); return rows


def summarize(rows: list[dict]) -> list[dict]:
    summary=[]
    for q,rho,method in sorted({(r["q"],r["rho"],r["method"]) for r in rows}):
        values=np.array([r["error"] for r in rows if (r["q"],r["rho"],r["method"])==(q,rho,method)])
        raw=np.array([r["error"] for r in rows if (r["q"],r["rho"],r["method"])==(q,rho,"raw")])
        paired=values-raw
        summary.append({"q":q,"rho":rho,"method":method,"replications":len(values),"mean_error":float(values.mean()),"sd_error":float(values.std(ddof=1)),"mean_error_ci_low":float(values.mean()-1.96*values.std(ddof=1)/np.sqrt(len(values))),"mean_error_ci_high":float(values.mean()+1.96*values.std(ddof=1)/np.sqrt(len(values))),"paired_difference_vs_raw":float(paired.mean())})
    return summary


def mp_diagnostic(config: dict, output_path: Path) -> np.ndarray:
    rng=np.random.default_rng(config["seed"]+1); n=config["dimensions"]; q=.5; t=round(n/q); values=[]
    for _ in range(config["mp_diagnostic_replications"]):
        x=rng.normal(size=(t,n)); covariance=x.T@x/t; values.extend(np.linalg.eigvalsh(covariance))
    values=np.asarray(values); output_path.write_text(json.dumps({"q":q,"dimensions":n,"observations":t,"eigenvalues":values.tolist()},indent=2)+"\n"); return values

