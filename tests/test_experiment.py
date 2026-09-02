import json
from pathlib import Path
from spectral_recovery.experiment import run,summarize

def test_seeded_run_is_deterministic_and_complete(tmp_path:Path):
    config={"seed":2,"dimensions":8,"replications":2,"q_values":[.5],"factor_correlations":[0,.2],"mp_diagnostic_replications":2}; path=tmp_path/'config.json'; path.write_text(json.dumps(config)); first=run(path,tmp_path/'a.json'); second=run(path,tmp_path/'b.json'); assert first==second; assert len(first)==12; summary=summarize(first); assert len(summary)==6; assert all(row['replications']==2 for row in summary)

