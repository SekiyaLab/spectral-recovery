from __future__ import annotations
import argparse,json
from pathlib import Path
from .experiment import mp_diagnostic,run,summarize
from .plotting import heatmap,mechanism,mp_histogram,recovery

def main()->None:
    parser=argparse.ArgumentParser(); parser.add_argument('--config',type=Path,default=Path('configs/primary.json')); args=parser.parse_args(); raw=Path('data/raw/primary-results.json'); rows=run(args.config,raw); summary=summarize(rows); Path('data/summary.json').write_text(json.dumps(summary,indent=2)+'\n'); config=json.loads(args.config.read_text()); values=mp_diagnostic(config,Path('data/raw/mp-null-eigenvalues.json')); figures=Path('figures'); figures.mkdir(exist_ok=True); mechanism(figures/'mechanism.png'); recovery(summary,figures/'recovery-curves.png'); heatmap(summary,figures/'improvement-heatmap.png'); mp_histogram(values,figures/'mp-null-sanity.png')
    for row in summary: print(json.dumps(row,sort_keys=True))

if __name__=='__main__': main()

