from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def style(): plt.rcParams.update({"font.family":"DejaVu Sans","axes.spines.top":False,"axes.spines.right":False,"figure.dpi":160,"axes.titleweight":"bold"})

def mechanism(path: Path) -> None:
    style(); fig,axes=plt.subplots(1,3,figsize=(12,3.8),layout="constrained"); n=32; rho=.3; truth=(1-rho)*np.eye(n)+rho*np.ones((n,n)); rng=np.random.default_rng(4); raw=np.corrcoef(rng.normal(size=(70,n)) @ np.linalg.cholesky(truth).T,rowvar=False)
    for ax,matrix,title in zip(axes,[truth,raw,np.where(np.abs(raw)>.18,raw,0)], ["Known population correlation","Finite-sample correlation","Spectral idea: retain structure"]):
        image=ax.imshow(matrix,cmap="coolwarm",vmin=-.5,vmax=1); ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(image,ax=axes,shrink=.75,label="correlation"); fig.suptitle("One-factor signal is observed through finite-sample noise",weight="bold"); fig.savefig(path,bbox_inches="tight"); plt.close(fig)

def recovery(summary:list[dict],path:Path)->None:
    style(); fig,axes=plt.subplots(1,3,figsize=(12,3.8),sharey=True); qs=sorted({r['q'] for r in summary}); colours={"raw":"#D94C3D","mp_clip":"#245EAA","constant_shrinkage":"#16806A"}
    for ax,q in zip(axes,qs):
        for method in colours:
            rows=sorted([r for r in summary if r['q']==q and r['method']==method],key=lambda x:x['rho']); x=np.array([r['rho'] for r in rows]); y=np.array([r['mean_error'] for r in rows]); err=np.array([(r['mean_error_ci_high']-r['mean_error_ci_low'])/2 for r in rows]); ax.errorbar(x,y,yerr=err,marker='o',label=method.replace('_',' '),color=colours[method],capsize=3)
        ax.set_title(f"q = {q}"); ax.set_xlabel("factor correlation ρ")
    axes[0].set_ylabel("off-diagonal recovery error"); axes[0].legend(frameon=False,fontsize=8); fig.suptitle("Recovery error by signal strength and matrix aspect ratio",weight="bold"); fig.tight_layout(); fig.savefig(path,bbox_inches="tight"); plt.close(fig)

def heatmap(summary:list[dict],path:Path)->None:
    style(); fig,axes=plt.subplots(1,2,figsize=(10,4),layout="constrained"); methods=["mp_clip","constant_shrinkage"]; qs=sorted({r['q'] for r in summary}); rhos=sorted({r['rho'] for r in summary})
    for ax,method in zip(axes,methods):
        data=np.array([[next(r['paired_difference_vs_raw'] for r in summary if r['q']==q and r['rho']==rho and r['method']==method) for rho in rhos] for q in qs]); im=ax.imshow(data,cmap="RdBu_r",vmin=-max(abs(data.min()),abs(data.max())),vmax=max(abs(data.min()),abs(data.max())))
        ax.set(title=method.replace('_',' ') + " minus raw",xticks=range(len(rhos)),xticklabels=rhos,yticks=range(len(qs)),yticklabels=qs,xlabel="factor correlation ρ",ylabel="q = N/T")
        for i in range(len(qs)):
            for j in range(len(rhos)): ax.text(j,i,f"{data[i,j]:+.3f}",ha="center",va="center",fontsize=8)
    fig.colorbar(im,ax=axes,label="paired mean error difference (negative favours method)"); fig.suptitle("Decision map: improvement relative to raw sample correlation",weight="bold"); fig.savefig(path,bbox_inches="tight"); plt.close(fig)

def mp_histogram(values:np.ndarray,path:Path)->None:
    style(); q=.5; lower=(1-np.sqrt(q))**2; upper=(1+np.sqrt(q))**2; x=np.linspace(lower+.001,upper-.001,400); density=np.sqrt((upper-x)*(x-lower))/(2*np.pi*q*x)
    fig,ax=plt.subplots(figsize=(8,4.5)); ax.hist(values,bins=65,density=True,alpha=.6,color="#245EAA",label="pooled null sample-covariance eigenvalues"); ax.plot(x,density,color="#D94C3D",lw=2,label="MP density, q=0.5"); ax.axvline(upper,color="#172033",ls='--',label=f"upper edge = {upper:.2f}"); ax.set(xlabel="eigenvalue",ylabel="density",title="MP null sanity diagnostic (finite sample)"); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(path,bbox_inches="tight"); plt.close(fig)
