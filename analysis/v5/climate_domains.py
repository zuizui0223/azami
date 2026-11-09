"""Spatial holdout comparison of climate domains on the original V5 S1 cells."""
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

DOMAINS = {
    "temperature": ["bio01", "bio04", "bio05", "bio06"],
    "precipitation": ["bio12", "bio13", "bio14", "bio15", "bio17", "bio18"],
    "radiation": ["rsds_mean"],
}


def predict(x, y, weights, train, test, penalty):
    w = weights[train] / weights[train].sum()
    mean = np.average(x[train], weights=w, axis=0)
    sd = np.sqrt(np.average((x[train]-mean)**2, weights=w, axis=0))
    sd[sd < 1e-10] = 1
    fit = Ridge(alpha=penalty).fit((x[train]-mean)/sd, y[train], sample_weight=w)
    return fit.predict((x[test]-mean)/sd)


def evaluate(cells, label, blocks, weights):
    domains = list(DOMAINS)
    models = {"annual": ["bio01","bio12","rsds_mean"], "warm": ["bio05","bio18","rsds_mean"]}
    for omit in [None, *domains]:
        models["all" if omit is None else "without_"+omit] = [
            axis for domain in domains if domain != omit for axis in DOMAINS[domain]]
    y = cells[label].to_numpy()
    groups = np.asarray(blocks)
    rows = []
    for model, axes in models.items():
        means = cells.groupby("accepted_key")[axes].transform("mean")
        x = np.column_stack([means, cells[axes]-means])
        for fold,(tr,te) in enumerate(GroupKFold(5).split(x,groups=groups)):
            penalties = [0.0001,0.001,0.01,0.1,1,10]
            scores = []
            inner = list(GroupKFold(4).split(tr, groups=groups[tr]))
            for penalty in penalties:
                score = 0
                for a,b in inner:
                    err = y[tr[b]] - predict(x,y,weights,tr[a],tr[b],penalty)
                    score += np.sum(weights[tr[b]]*err**2)
                scores.append(score)
            penalty = penalties[np.argmin(scores)]
            pred = predict(x,y,weights,tr,te,penalty)
            baseline = np.average(y[tr],weights=weights[tr])
            rows.append(dict(trait=label,model=model,fold=fold,n_axes=len(axes),penalty=penalty,
                             error=np.sum(weights[te]*(y[te]-pred)**2),
                             null_error=np.sum(weights[te]*(y[te]-baseline)**2)))
    return rows
