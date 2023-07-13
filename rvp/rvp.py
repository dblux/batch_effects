import numpy as np
import pandas as pd

from dataclasses import dataclass
from typing import Optional, Union


@dataclass
class RVP:
    percentage: float
    sum_squares: Optional[pd.DataFrame]

def rvp(
    X: Union[np.ndarray, pd.DataFrame],
    batch: pd.Series,
    cls: Optional[pd.Series]=None,
    return_obj: bool=False
):
    '''Calculates percentage of variance in data due to batch effects

    Args:
        X: {array-like, sparse matrix} of shape (n_samples, n_features)
        batch: array-like of shape (n_samples, )
        cls: array-like of shape (n_samples, )
        return_obj: bool indicating whether to return RVP dataclass
    
    Returns:
        Either a float indicating total percentage of variance in data due to
        batch effects or a dataclass RVP.
    '''

    if X.shape[0] != batch.shape[0]:
        raise Exception("Length of batch does not match number of rows in X!")
    
    if len(np.unique(batch)) == 1:
        return RVP(percentage=0, sum_squares=None)

    X = np.nan_to_num(X)

    # COMPUTE RVP
    if cls is None:
        feature_means = np.mean(X, axis=0)
        ss_total = np.sum((X - feature_means) ** 2, axis=0)
        X_batches = pd.DataFrame(X).groupby(batch)
        squares = (X_batches.mean() - feature_means) ** 2
        nperbatches = X_batches.size().to_numpy()[:, np.newaxis]
        ss_batch = np.sum(squares * nperbatches, axis=0)
        assert ss_batch.shape[0] == X.shape[1]
        pct_batch = np.sum(ss_batch) / np.sum(ss_total)
        
        if return_obj:
            sum_squares = pd.DataFrame({
                "ss_batch": ss_batch,
                "ss_total": ss_total
            })
            return RVP(percentage=pct_batch, sum_squares=sum_squares)
        else:
            return pct_batch
    else:
        ss_total = np.sum((X - X.mean(axis=0)) ** 2, axis=0)
        idx = pd.MultiIndex.from_frame(
            pd.DataFrame({"batch": batch, "class": cls})
        )
        if isinstance(X, pd.DataFrame):
            X.index = idx
            X_classes = X.groupby(level="class")
        else:
            X_classes = pd.DataFrame(X, index=idx).groupby(level="class")

        assert all(X_classes.size() != 0)
        
        objs = []
        for _, X_class in X_classes:
            objs.append(rvp(
                X_class,
                X_class.index.get_level_values("batch"),
                return_obj=True
            ))
        ss_batch_classes = [
            obj.sum_squares["ss_batch"]
            for obj in objs if obj.sum_squares is not None
        ]
        ss_batch = pd.DataFrame(ss_batch_classes).sum()
        assert ss_batch.shape[0] == ss_total.shape[0]
        pct_batch = np.sum(ss_batch) / np.sum(ss_total)

        if return_obj:
            sum_squares = pd.DataFrame({
                "ss_batch": ss_batch,
                "ss_total": ss_total
            })
            return RVP(percentage=pct_batch, sum_squares=sum_squares)
        else:
            return pct_batch
