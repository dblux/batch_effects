import unittest
import pandas as pd
import rvp


class TestComplexData(unittest.TestCase):
    def setUp(self):
        # Data with missing values and an entire class from the same batch
        file = "data/rvp_data.csv"
        self.X = pd.read_csv(file, index_col=[0, 1])
        self.batch = self.X.index.get_level_values("batch")
        self.cls = self.X.index.get_level_values("class")

    def test_rvp(self):
        pct_batch = rvp.rvp(self.X, self.batch)
        pct_batch = rvp.rvp(self.X, self.batch, self.cls)
        RVP = rvp.rvp(self.X, self.batch, return_obj=True)
        RVP = rvp.rvp(self.X, self.batch, self.cls, return_obj=True)
