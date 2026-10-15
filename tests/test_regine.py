"""Régine 2026, Python version — tests. Run from the repository root:  python3 tests/test_regine.py"""
import itertools
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import regine  # noqa: E402

DATA = regine.load(os.path.join(os.path.dirname(__file__), "..", "data"))
NITRATE = regine.read_genes(os.path.join(os.path.dirname(__file__), "..", "data", "example_nitrate_list.txt"))


class TestRegine(unittest.TestCase):
    def test_data(self):
        self.assertEqual(len(NITRATE), 13)
        self.assertEqual(sum(len(s) for s in DATA["dap"].targets), 2848928)
        self.assertEqual(len(DATA["target"].cols), 65)
        self.assertEqual(len({c["agi"] for c in DATA["target"].cols}), 57)
        self.assertEqual(sum(len(s) for s in DATA["target_up"].targets) + sum(len(s) for s in DATA["target_down"].targets), 139241)
        self.assertEqual(len(DATA["target"].universe), 18727)

    def test_statistics(self):
        # R: phyper(9, 829, 32606 - 829, 13, lower.tail = FALSE)
        self.assertAlmostEqual(regine.hyper_upper(10, 32606, 829, 13) / 2.856e-14, 1, delta=0.01)
        self.assertEqual([round(x, 6) for x in regine.bh([0.01, 0.04, 0.03])], [0.03, 0.04, 0.04])
        self.assertEqual(regine.auc([1, 1], [1, 1]), 0.5)
        self.assertEqual(regine.auc([2, 3], [0, 1]), 1)
        self.assertAlmostEqual(regine.fisher2(0.05, 0.05), 0.0175, delta=5e-4)

    def test_real_is_exact(self):
        net, k = DATA["dap"], 4
        B = regine.background(DATA, net, "network") | set(NITRATE)
        rows = regine.enrich(net, NITRATE, B, k)
        kernels = list(itertools.combinations(NITRATE, k))
        self.assertEqual(len(kernels), 715)
        for r in rows[:60]:
            s = net.targets[r["col"]]
            self.assertAlmostEqual(r["real"], sum(all(g in s for g in kk) for kk in kernels) / len(kernels), places=12)

    def test_nitrate(self):
        res = regine.run(DATA, NITRATE, reps=500)
        dap = min(res["dap"], key=lambda r: r["p"])
        self.assertEqual((dap["symbol"], dap["kt"], dap["Kt"]), ("NLP7", 10, 829))
        self.assertLess(dap["fdr"], 1e-10)
        tgt = min(res["target"], key=lambda r: r["p"])
        self.assertEqual((tgt["symbol"], tgt["effect"], tgt["k_up"], tgt["k_down"]), ("NLP7", "activates", 5, 0))
        self.assertEqual(res["both"][0]["symbol"], "NLP7")
        self.assertEqual(res["both"][0]["bound_and_regulated"], 4)
        self.assertTrue(all(p["auc"] > 0.5 for p in res["profiles"]["dap"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
