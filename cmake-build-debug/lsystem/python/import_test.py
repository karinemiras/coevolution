#!/usr/bin/env python3

import lsystem_python
import unittest


class ImportTest(unittest.TestCase):
    def test_import(self):
        self.assertIsNotNone(revolve_brain_python.Evaluator)
        self.assertIsNotNone(revolve_brain_python.Brain)
        self.assertIsNotNone(revolve_brain_python.NeuralNetwork)
        self.assertIsNotNone(revolve_brain_python.RLPower)
        self.assertIsNotNone(revolve_brain_python.CPGBrain)


if __name__ == "__main__":
    unittest.main()
