import unittest
from your_project import your_project
import os
import sys
from functools import partial
if sys.version[0] == '2':
        import __builtin__ as builtins  # pylint:disable=import-error
else:
        import builtins  # pylint:disable=import-error

THISD = os.path.dirname(os.path.abspath(__file__))
here = partial(os.path.join, THISD)
class Test(unittest.TestCase):
    def test_returns_true(self):
        self.assertEquals(True, your_project.returns_true())
        self.assertTrue(your_project.returns_true())

    def test_will_pass(self):
        self.assertTrue(True)
