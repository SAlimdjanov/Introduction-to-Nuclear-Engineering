"""
answer_format.py

Template objects for solution formats, as a way to clean up problem type annotations

"""

from typing import List, Tuple, Union


class SingleAnswer:
    """Object to store a single numeric answer with its associated units and significant figures"""

    def __init__(self, ans: float, units: Union[str, None], sig_figs: int) -> None:
        self.ans = ans
        self.units = units
        self.sig_figs = sig_figs

    def format(self) -> Tuple[float, str, int]:
        """Returns a tuple of the form (ans, units, sig_figs)"""
        if self.units is None:
            return self.ans, "", self.sig_figs
        return self.ans, self.units, self.sig_figs


class MultiPartAnswer:
    """Object to store multiple numeric answers The units and significant figures args can be lists
    in as required"""

    def __init__(
        self,
        ans: List[float],
        units: Union[str, List[str]],
        sig_figs: Union[int, List[int]],
    ) -> None:
        self.ans = ans
        self.units = units
        self.sig_figs = sig_figs

    def format(self) -> List[Tuple[float, str, int]]:
        """Returns a list of tuples of the form (ans, units, sig_figs)"""
        n = len(self.ans)
        if not isinstance(self.units, list):
            units = [self.units] * n
        else:
            units = self.units
        if not isinstance(self.sig_figs, list):
            sig_figs = [self.sig_figs] * n
        else:
            sig_figs = self.sig_figs
        return list((a, u, s) for a, u, s in zip(self.ans, units, sig_figs))
