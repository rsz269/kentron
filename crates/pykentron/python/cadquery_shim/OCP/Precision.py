"""OCP.Precision: OCCT's tolerances."""


class Precision:
    @staticmethod
    def Confusion_s() -> float:
        return 1e-7

    @staticmethod
    def Angular_s() -> float:
        return 1e-12

    @staticmethod
    def Infinite_s() -> float:
        return 2e100

    @staticmethod
    def IsInfinite_s(r: float) -> bool:
        return abs(r) >= 0.5 * 2e100
