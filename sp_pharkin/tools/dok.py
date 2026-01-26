"""Quick demo script for sp-pharkin calculations."""

from sp_pharkin import clearance as sppkcl
from sp_pharkin import expo


def main() -> None:
    v = sppkcl.average_clearance_weight(
        average_clearance="0.04L/(hour*kilogram)", weight="75 kilogram"
    )
    print(v)

    expo.solve_for_k(c_0=100, c_t=50, t=20)


if __name__ == "__main__":
    main()
