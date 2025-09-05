"""matchmaps.version returns package versions"""
import subprocess
import os
from importlib.metadata import version, PackageNotFoundError
import time

from matchmaps._utils import _detect_phenix_version

def main():

    print( "========================")
    print(f"matchmaps.version run on {time.strftime('%c')}")
    print( "========================")

    print()
    print(f"versions of matchmaps and relevant python dependencies:")
    print( "-------------------------------------------------------")
    print()

    packages = ("matchmaps", "gemmi", "reciprocalspaceship", "rs-booster", "pandas", "numpy")
    longest_package_name = max(len(package) for package in packages)

    for package in packages:
        try:
            print(f"{package:>{longest_package_name}}: v{version(package)}")
        except PackageNotFoundError:
            print(package)
            raise PackageNotFoundError(f"{package}, which is weird, because it is a dependency of matchmaps!")

    print()
    print(f"external dependencies:")
    print( "----------------------")
    print()

    try:
        phenix = f"phenix found in environment at: {os.environ['PHENIX']}"

    except KeyError:
        phenix = "phenix not found in environment;\n       you will need this to run matchmaps, matchmaps.mr, or matchmaps.ncs"

    try:
        ccp4 = f"  ccp4 found in environment at: {os.environ['CCP4']}"

    except KeyError:
        ccp4 = "ccp4 not found in environment;\n       you will need this to run matchmaps;\n       matchmaps.mr and matchmaps.ncs do not require ccp4"

    print(phenix)
    print(ccp4)

    return


if __name__ == "__ main__":
    main()
