"""Make unbiased difference maps even for non-isomorphous inputs"""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("matchmaps")
except PackageNotFoundError:
    __version__ = "uninstalled"

__author__ = "Dennis Brookner"
__email__ = "debrookner@gmail.com"

from matchmaps._utils import make_floatgrid_from_mtz, rigid_body_refinement_wrapper, align_grids_from_model_transform
