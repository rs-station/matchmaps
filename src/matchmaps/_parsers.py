import argparse

from matchmaps._args import (
    matchmaps_description,
    matchmaps_mr_description,
    matchmaps_ncs_description,
    base_and_mr_args,
    ncs_args,
    common_args
)

matchmaps_parser = argparse.ArgumentParser(description=matchmaps_description)
matchmaps_mr_parser = argparse.ArgumentParser(description=matchmaps_mr_description)
matchmaps_ncs_parser = argparse.ArgumentParser(description=matchmaps_ncs_description)

for args, kwargs in base_and_mr_args:
    matchmaps_parser.add_argument(*args, **kwargs)
    matchmaps_mr_parser.add_argument(*args, **kwargs)

for args, kwargs in ncs_args:
    matchmaps_ncs_parser.add_argument(*args, **kwargs)

for args, kwargs in common_args:
    matchmaps_parser.add_argument(*args, **kwargs)
    matchmaps_mr_parser.add_argument(*args, **kwargs)
    matchmaps_ncs_parser.add_argument(*args, **kwargs)

