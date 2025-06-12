import argparse

from matchmaps._args import (
    matchmaps_description,
    matchmaps_mr_description,
    matchmaps_ncs_description,
    matchmaps_diagnose_description,
    base_and_mr_args,
    ncs_args,
    common_matchmaps_args,
    global_args,
    diagnose_args,
)

matchmaps_parser = argparse.ArgumentParser(description=matchmaps_description)
matchmaps_mr_parser = argparse.ArgumentParser(description=matchmaps_mr_description)
matchmaps_ncs_parser = argparse.ArgumentParser(description=matchmaps_ncs_description)

matchmaps_diagnose_parser = argparse.ArgumentParser(
    description=matchmaps_diagnose_description
)


for args, kwargs in base_and_mr_args:
    matchmaps_parser.add_argument(*args, **kwargs)
    matchmaps_mr_parser.add_argument(*args, **kwargs)

for args, kwargs in ncs_args:
    matchmaps_ncs_parser.add_argument(*args, **kwargs)

for args, kwargs in diagnose_args:
    matchmaps_diagnose_parser.add_argument(*args, **kwargs)

for args, kwargs in common_matchmaps_args:
    matchmaps_parser.add_argument(*args, **kwargs)
    matchmaps_mr_parser.add_argument(*args, **kwargs)
    matchmaps_ncs_parser.add_argument(*args, **kwargs)

for args, kwargs in global_args:
    matchmaps_parser.add_argument(*args, **kwargs)
    matchmaps_mr_parser.add_argument(*args, **kwargs)
    matchmaps_ncs_parser.add_argument(*args, **kwargs)
    matchmaps_diagnose_parser.add_argument(*args, **kwargs)
