"""Roof-embedded PCM study (Abass & Muthulingam 2024): OpenFOAM case generation for a
concrete roof slab unit cell. Standard library + PyYAML only, so it runs on any HPC
login or compute node.

    python3 -m roofpipe baseline     Phase 1: pure concrete slab case
    python3 -m roofpipe check RUN    Phase 1 gates: periodicity, 1D reference, energy
"""
__version__ = "0.1.0"
