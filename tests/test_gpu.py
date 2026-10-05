"""GPU smoke tests: OpenMM must run on CUDA and PyTorch must see the RTX 3090."""
import pytest

pytestmark = pytest.mark.gpu


def test_openmm_cuda_runs():
    import openmm as mm
    from openmm import unit

    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    assert "CUDA" in names, f"OpenMM CUDA platform missing; available: {names}"

    # Tiny Lennard-Jones-free system: free particles under a Langevin integrator.
    system = mm.System()
    n = 64
    for _ in range(n):
        system.addParticle(1.0)
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = mm.Platform.getPlatformByName("CUDA")
    ctx = mm.Context(system, integrator, platform, {"Precision": "mixed"})
    assert ctx.getPlatform().getName() == "CUDA"

    import numpy as np
    rng = np.random.default_rng(0)
    ctx.setPositions(rng.random((n, 3)))
    ctx.setVelocitiesToTemperature(300 * unit.kelvin, 1)
    integrator.step(100)
    pos = ctx.getState(getPositions=True).getPositions(asNumpy=True)
    assert np.isfinite(pos.value_in_unit(unit.nanometer)).all()

    # Confirm it is the 3090 and not some other CUDA device.
    device = platform.getPropertyValue(ctx, "DeviceIndex")
    print(f"OpenMM CUDA device index: {device}, precision: {platform.getPropertyValue(ctx, 'Precision')}")


def test_polychrom_importable():
    import polychrom  # noqa: F401
    from polychrom import simulation, forces  # noqa: F401


def test_torch_sees_gpu():
    import torch

    assert torch.cuda.is_available(), "torch.cuda.is_available() is False (CPU-only torch build?)"
    name = torch.cuda.get_device_name(0)
    assert "3090" in name, f"Expected RTX 3090, got {name}"
    x = torch.randn(1024, 1024, device="cuda")
    assert torch.isfinite(x @ x).all()


def test_cupy_sees_gpu():
    import cupy as cp

    assert cp.cuda.runtime.getDeviceCount() >= 1
    a = cp.arange(10, dtype=cp.float32)
    assert float(a.sum().get()) == 45.0


def test_trackastra_importable():
    import trackastra  # noqa: F401
    from trackastra.model import Trackastra  # noqa: F401
