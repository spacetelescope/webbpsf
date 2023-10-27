import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "webbpsf: mark test as requiring webbpsf data to run"
    )


def pytest_addoption(parser):
    parser.addoption(
        "--webbpsf", action="store_true", default=False, help="run webbpsf tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--webbpsf"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_webbpsf = pytest.mark.skip(reason="need --webbpsf option to run")
    for item in items:
        if "webbpsf" in item.keywords:
            item.add_marker(skip_webbpsf)
