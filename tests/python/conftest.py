import os, pytest

def pytest_addoption(parser):
    parser.addoption(
        "--datadir", action="store", default="tests/data", help="location of data dir"
    )

def pytest_sessionstart(session):
    os.environ['VRNA_TEST_DATA'] = session.config.getoption("--datadir")
