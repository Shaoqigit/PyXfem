import sys

from SAcouS.interface.sol_setup import PyAcousiXSetuper


def main():
    file_path = sys.argv[2]
    sol_setuper = PyAcousiXSetuper(file_path)
    sol_setuper.welcome()