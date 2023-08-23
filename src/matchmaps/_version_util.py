"""matchmaps.version returns package version"""

from importlib.metadata import version

def main():
    print(version("matchmaps"))
    return

if __name__ == "__ main__": 
    main()
