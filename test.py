import ete3

try:
    from ete3 import TreeStyle
    print("TreeStyle imported successfully!")
except ImportError as e:
    print(f"ImportError: {e}")