from pytokr import items

class Node:
    def __init__(self, value):
        self.value = value
        self.left = None
        self.right = None

# Insert new node into corresponding BST root
def insert(root, new_node):
    # If root is empty
    if not root:
        return Node(new_node)
    if new_node < root.value:
        root.left = insert(root.left, new_node)
    elif new_node > root.value:
        root.right = insert(root.right, new_node)
    return root

# Search for a value in the BST
def search(root, target):
    if not root:
        return 0  # Value not found
    if root.value == target:
        return 1  # Value found
    elif target < root.value:
        return search(root.left, target)
    else:
        return search(root.right, target)

# Read input items
def items():
    while True:
        line = input().strip()
        if line:
            yield line
        else:
            break

# Read description of the BST
_ = int(input())  # Ignore the number of internal nodes

# Read the description of the BST as a list of integers
tree_desc = list(map(int, input().split()))

# Construct the BST
root = None
for value in tree_desc:
    root = insert(root, value)

# Read natural numbers and check if they are in the tree
while True:
    try:
        natural = int(input())
        print(f"{natural} {search(root, natural)}")
    except EOFError:
        break
