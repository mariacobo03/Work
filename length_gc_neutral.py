import sys
s = sys.stdin.readline().strip()


def gc_neutral(s):
    n = len(s)
    if n % 2 != 0:
        return False
    else:
        gc = 0
        for i in range(n):
            if s[i] == 'G' or s[i] == 'C':
                gc += 1
    return gc == n/2


print(gc_neutral(s))


def longest_subsequence(s):
    longest = 0
    for i in range(len(s) - 1):
        for j in range(i + 1, len(s)):
            if len(s[i:j + 1]) >= longest and gc_neutral(s[i:j + 1]):
                longest = len(s[i:j + 1])
    return longest


print(longest_subsequence(s))

