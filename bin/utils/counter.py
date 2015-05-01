from time import time
counters = {}

def counter_inc(flag):
    if not (flag in counters):
        counters[flag] = 1
    else:
        counters[flag] += 1

def counter_max(flag, value):
    if not (flag in counters):
        counters[flag] = value
    elif counters[flag] < value:
        counters[flag] = value

def counter_print_counters():
    for flag in counters:
        print '[COUNTER] %s: %d' % (flag, counters[flag])


