import resource
import platform
import sys
import time
import logging as log


def memory_limit(percentage: float):
    """
    只在linux操作系统起作用
    """
    if platform.system() != "Linux":
        print("Only works on linux!")
        return
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    limit = int(get_memory() * 1024 * percentage)
    print(limit, hard)
    resource.setrlimit(resource.RLIMIT_AS, (limit, hard))


def get_memory():
    with open("/proc/meminfo", "r") as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) == "MemAvailable:":
                free_memory = int(sline[1])
                break
    return free_memory


def memory(percentage=0.8):
    def decorator(function):
        def wrapper(*args, **kwargs):
            memory_limit(percentage)
            try:
                return function(*args, **kwargs)
            except MemoryError:
                mem = get_memory() / 1024 / 1024
                print("Remain: %.2f GB" % mem)
                sys.stderr.write("\n\nERROR: Memory Exception\n")
                sys.exit(1)

        return wrapper

    return decorator


def measure_time(func):
    def wrapper(*arg):
        t = time.time()
        res = func(*arg)
        log.DEBUG("Function took " + str(time.time() - t) + " seconds to run")
        return res

    return wrapper
