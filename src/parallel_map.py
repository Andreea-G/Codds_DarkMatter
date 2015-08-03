# Module providing the 'parmap' function for multiprocessing, which is used instead of
# multiprocessing.Pool.map.
# Original code written by klaus se
# see http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class

from __future__ import print_function
import multiprocessing


def fun(f, q_in, q_out, verbose):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        if verbose:
            print("process #", i)
        if isinstance(x, tuple):
            q_out.put((i, f(*x)))
        elif isinstance(x, dict):
            q_out.put((i, f(**x)))
        else:
            q_out.put((i, f(x)))


def parmap(f, X, processes=multiprocessing.cpu_count(), verbose=True):
    ''' Used instead of Pool.map function for parallel programming.
        Works with lambdas and with member functions of classes.
    Input:
        f: callable
            Function to be called.
        X: iterable
            Apply f to each element in X, collecting the results in a list
            which is returned.
        processes: int, optional
            Number of processes running at one time. If not given, it is set to
            the cpu_count of the computer.
        verbose: bool, optional
            If true, it prints the process # when it starts.
    Returns:
        res: list
            List of f(x) for each x in X.
    '''
    if processes is None:
        processes = multiprocessing.cpu_count()
    print("Number of processes = ", processes)
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out, verbose))
            for _ in range(processes)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(processes)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]


if __name__ == '__main__':
    import time

    args = [(i, i + 1) for i in range(8)]
    print(isinstance(args[0], tuple))

    def func(i, j):
        time.sleep(1)
        return 2 * i + j

    def wrapfunc(args):
        print(args)
        return func(*args)

    print(parmap(func, args))
    print("done")

    args1 = list(range(10))
    print(parmap(lambda x: x**2, args1))
    print("done")

    kwargs = [{'c': 1, 'b': 2, 'a': 3}] * 8

    def func2(a, b, e=0, c=0):
        return a + 10 * b + 100 * c + 1000 * e

    print(parmap(func2, kwargs))
