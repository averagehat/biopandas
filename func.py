from functools import partial
from itertools import izip
import string
#notin = compose(_not, operator.methodcaller('__contains__'))
#notin = compose(_not, attr('__contains__'))
#mismatches = pfilter(notin('M='))

def _not(x):
    return not x
def partial2(method, param):
      def t(x):
              return method(x, param)
      return t

def _id(x): return x

def apply_key_func(k, v, funcdict):
    return funcdict.get(k, _id)(v)
def compose_all(*funcs):
    return reduce(compose, funcs)

#def k_compose(outer, **okwargs):
#    ''' compose(f, g)(x) == f(g(x)) '''
#    def newfunc(*args, **ikwargs):
#        _kwargs = dict( (k, apply_key_func(k, v, okwargs)) for k, v in ikwargs.items())
#        return outer(*args, **_kwargs)
#    return newfunc

def compose(outer, inner):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc

def fzip(funcs, args):
    for func, arg in izip(funcs, args):
        yield func(arg)

def dictmap(func, _dict):
    return dict( (key, func(val)) for key, val in _dict.items())

def reverse(collection): return collection[::-1]

compose_list = partial(reduce, compose)
compose_all = compose(compose_list, lambda *a : a)
pmap = partial(partial, map)
pfilter = partial(partial, filter)
#TODO: could use partial2 instead
pstrip = lambda x: partial(string.split, chars=x)
psplit = lambda x: partial(string.split, sep=x)
