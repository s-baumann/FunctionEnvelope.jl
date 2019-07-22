# FunctionEnvelope

This is a small helper package that efficiently maps out the upper envelope of a set of nondecreasing functions.

This is incredibly rough. I do not know of any established algorithms for mapping out the upper envelope of a set of functions whilst minimising the number of function evaluations. Thus I had to invent one and it is rough as I have an application and am not really looking at this problem for its own sake. I only consider 1d and continuous functions. The algorithm is basically:
* Evaluate each function at each endpoint of an interval.
* Then we move left slightly to some point x. We suspect some function g is the highest as it is the highest function at a neighbouring point. We evaluate g and get a value of g'. We then look at all of the other functions and see how high they could potentially be. For increasing functions they cannot be higher than they are in the nearest value we have for them to the right. For a decreasing function they cannot be higher than the nearest value we have for them to the left. For a constant function we know what the constant is. If we cannot say what type of function it is we have to evaluate at every point.
* Then we keep applying this logic a bit to the left repeatedly until we have covered the interval.
* We now have the highest functions at each point in a grid covering the space. We just need to solve for the exact locations of any switchovers between two functions. We can do this with a bisection algorithm.

Of course we will be getting functions with big but narrow high points but there is no way I can think to avoid this (and in the economics problems that I care about this should not happen).

This package does not yet support lower envelopes. It should be relatively easy to do however by flipping each function and constant.

So it is pretty rough but it seems to work. If anyone comes across a better algorithm I would interested to hear of it.


# TODO

I reckon I could get a 20%ish speed boost in the cases where we know functions are globally concave or globally convex if we get the first derivative at a point to the left and use this to get a more accurate functions_potentially_above_mark when we are crawling left looking for the highest function at each grid point.
