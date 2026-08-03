

def makewindows(start, end, windowsize):
    """create vectors of specified windows"""
    _end = min(end, windowsize)
    starts = [start]
    ends = [_end]
    while _end < end:
        _end = min(_end + windowsize, end)
        ends.append(_end)
        start += windowsize
        starts.append(start)
    return starts, ends

