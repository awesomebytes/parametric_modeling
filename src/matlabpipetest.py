from matlabpipe import MatlabPipe


def test_t_tide_demo(matlab):
    """Call Matlab and run t_tide demos."""
    cmd = """[b,a] = butter(6,0.2);              % Butterworth filter design
    h = filter(b,a,[1 zeros(1,100)]);   % Filter data using above filter
    [bb,aa] = stmcb(h,4,4);
    """
    out = matlab.eval(cmd)
    print "Printing output"
    print out
    print "bb is:"
    bb = matlab.get('bb')
    print bb
    del(out)

if __name__ == '__main__':
    matlab = MatlabPipe(matlab_version='2013a')
    matlab.open()
    test_t_tide_demo(matlab)
    
    