from __future__ import print_function, division
import gnomad_utils
from optparse import OptionParser

'''
main
'''
if __name__ == "__main__":
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage=usage)

    parser.add_option("--config",
                      dest="config", default='config.json',
                      help="config everything there [default: %default]")
    (options, args) = parser.parse_args()

    # load options
    opt = json.load(open(options.config,'r'))
    
    r = report(opt)
    r.run()
