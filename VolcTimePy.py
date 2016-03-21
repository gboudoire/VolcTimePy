#!/usr/bin/env python
# encoding: utf-8
'''
VolcTimePy -- Timescales calculation from diffusion profiles in zoned crystals.

VolcTimePy is a timescales calculation from diffusion profiles in zoned
crystals.

@author:     Guillaume Boudoire (LGSR/OVPF/IPGP)

@author:     Patrice Boissier (OVPF/IPGP)

@copyright:      Guillaume Boudoire (guillaume.boudoire@gmail.com)
                    Laboratoire GeoSciences Réunion
                    Université de La Reunion
                    Observatoire Volcanologique du Piton de La Fournaise
                    Institut de Physique du Globe de Paris
                Patrice Boissier (boissier@ipgp.fr)
                    Observatoire Volcanologique du Piton de La Fournaise
                    Institut de Physique du Globe de Paris

@license:    GNU Lesser General Public License, Version 3
            (http://www.gnu.org/copyleft/lesser.html)

@contact:    guillaume.boudoire@gmail.com
'''

import sys
import os

from optparse import OptionParser

__all__ = []
__version__ = 0.1
__date__ = '2016-03-21'
__updated__ = '2016-03-21'

DEBUG = 1
TESTRUN = 0
PROFILE = 0


def main(argv=None):
    '''Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v0.1"
    program_build_date = "%s" % __updated__

    program_version_string = '%%prog %s (%s)' % (program_version,
                                                 program_build_date)
    program_longdesc = ''''''
    program_license = "Copyright 2016 Boudoire / Boissier (OVPF/IPGP)\n\
                Licensed under the GNU Lesser General Public License, V3"

    if argv is None:
        argv = sys.argv[1:]
    try:
        # setup option parser
        parser = OptionParser(version=program_version_string,
                              epilog=program_longdesc,
                              description=program_license)
        parser.add_option("-i", "--in",
                          dest="infile",
                          help="set input path [default: %default]",
                          metavar="FILE")
        parser.add_option("-o", "--out",
                          dest="outfile",
                          help="set output path [default: %default]",
                          metavar="FILE")
        parser.add_option("-v", "--verbose",
                          dest="verbose",
                          action="count",
                          help="set verbosity level [default: %default]")

        # set defaults
        parser.set_defaults(outfile="./out.txt", infile="./in.txt")

        # process options
        (opts, args) = parser.parse_args(argv)

        if opts.verbose > 0:
            print("verbosity level = %d" % opts.verbose)
        if opts.infile:
            print("infile = %s" % opts.infile)
        if opts.outfile:
            print("outfile = %s" % opts.outfile)

        # MAIN BODY #

    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'VolcTimePy_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
