#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse


__author__ = "Joerg Weingrill"
__copyright__ = "Copyright 2021 Leibniz-Insitute for Astrophysics Potsdam (AIP)"
__credits__ = ["Joerg Weingrill"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Joerg Weingrill"
__email__ = "jweingrill@aip.de"
__status__ = "Development"
__date__ = "12/1/21"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='version counter')
    parser.add_argument('-b', '--build', action='store_true', help='increase build')
    parser.add_argument('-r', '--release', action='store_true', help='increase release')
    parser.add_argument('-v', '--version', action='store_true', help='increase version')

    parser.add_argument('filename', help='python file to pase for version string')

    args = parser.parse_args()
    lines = ''
    with open(args.filename, "rt") as file:
        for line in file:
            if line.startswith('__version__'):
                _, separator, versionstring = line.rpartition('=')
                versionstring = versionstring.rstrip("\n").strip('" ')
                version, release, build = versionstring.split('.')
                if args.build:
                    build = int(build) + 1
                if args.release:
                    release = int(release) + 1
                if args.version:
                    version = int(version) + 1
                new_version_string = f'__version__ = "{version}.{release}.{build}"'
                print(new_version_string)
                lines = lines + new_version_string
            else:
                print(f'{line}', end='')
                lines = lines + line
    with open(args.filename, "wt") as file:
        file.writelines(lines)
