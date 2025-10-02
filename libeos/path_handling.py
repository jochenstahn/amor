"""
Defines how file paths are resolved from short_notation, year and number to filename.
"""
import os
from typing import List


class PathResolver:
    def __init__(self, year, rawPath):
        self.year = year
        self.rawPath = rawPath

    def resolve(self, short_notation):
        return list(map(self.get_path, self.expand_file_list(short_notation)))

    @staticmethod
    def expand_file_list(short_notation)->List[int]:
        """Evaluate string entry for file number lists"""
        file_list = []
        for i in short_notation.split(','):
            if '-' in i:
                if ':' in i:
                    step = i.split(':', 1)[1]
                    file_list += range(int(i.split('-', 1)[0]),
                                       int((i.rsplit('-', 1)[1]).split(':', 1)[0])+1,
                                       int(step))
                else:
                    step = 1
                    file_list += range(int(i.split('-', 1)[0]),
                                       int(i.split('-', 1)[1])+1,
                                       int(step))
            else:
                file_list += [int(i)]
        return list(sorted(file_list))

    def get_path(self, number):
        fileName = f'amor{self.year}n{number:06d}.hdf'
        path = ''
        for rawd in self.rawPath:
            if os.path.exists(os.path.join(rawd, fileName)):
                path = rawd
                break
        if not path:
            if os.path.exists(
                    f'/afs/psi.ch/project/sinqdata/{self.year}/amor/{int(number/1000)}/{fileName}'):
                path = f'/afs/psi.ch/project/sinqdata/{self.year}/amor/{int(number/1000)}'
            else:
                raise FileNotFoundError(f'# ERROR: the file {fileName} can not be found in {self.rawPath}')
        return os.path.join(path, fileName)
