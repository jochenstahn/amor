"""
Defines how file paths are resolved from short_notation, year and number to filename.
"""
import logging
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
            if '-' in i and not i.startswith('-'):
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
        if number<=0:
            number = self.search_latest(number)
        fileName = f'amor{self.year}n{number:06d}.hdf'
        path = ''
        for rawd in self.rawPath:
            if os.path.exists(os.path.join(rawd, fileName)):
                path = rawd
                break
        if not path:
            from glob import glob
            potential_file = glob(f'/home/amor/data/{self.year}/*/{fileName}')
            if len(potential_file)>0:
                path = os.path.dirname(potential_file[0])
            else:
                raise FileNotFoundError(f'# ERROR: the file {fileName} can not be found '
                                        f'in {self.rawPath+["/home/amor/data"]}')
        return os.path.join(path, fileName)

    def search_latest(self, number):
        if number>0:
            raise ValueError('number needs to be relative index (negative)')
        if os.path.exists(f'/home/amor/data/{self.year}/DataNumber'):
            try:
                with open(f'/home/amor/data/{self.year}/DataNumber', 'r') as fh:
                    current_index = int(fh.readline())-1
            except Exception:
                logging.error('Can not access DataNumber', exc_info=True)
            else:
                return current_index+number
        # find all files from the given year, convert to number and return latest
        from glob import glob
        possible_files = []
        for rawd in self.rawPath:
            possible_files += glob(os.path.join(rawd, f'amor{self.year}n??????.hdf'))

        possible_files += glob(f'/home/amor/data/{self.year}/*/amor{self.year}n??????.hdf')
        possible_indices = list(set([int(os.path.basename(fi)[9:15]) for fi in possible_files]))
        possible_indices.sort()
        try:
            return possible_indices[number-1]
        except IndexError:
            raise FileNotFoundError(f'# Could not find suitable file for relative index {number} '
                                    f'in {self.rawPath+["/home/amor/data"]}, '
                                    f'possible indices {possible_indices}')

