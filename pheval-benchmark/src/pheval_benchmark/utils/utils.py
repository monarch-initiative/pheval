import os
from pheval_benchmark.custom_exceptions import IncorrectFileFormatError


class DirectoryFiles:
    """ Class that retrieves an ordered list of relevant files from a directory"""

    def __init__(self, directory: str, file_suffix: str):
        self.directory = directory
        self.file_suffix = file_suffix

    def obtain_all_files(self) -> list:
        file_list = []
        file_directory = os.fsdecode(os.path.join(self.directory, ''))
        for file in os.listdir(file_directory):
            filename = os.fsdecode(file)
            file_list.append(os.path.join(file_directory, filename))
        file_list.sort()
        return file_list

    def obtain_files_suffix(self) -> list:
        file_list = []
        directory_ = os.fsencode(os.path.join(self.directory, ''))
        for file in os.listdir(directory_):
            filename = os.fsdecode(file)
            if filename.endswith(self.file_suffix):
                file_list.append(filename)
        file_list.sort()
        return file_list

    def obtain_files_full_path(self) -> list:
        file_list = []
        file_directory = os.fsdecode(os.path.join(self.directory, ''))
        for file in os.listdir(file_directory):
            filename = os.fsdecode(file)
            if not filename.endswith(self.file_suffix):
                raise IncorrectFileFormatError(filename, self.file_suffix)
            file_list.append(os.path.join(file_directory, filename))
        file_list.sort()
        return file_list
