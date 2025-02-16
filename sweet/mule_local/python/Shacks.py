import os, sys


class Shacks:
    def __init__(self, i_directory: str):
        """
        First, we will setup a list of all shack modules in the directory 'shacksShared'
        """
        import pkgutil
        from importlib.machinery import SourceFileLoader
        
        pkgpath = os.path.dirname(__file__)+"/"+i_directory
        shackFiles = [shackFile for _, shackFile, _ in pkgutil.iter_modules([pkgpath])]

        self.shacks = {}
        for shackFile in shackFiles:
            if shackFile[0] == '_':
                continue
            
            fullpath = pkgpath+"/"+shackFile+".py"
            module = SourceFileLoader("",fullpath).load_module()
            
            # Strip leading numbers and "_"
            shackName = shackFile[:]
            while shackName[0] in "0123456789_":
                shackName = shackName[1:]
            
            self.shacks[shackFile] = module.__dict__[shackName]


    def keys(self):
        return self.shacks.keys()

    def values(self):
        return self.shacks.values()
