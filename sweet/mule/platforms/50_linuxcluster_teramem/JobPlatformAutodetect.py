
import platform
import os

def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    #
    # There's no autodetection for this platform since this is a special type on the linux cluster
    #

    #if platform.node()[:8] == 'cm2login':
    #	return True

    return False


if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

