################################################################################
# displayConfig.py
# By Mathurin Dorel
# Configuration of the color gradients for NaviCell
################################################################################

def rgbValid(color):
    """
    Check whether a color string is a valid rgb color
    """
    if (len(color) != 6):
        if (len(color) == 7 and color[0]=="#"):
            color = color[1:]
        elif (len(color) == 8 and color[1]=="x"):
            color = color[2:]
        else:
            raise ValueError(str(color) + " is not a valid RGB color (does not comply to the format (#|0x)?RRGGBB)")
    for cc in color:
        if (not cc.lower() in [str(el) for el in range(10)]+[chr(aa) for aa in range(ord("a"), ord("f")+1)]):
            raise ValueError(str(color) + " is not a valid RGB color (" + str(cc) + " is not an hexadecimal digit)")
    return(color)

def getGradient(color0, color1, steps):
    """
    Create a gradient between two colors.
    Args:
        color0 (str): first RGB color of the gradient
        color1 (str): second RGB color of the gradient
        steps (int): number of colors in the final gradient
    Return:
        gradient (str): a list of RGB colors (str)
    """
    assert steps >= 2, ValueError("Cannot make a gradient with only one value")
    def hex2(nb):
        nb = hex(nb)[2:]
        if (len(nb) == 1):
            nb = "0"+nb
        return(nb)
    c0 = [int(color0[ii:ii+2], 16) for ii in range(0, 6, 2)]
    c1 = [int(color1[ii:ii+2], 16) for ii in range(0, 6, 2)]
    gradient=list()
    for ii in range(steps):
        ccolor = [c0[jj]+ii*(c1[jj]-c0[jj])//(steps-1) for jj in range(3)]
        ccolor = "".join([hex2(el) for el in ccolor])
        gradient.append(ccolor)
    return(gradient)

class DisplayConfig():
    """
    DisplayConfig class to set the color gradients configuration in NaviCell
    """

    def __init__(self, step_count=3, color_gradient=["00FF00", "FF0000"], zero_color="ffffff", na_color="ffffff", zero_size=0, na_size=0, use_absolute_values = False):
        """
        Initialise a color gradient configuration
        Args:
            step_count (int): number of steps for the color gradients. A step for NAs is automatically attributed.
            color_gradient (list): a list of colors of length 2 or step_count. If length 2 a gradient is built, if the length is step_count the list is used for the colors.
            zero_color (str): an hexadecimal string for the color of the zero, only visible if step_count is odd
        """
        assert isinstance(step_count, int), ValueError("'step_count' must be an integer")
        assert isinstance(color_gradient, list) and len(color_gradient)>=2, ValueError("'color_gradient' must be list of 2 or more elements")
        for ii in range(len(color_gradient)):
            color_gradient[ii] = rgbValid(color_gradient[ii])
        assert isinstance(zero_color, str), ValueError("'zero_color' must be a string")
        zero_color = rgbValid(zero_color)
        assert isinstance(na_color, str) and len(na_color) > 3, ValueError("'na_color' must be a string of at least 3 characters representing a RGB color")
        self.na_color = rgbValid(na_color)
        self.zero_color = zero_color
        self.na_size = na_size
        self.zero_size = zero_size

        self.use_absolute_values = use_absolute_values

        self.step_count = step_count
        if (len(color_gradient) == step_count):
            self.colors = color_gradient
        elif (len(color_gradient) == 2):
            self.color = list()
            if (zero_color == ""):
                self.colors = getGradient(color_gradient[0], color_gradient[1], self.step_count)
            else:
                self.colors = getGradient(color_gradient[0], zero_color, self.step_count//2+1)[:-1]
                if (self.step_count%2 == 1):
                    self.colors += [zero_color]
                self.colors += getGradient(zero_color, color_gradient[1], self.step_count//2+1)[1:]
        elif (len(color_gradient) == step_count-1 and zero_color != ""):
            self.colors = color_gradient[:step_count//2+1]
            if (step_count%2 == 1):
                self.colors += [zero_color]
            self.colors += color_gradient[step_count//2+1:]
        else:
            raise ValueError("The length of 'color_gradient' must be 2 or equal to step_count")

    def __repr__(self):
        rpr = "DisplayConfig object :\n"
        rpr+= "\tGradient colors: " + str(self.colors) + "\n"
        rpr+= "\tNA color: " + str(self.na_color)
        return(rpr)

