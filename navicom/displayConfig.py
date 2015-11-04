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

def hex2(nb):
    nb = hex(nb)[2:]
    if (len(nb) == 1):
        nb = "0"+nb
    return(nb)

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
    c0 = [int(color0[ii:ii+2], 16) for ii in range(0, 6, 2)]
    c1 = [int(color1[ii:ii+2], 16) for ii in range(0, 6, 2)]
    gradient=list()
    for ii in range(steps):
        ccolor = [c0[jj]+ii*(c1[jj]-c0[jj])//(steps-1) for jj in range(3)]
        ccolor = "".join([hex2(el) for el in ccolor])
        gradient.append(ccolor)
    return(gradient)

def numberToColor(color_numbers):
    """
        Convert a triplet of base16 int to a color string

        Args:
            color_numbers (list): 3 base16 integers representing red, green and blue colors
        Return:
            color_string (str): A string representing the color
    """
    if (isinstance(color_numbers, str)): return(color_numbers)
    return( "".join([hex2(el) for el in color_numbers]) )

def colorToNumber(color_string):
    """
        Convert a color string to a base16 triplet

        Args:
            color_string (str): A string representing the color
        Return:
            color_numbers (list): 3 base16 integers representing red, green and blue colors
    """
    if (isinstance(color_string, list)): return(color_string)
    return( [int(color_string[ii:ii+2], 16) for ii in range(0, 6, 2)] )

def addColors(c0, c1, ratio=0.5):
    """
        Mix two colors, like in painting

        Args:
            c0 (str): The first color
            c1 (str): The second color
            ratio (float): Proportion of the first color in the mix
    """
    assert ratio >= 0 and ratio <= 1, ValueError("Ratio should be between 0 and 1")
    c0 = colorToNumber(c0)
    c1 = colorToNumber(c1)
    cc = [int(ratio*c0[ii]+(1-ratio)*c1[ii]) for ii in range(3)]
    return(numberToColor(cc))


class DisplayConfig():
    """
    DisplayConfig class to set the color gradients configuration in NaviCell
    """

    def __init__(self, step_count=5, color_gradient=["00FF00", "FF0000"], zero_color="ffffff", na_color="ffffff", uniform_color="0000ff", zero_size=0, na_size=0, use_absolute_values = False, excluded="10%", groups_sharpening=10):
        """
        Initialise a color gradient configuration
        Args:
            step_count (int): number of steps for the color gradients. A step for NAs is automatically attributed.
            color_gradient (list): a list of colors of length 2 or step_count. If length 2 a gradient is built, if the length is step_count the list is used for the colors.
            zero_color (str): an hexadecimal string for the color of the zero, only visible if step_count is odd
            zero_size (int): size of the glyphs for values close to 0
            excluded (str or float): percentage of values to exclude for the color gradient, usefull to avoid a distortion of the colorscale by the extreme values. If str, must match '\d\d?%'.
            groups_sharpening (int or float): intensity of shrinkage of the color scales for groups compared to samples
        """
        assert isinstance(step_count, int), ValueError("'step_count' must be an integer")
        assert isinstance(color_gradient, list) and len(color_gradient)>=2, ValueError("'color_gradient' must be list of 2 or more elements")
        for ii in range(len(color_gradient)):
            color_gradient[ii] = rgbValid(color_gradient[ii])
        assert isinstance(zero_color, str), ValueError("'zero_color' must be a string")
        zero_color = rgbValid(zero_color)
        assert isinstance(na_color, str) and len(na_color) > 3, ValueError("'na_color' must be a string of at least 3 characters representing a RGB color")
        self.na_color = rgbValid(na_color)
        self._zero_color = zero_color
        self.na_size = na_size
        self.zero_size = zero_size
        self.uniform_color = rgbValid(uniform_color)
        
        assert isinstance(groups_sharpening, int) or isinstance(groups_sharpening, float)
        self._groups_sharpening = groups_sharpening

        self.use_absolute_values = use_absolute_values

        self.step_count = step_count
        if (len(color_gradient) == step_count):
            self._colors = color_gradient
        elif (len(color_gradient) == 2):
            self._colors = list()
            if (zero_color == ""):
                self._colors = getGradient(color_gradient[0], color_gradient[1], self.step_count)
            else:
                self._colors = getGradient(color_gradient[0], zero_color, self.step_count//2+1)[:-1]
                if (self.step_count%2 == 1):
                    self._colors += [zero_color]
                self._colors += getGradient(zero_color, color_gradient[1], self.step_count//2+1)[1:]
        elif (len(color_gradient) == step_count-1 and zero_color != ""):
            self._colors = color_gradient[:step_count//2+1]
            if (step_count%2 == 1):
                self._colors += [zero_color]
            self._colors += color_gradient[step_count//2+1:]
        else:
            raise ValueError("The length of 'color_gradient' must be 2 or equal to step_count")

        # Fraction of excluded values
        if (isinstance(excluded, str)):
            assert excluded[-1] == "%", "Invalid percentage format"
            self._excluded = float(excluded[:-1]) / 100
        elif (isinstance(excluded, float)):
            self._excluded = excluded
        else:
            raise TypeError("'excluded' must be a float or a str (X%)")
        assert self._excluded < 1., "Cannot exclude more than 100% of the extreme values"

    def __repr__(self):
        rpr = "DisplayConfig object :\n"
        rpr+= "\tGradient colors: " + str(self._colors) + "\n"
        rpr+= "\tNA color: " + str(self.na_color)
        return(rpr)

SHAPE_ID = dict()
SHAPE_ID["triangle"] = 0
SHAPE_ID["square"] = 1
SHAPE_ID["rectangle"] = 2
SHAPE_ID["diamond"] = 3
SHAPE_ID["hexagon"] = 4
SHAPE_ID["circle"] = 5
class GlyphConfig():
    """
        GlyphConfig class to specify the color and shape configuration of glyphs for continuous data
    """

    def __init__(self, color="0000ff", shape="triangle", min_size=0, na_size=0):
        """
            Initialise a color and a shape for glyph that will be used for all values of a continuous data

            Args:
                color (str): Color to use for the glyph
                shape (str): Shape to use for the glyph
                min_size (int): Size of the glyph for the smallest value
                na_color (int): Size of the glyph for the NA value
        """
        assert isinstance(min_size, int), ValueError("'min_size' must be an integer")
        assert isinstance(na_size, int), ValueError("'na_size' must be an integer")
        assert isinstance(color, str), ValueError("'color' must be an string")
        assert isinstance(shape, str), ValueError("'shape' must be an string")
        self.min_size = min_size
        self.na_size = na_size
        self.color = color
        if (shape in SHAPE_ID):
            self.shape = SHAPE_ID[shape]
        else:
            warn("'" + shape + "' is not a valid shape")
            self.shape = 0

