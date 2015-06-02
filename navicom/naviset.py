################################################################################
# naviset.py
# By Mathurin Dorel
# Manage a set of datasets in navicom objects and perform comparison and 
# linked displays with them
################################################################################

from .navicom import *
from copy import deepcopy

class NaviSet():
    """
    Manage a set of datasets in navicom objects and perform comparison and linked displays with them
    """

    def __init__(self, map_url='https://navicell.curie.fr/navicell/maps/cellcycle/master/index.php', browser_command="firefox %s", display_config = DisplayConfig()):
        """
        Initialise a NaviSet object.

        Args:
            display_config (DisplayConfig): The display configuration to apply to all objects managed by the NaviSet
        """
        self._ncs = dict()
        assert isinstance(display_config, DisplayConfig), "The type of 'display_config' must be DisplayConfig"
        self._display_config = display_config
        self._map_url = map_url
        self._browser_command = browser_command

    def addNaviCom(self, fname="", modules_dict="", name="no name"):
        """
        Create a new NaviCom object bound to the NaviSet. NaviCom objects created in this way can be called directly by the NaviSet object.
        """
        nc = NaviCom(self._map_url, fname, modules_dict, self._browser_command, self._display_config, name)
        if (nc.name in self._ncs):
            warn("Replacing NaviCom object at " + name)
        self._ncs[nc.name] = nc

    def bindNaviCom(self, navicom):
        """
        Add a new NaviCom object to manage. A deep copy of the object is created, so a new NaviCell window will be used.
        """
        if (navicom.name in self._ncs):
            warn("Replacing NaviCom object at " + name)
        self._ncs[navicom.name] = deepcopy(navicom)
        self._ncs[navicom.name]._display_config = self._display_config

    def compareDatasets(self, ds1, ds2):
        """
            Create a differential dataset, where the datatables are the difference ds2 - ds1 between the corresponding datatables.
        """
        nname = ds2.name + "-" + ds1.name
        self.addNaviCom(name=nname)
        #for processing in ds1._data:
            #for method in ds1._data[processing]:
                #if (method in ds2._data[processing]):
                    #nd = NaviData(ds2._data - ds1._data, )

