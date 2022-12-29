"""This module handles zone data from `webnucleo <https://webnucleo.readthedocs.io>`_ files."""

import wnutils.xml as wx


class Zones_Xml:
    """A class for handling webnucleo zones.

    Args:
        ``file`` (:obj:`str`):  A string giving the XML file name with the zone data.

    """

    def __init__(self, file):
        self.xml = wx.Xml(file)
        self.zones = self.xml.get_zone_data()

    def get_zones(self, zone_xpath=""):
        """Method to return zones.

        Args:
            ``zone_xpath`` (:obj:`str`, optional):  An XPath expression to select zones.  Default is all zones.

        Returns:
            A :obj:`dict` of `wnutils <https://wnutils.readthedocs.io>`_ zone data objects.

        """

        if not zone_xpath:
            return self.zones
        else:
            return self.xml.get_zone_data(zone_xpath=zone_xpath)
