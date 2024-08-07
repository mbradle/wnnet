"""This module handles zone data from `webnucleo <https://webnucleo.readthedocs.io>`_ files."""

import wnutils.xml as wx


class Zones:
    """A class for handling webnucleo zones.

    Args:
        ``file`` (:obj:`str`):  A string giving the XML file name with the
        zone data.

    """

    def __init__(self):
        self.xml = None

    def store_data_from_xml(self, xml_file):
        """Method to retrieve zone data from an xml file.

        Args:
            ``xml_file`` (:obj:`str`):  The name of the XML file with zone
            data.

        Returns:
            On successful return, the data have been stored.

        """
        self.xml = wx.Xml(xml_file)

    def get_zones(self, zone_xpath=""):
        """Method to return zones.

        Args:
            ``zone_xpath`` (:obj:`str`, optional):  An XPath expression to
            select zones.  Default is all zones.

        Returns:
            A :obj:`dict` of `wnutils <https://wnutils.readthedocs.io>`_ zone
            data objects.

        """

        return self.xml.get_zone_data(zone_xpath=zone_xpath)
