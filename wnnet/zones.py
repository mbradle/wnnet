"""This module handles zone data from `webnucleo <https://webnucleo.readthedocs.io>`_ files."""

import wnutils.xml as wx


class Zones:
    """A class for handling webnucleo zones.

    Args:
        ``file`` (:obj:`str`):  A string giving the XML file name with the\
         zone data.

    """

    def __init__(self, file):
        self.xml = wx.Xml(file)

    def get_zones(self, zone_xpath=""):
        """Method to return zones.

        Args:
            ``zone_xpath`` (:obj:`str`, optional):  An XPath expression to\
             select zones.  Default is all zones.

        Returns:
            A :obj:`dict` of `wnutils <https://wnutils.readthedocs.io>`_\
             zone data objects.

        """

        return self.xml.get_zone_data(zone_xpath=zone_xpath)

    def write(self, file, pretty_print=True, zone_xpath=""):
        """Method to write the zones to an xml file.

        Args:

           ``file`` (:obj:`str`): A string giving the name of output\
            xml file.

           ``pretty_print`` (:obj:`bool`, optional): If set to True,\
           routine outputs the xml in nice indented format.

           ``zone_xpath`` (:obj:`str`, optional):  An XPath expression to\
             select zones.  Default is all zones.

        Returns:
            On successful return, the underlying xml has been written\
            to ``file``.

        """

        zones = self.get_zones(zone_xpath=zone_xpath)

        new_xml = wx.New_Xml("zone_data")

        new_xml.set_zone_data(zones)
        new_xml.write(file, pretty_print=pretty_print)
