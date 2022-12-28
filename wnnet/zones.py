import wnutils.xml as wx


class Zones_Xml:
    """A class for handling webnucleo zones."""

    def __init__(self, file):
        self.xml = wx.Xml(file)
        self.zones = self.xml.get_zone_data()

    def get_zones(self, zone_xpath=""):
        if not zone_xpath:
            return self.zones
        else:
            return self.xml.get_zone_data(zone_xpath=zone_xpath)
