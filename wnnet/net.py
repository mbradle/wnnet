import wnutils.xml as wx
import wnnet.base as fb
import numpy as np

class Net(fb.Base):
    """A class for to store webnucleo networks."""

    def __init__(self, file):
        self.xml = wx.Xml(file)
        self.nuclides = self.xml.get_nuclide_data()
        self.reactions = self.xml.get_reaction_data()

    def get_nuclides(self, nuc_xpath = ""):
        if not nuc_xpath:
            return self.nuclides
        else:
            return self.xml.get_nuclide_data(nuc_xpath = nuc_xpath)

    def get_reactions(self, reac_xpath = ""):
        if not reac_xpath:
            return self.reactions
        else:
            return self.xml.get_reaction_data(reac_xpath = reac_xpath)

    def compute_Q_values(self, nuc_xpath = "", reac_xpath = ""):
        result = {}
        nuclides = self.get_nuclides(nuc_xpath = nuc_xpath)
        reactions = self.get_reactions(reac_xpath = reac_xpath)
        for r in reactions:
            tmp = self.compute_reaction_Q_value(nuclides, reactions[r])
            if tmp:
                result[r] = tmp

        return result

    def get_valid_reactions(self, nuc_xpath = "", reac_xpath = ""):
        result = {}
        nuclides = self.get_nuclides(nuc_xpath = nuc_xpath)
        reactions = self.get_reactions(reac_xpath = reac_xpath)
        for r in reactions:
            if self.is_valid_reaction(nuclides, reactions[r]):
                result[r] = reactions[r]
        return result
