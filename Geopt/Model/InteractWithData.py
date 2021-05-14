import xml.etree.ElementTree as ET


def GetXML():
    treeMain = ET.parse('Data/mainblocks.xml')
    treeF = ET.parse('Data/fblock.xml')
    treerootMain = treeMain.getroot()
    treerootF = treeF.getroot()
    return treerootMain, treerootF
