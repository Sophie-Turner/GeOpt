# Functions associated with the BuildNew View
import xml.etree.ElementTree as ET

def GetXML():
    treeMain = ET.parse('../Model/mainblocks.xml')
    treeF = ET.parse('../Model/fblock.xml')
    treerootMain = treeMain.getroot()
    treerootF = treeF.getroot()
    return treerootMain, treerootF
    # how to get certain XML elements:
    #for element in treeroot:
        #print(element.tag, element.attrib)
        #print(element[2].text)
    #print(root[0][2].text)
    #print("len(treeroot):", len(treeroot))
    #print("len(treeroot[0]):", len(treeroot[0]))

