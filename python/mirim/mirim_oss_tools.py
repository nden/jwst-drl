import xml.dom.minidom as minidom

def testoss():

    doc = minidom.parse('MIRI_SIAF_2017-08-28.xml')
    node = doc.documentElement
    siaf = doc.getElementsByTagName("SiafEntry")
 
    apername = []
    xdetref = []
    ydetref = []
    v2ref = []
    v3ref = []
    x1 = []
    y1 = []
    for entry in siaf:
        AperNameObj = entry.getElementsByTagName("AperName")[0]
        nodes = AperNameObj.childNodes
        for node in nodes:
            if node.nodeType == node.TEXT_NODE:
                apername.append(node.data)

        TempObj = entry.getElementsByTagName("XDetRef")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            xdetref.append(node.data)

        TempObj = entry.getElementsByTagName("YDetRef")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            ydetref.append(node.data)
              
        TempObj = entry.getElementsByTagName("V2Ref")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            v2ref.append(node.data)

        TempObj = entry.getElementsByTagName("V3Ref")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            v3ref.append(node.data)

    nelem=len(siaf)

    totest=['MIRIM_TA1550_UL']
    ntest=len(totest)

    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        print(apername[indx],v2ref[indx])


# Remember that the vertices defined in SIAF are pixel EDGES not center, so 0.5 off
# XIDL,YIDL vertices in SIAF are offsets rotated by V3IDLYANGLE from v2,v3 system
dx=XIdlVert
dy=YIdlVert
angle=V3IdlYAngle
dv2=dx*cos(-angle)+dy*sin(-angle)
dv3=dx*sin(-angle)-dy*cos(-angle)
newv2=v2ref-dv2
newv3=v3ref-dv3
# Convert from v2v3 to x,y
# Add 1 to go to SIAF reference frame
# Add 0.5 depending for the lower-left corner


    return 0
