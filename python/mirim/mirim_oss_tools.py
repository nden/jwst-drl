import xml.dom.minidom as minidom
import math
import mirim_tools as mirim_tools

# Note that this isn't really 'ideal' coordinates per se
def siaftov2v3(idlx,idly,angle,v2ref,v3ref):
    thisangle=-angle*3.1415926536/180.
    dv2=idlx*math.cos(thisangle)+idly*math.sin(thisangle)
    dv3=idlx*math.sin(thisangle)-idly*math.cos(thisangle)
    v2=v2ref-dv2
    v3=v3ref-dv3

    return v2,v3

def testoss():

    doc = minidom.parse('MIRI_SIAF_2017-08-28.xml')
    node = doc.documentElement
    siaf = doc.getElementsByTagName("SiafEntry")
    
    xform=mirim_tools.xytov2v3('F770W')
    
    apername = []
    xdetref = []
    ydetref = []
    v2ref = []
    v3ref = []
    # Lower-left corner
    x1 = []
    y1 = []
    # Upper-right corner
    x3 = []
    y3 = []
    angle = []
    for entry in siaf:
        AperNameObj = entry.getElementsByTagName("AperName")[0]
        nodes = AperNameObj.childNodes
        for node in nodes:
            if node.nodeType == node.TEXT_NODE:
                apername.append(node.data)

        TempObj = entry.getElementsByTagName("XDetRef")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            xdetref.append(float(node.data))

        TempObj = entry.getElementsByTagName("YDetRef")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            ydetref.append(float(node.data))
              
        TempObj = entry.getElementsByTagName("V2Ref")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            v2ref.append(float(node.data))

        TempObj = entry.getElementsByTagName("V3Ref")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            v3ref.append(float(node.data))
            
        TempObj = entry.getElementsByTagName("V3IdlYAngle")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            angle.append(float(node.data))
            
        TempObj = entry.getElementsByTagName("XIdlVert1")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            x1.append(float(node.data))

        TempObj = entry.getElementsByTagName("YIdlVert1")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            y1.append(float(node.data))

        TempObj = entry.getElementsByTagName("XIdlVert3")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            x3.append(float(node.data))

        TempObj = entry.getElementsByTagName("YIdlVert3")[0]
        nodes = TempObj.childNodes
        for node in nodes:
            y3.append(float(node.data))           
            
    nelem=len(siaf)

    # Add 1 to go to SIAF reference frame, add 0.5 to account
    # for the fact that SIAF IdlVert1 is the lower-left corner
    # of the lower-left pixel rather than the midpoint.  -0.5 for
    # IdlVert 3 since this is upper right
    

    # Define subarray locations
    totest=['MIRIM_SLITLESSPRISM']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x_sp,y_sp=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))

    totest=['MIRIM_MASKLYOT']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x_lyot,y_lyot=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))

    totest=['MIRIM_MASK1550']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x_1550,y_1550=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))

    totest=['MIRIM_MASK1140']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x_1140,y_1140=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))

    totest=['MIRIM_MASK1065']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x_1065,y_1065=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))        
        
    #############################

    # Do MIRIM_TAMRS
    totest=['MIRIM_TAMRS']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))

        print(totest[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')

    # Do MIRIM_TALRS
    totest=['MIRIM_TALRS']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))

        print(totest[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')
        
    # Do MIRIM_TASLITLESS
    totest=['MIRIM_TASLITLESSPRISM']
    ntest=len(totest)
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))
        # Adjust to slitless prism subarray
        x,y=(x-x_sp+1,y-y_sp+1)
        xu,yu=(xu-x_sp+1,yu-y_sp+1)

        print(totest[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')

    # Do Lyot TA
    totest=['MIRIM_TALYOT_UR','MIRIM_TALYOT_CUR','MIRIM_TALYOT_UL','MIRIM_TALYOT_CUL','MIRIM_TALYOT_LL','MIRIM_TALYOT_CLL','MIRIM_TALYOT_LR','MIRIM_TALYOT_CLR']
    ossnames=['1P','1S','2P','2S','3P','3S','4P','4S']
    ntest=len(totest)
    test_x = []
    test_y = []
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))
        # Adjust to Lyot subarray
        x,y=(x-x_lyot+1,y-y_lyot+1)
        xu,yu=(xu-x_lyot+1,yu-y_lyot+1)
        # Append info to vectors for IDL testing purposes
        test_x.append(x)
        test_y.append(y)
        
        print('SIAF name',totest[i])
        print('OSS region',ossnames[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')
        
    #print(test_x)
    #print(test_y)

    # Do 1550 TA
    totest=['MIRIM_TA1550_UR','MIRIM_TA1550_CUR','MIRIM_TA1550_UL','MIRIM_TA1550_CUL','MIRIM_TA1550_LL','MIRIM_TA1550_CLL','MIRIM_TA1550_LR','MIRIM_TA1550_CLR']
    ossnames=['1P','1S','2P','2S','3P','3S','4P','4S']
    ntest=len(totest)
    test_x = []
    test_y = []
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))
        # Adjust to 1550 subarray
        x,y=(x-x_1550+1,y-y_1550+1)
        xu,yu=(xu-x_1550+1,yu-y_1550+1)
        # Append info to vectors for IDL testing purposes
        test_x.append(x)
        test_y.append(y)

        print('SIAF name',totest[i])
        print('OSS region',ossnames[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')
        
    #print(test_x)
    #print(test_y)

    # Do 1140 TA
    totest=['MIRIM_TA1140_UR','MIRIM_TA1140_CUR','MIRIM_TA1140_UL','MIRIM_TA1140_CUL','MIRIM_TA1140_LL','MIRIM_TA1140_CLL','MIRIM_TA1140_LR','MIRIM_TA1140_CLR']
    ossnames=['1P','1S','2P','2S','3P','3S','4P','4S']
    ntest=len(totest)
    test_x = []
    test_y = []
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))
        # Adjust to 1140 subarray
        x,y=(x-x_1140+1,y-y_1140+1)
        xu,yu=(xu-x_1140+1,yu-y_1140+1)
        # Append info to vectors for IDL testing purposes
        test_x.append(x)
        test_y.append(y)

        print('SIAF name',totest[i])
        print('OSS region',ossnames[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')
        
    #print(test_x)
    #print(test_y)
        
    # Do 1065 TA
    totest=['MIRIM_TA1065_UR','MIRIM_TA1065_CUR','MIRIM_TA1065_UL','MIRIM_TA1065_CUL','MIRIM_TA1065_LL','MIRIM_TA1065_CLL','MIRIM_TA1065_LR','MIRIM_TA1065_CLR']
    ossnames=['1P','1S','2P','2S','3P','3S','4P','4S']
    ntest=len(totest)
    test_x = []
    test_y = []
    for i in range(0,ntest):
        # Find the index of where this aperture name is
        indx=apername.index(totest[i])
        v2,v3=siaftov2v3(x1[indx],y1[indx],angle[indx],v2ref[indx],v3ref[indx])
        x,y=xform.inverse(v2,v3)
        x,y=(int(round(x+1.0+0.5)),int(round(y+1.0+0.5)))
        v2,v3=siaftov2v3(x3[indx],y3[indx],angle[indx],v2ref[indx],v3ref[indx])
        xu,yu=xform.inverse(v2,v3)
        xu,yu=(int(round(xu+1.0-0.5)),int(round(yu+1.0-0.5)))
        # Adjust to 1065 subarray
        x,y=(x-x_1065+1,y-y_1065+1)
        xu,yu=(xu-x_1065+1,yu-y_1065+1)
        # Append info to vectors for IDL testing purposes
        test_x.append(x)
        test_y.append(y)

        print('SIAF name',totest[i])
        print('OSS region',ossnames[i])
        print('regionXsize = ',xu-x+1,';',sep='')
        print('regionYsize = ',yu-y+1,';',sep='')
        print('regionXcorner = ',x,';',sep='')
        print('regionYcorner = ',y,';',sep='')
        print('V2RefVal = ',v2ref[indx],';',sep='')
        print('V3RefVal = ',v3ref[indx],';',sep='')
        print('')
        
    #print(test_x)
    #print(test_y)


    return 0
